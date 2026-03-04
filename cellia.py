from pathlib import Path
from typing import Optional

import scanpy as sc

def find_markers(
        adata: sc.AnnData, 
        groupby="cluster", 
        method="wilcoxon", 
        corr_method="bonferroni", 
        tie_correct=True
    ) -> sc.AnnData:
    """
    Identify marker genes for each cluster and store results in adata.uns['find_markers'].

    Parameters
    ----------
    adata : AnnData
        Single-cell AnnData object.
    groupby : str
        Key in adata.obs to group cells by.
    method : str
        Statistical test method for rank_genes_groups.
    corr_method : str
        Method for p-value correction.
    tie_correct : bool
        Whether to apply tie correction in statistical test.

    Returns
    -------
    adata : AnnData
        Updated AnnData object with marker gene results in adata.uns['find_markers'].
    """
    import scanpy as sc
    import pandas as pd

    if groupby not in adata.obs:
        raise ValueError(f"{groupby} not found in adata.obs")

    adata.obs[groupby] = (
        adata.obs[groupby]
        .astype(str)
        .astype("category")
    )

    # 1. Run statistical test
    sc.tl.rank_genes_groups(
        adata,
        use_raw=False,
        groupby=groupby,
        method=method,
        pts=True,
        corr_method=corr_method,
        tie_correct=tie_correct
    )

    # 2. Combine results into a DataFrame
    all_groups = adata.uns["rank_genes_groups"]["names"].dtype.names
    dfs = []
    for group in all_groups:
        df = sc.get.rank_genes_groups_df(adata, group=group)
        df["cluster"] = group
        dfs.append(df)

    combined_df = pd.concat(dfs, ignore_index=True)

    # 3. Rename columns for consistency
    rename_dict = {
        "pvals_adj": "p_val_adj",
        "logfoldchanges": "avg_log2FC",
        "pct_nz_group": "pct.1",
        "pct_nz_reference": "pct.2",
        "names": "gene"
    }
    combined_df = combined_df.rename(columns=rename_dict)

    # 4. Store results
    adata.uns["find_markers"] = combined_df

    return adata


def filter_markers(
        adata: sc.AnnData, 
        k=15,
        mode="db",  # "db", "subset_db"
        tissue_db=None,
        subset_db=None,
        db_path="./database/Marker_DB.csv",
        deg_mode: str = "major",  # "major" or "subset"
    ) -> sc.AnnData:
    import pandas as pd
    import scanpy as sc
    import numpy as np
    """
    Filter adata.uns['find_markers'] by thresholds and select markers per cluster.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with adata.uns['find_markers'].
    k : int
        Number of top markers to select.
    mode : str
        "db" → prioritize markers overlapping with tissue-specific DB, max k markers total
        "subset_db" → prioritize markers overlapping with tissue-specific DB, max k markers total in subset clustering
    tissue_db : str
        Tissue name for DB filtering, required if mode == "db"
    db_path : str
        Path to marker DB CSV, used if mode == "db"
    """
    if "find_markers" not in adata.uns:
        raise ValueError("adata.uns['find_markers'] not found. Run compute_find_markers first.")

    df = adata.uns["find_markers"].copy()
    marker_list = []
    
    if mode in ["db", "subset_db"]:
        #if db_path is None:
        #    raise ValueError("db_path must be provided for db or subset_db mode.")
        db_df = pd.read_csv(db_path)
        if mode == "db":
            db_filtered = db_df[db_df["tissue_type"].str.contains(tissue_db, case=False, na=False)]
        elif mode == "subset_db":
            db_filtered = db_df[
                db_df["tissue_type"].str.contains(tissue_db, case=False, na=False) &
                db_df["cell_type"].str.contains(subset_db, case=False, na=False)
            ]
        db_filtered = db_filtered.dropna(subset=["marker_genes"])

        db_genes = set()
        for gene_str in db_filtered["marker_genes"]:
            genes = [g.strip().upper() for g in gene_str.split(",") if g.strip() and g.strip() != ","]
            db_genes.update(genes)
    
    deg_mode = deg_mode.lower()
    if deg_mode not in ["major", "subset"]:
        raise ValueError("deg_mode must be either 'major' or 'subset'.")

    for cluster in df["cluster"].unique():
        sub_df = df[df["cluster"] == cluster].copy()

        sub_df["pct_diff"] = sub_df["pct.1"] - sub_df["pct.2"]

        if deg_mode == "major":
            sub_df = sub_df[
                (sub_df["p_val_adj"] < 0.05) &
                (sub_df['pct.1'] > 0.25) &
                (sub_df["pct.2"] < 0.35) &
                (sub_df["avg_log2FC"] > 0.3)
            ]
        elif deg_mode == "subset":
            sub_df = sub_df[
                (sub_df["p_val_adj"] < 0.1) &
                (sub_df['pct.1'] > 0.25) &
                (sub_df["pct.2"] < 0.7) &
                (sub_df["avg_log2FC"] > 0.3)
            ]

        sub_df = sub_df.sort_values(by="avg_log2FC", ascending=False)

        if mode in ["db", "subset_db"]:
            gene_list = sub_df["gene"].tolist()
            marker_list_clean = [gene.strip().upper() for gene in gene_list if gene.strip()]

            overlapping_genes = [gene for gene in marker_list_clean if gene in db_genes]
            non_overlapping_genes = [gene for gene in marker_list_clean if gene not in overlapping_genes]
            final_genes = overlapping_genes + non_overlapping_genes
            final_genes = final_genes[:k]

            # Rebuild sub_df with only these genes, keeping the original metrics
            sub_df = sub_df[sub_df["gene"].str.upper().isin(final_genes)]

        marker_list.append(sub_df[[
            "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene"
        ]])


    import pandas as pd
    # Combine and store
    final_df = pd.concat(marker_list, ignore_index=True)
    if mode == "db":
        key_name = "marker_list"
    elif mode == "subset_db":
        key_name = "marker_list_subset"
        
    adata.uns[key_name] = final_df

    return adata

def gpt_anno(   
        adata: sc.AnnData, 
        tissue_type: str, 
        api_key: str, 
        model="gpt-4.1-2025-04-14",
        mode="major", # major or subset
        parent_celltype=None,
        db_path="./database/Marker_DB.csv"
    ) -> sc.AnnData:
    """
    GPT-based cell type annotation from marker genes.
    
    Parameters
    ----------
    adata : AnnData
        Single-cell AnnData object with adata.uns['marker_list'] and adata.obs['cluster'].
    tissue_type : str
        Tissue name for context (e.g., "human PBMC").
    api_key : str
        OpenAI API key.
    model : str
        GPT model name (default: gpt-4o).
    mode : str
        "db" → top k marker genes from adata.uns['marker_list']
        "subset_db" → top k marker genes from adata.uns['marker_list_subset']
    db_path : str
        Path to marker DB CSV file (used if mode="db").
    """
    from openai import OpenAI
    from tqdm import tqdm
    import pandas as pd
    import json
    import os
    import scanpy as sc
    import matplotlib.pyplot as plt

    if mode not in ["major", "subset"]:
        raise ValueError("marker_mode must be one of: 'major', 'subset'")

    if mode == "major":
        key_name = "marker_list"
    elif mode == "subset":
        key_name = "marker_list_subset"

    # --- Input checks ---
    if mode == "major":
        if key_name not in adata.uns:
            raise ValueError(f"adata.uns['{key_name}'] not found. Please run marker filtering first.")
    elif mode == "subset":
        if key_name not in adata.uns:
            raise ValueError(f"adata.uns['{key_name}'] not found. Please run marker filtering first.")
    
    if "cluster" not in adata.obs:
        raise ValueError("adata.obs must contain a 'cluster' column.")

    
    if mode == "major":
        marker_df = adata.uns[key_name]
    elif mode == "subset":
        marker_df = adata.uns[key_name]

    client = OpenAI(api_key=api_key)

    # --- Init output containers ---
    cluster_annotations = {}
    marker_explanations = {}
    top_dict = {}
    unique_clusters = marker_df["cluster"].unique()

    output_dir = "./cellia_output"
    os.makedirs(output_dir, exist_ok=True)

    # --- Loop over clusters ---
    for cluster in tqdm(unique_clusters, desc="LLM annotation per cluster"):
        marker_list = marker_df[marker_df["cluster"] == cluster]["gene"].tolist()

        marker_genes_str = ", ".join(marker_list)
        top_dict[str(cluster)] = marker_genes_str  

        if mode == "major":
            # --- major cell type annotation ---
            messages = [
                {
                    "role": "system",
                    "content": "You are a helpful assistant that identifies cell types based on marker genes and explains the reasoning for cell type annotation based on marker genes."
                    },
                {
                    "role": "user",
                    "content": (
                        f"Using single-cell RNA-seq data, we performed clustering and filtered the marker genes in each cluster. "
                        f"We identified the filtered marker genes in a particular cluster as follows:\n{marker_genes_str}\n\n"
                        f"The tissue sample is from {tissue_type}.\n"
                        "Given this set of marker genes identified in a cluster, which cell type best matches this expression profile? "
                        "Only provide the cell type name.\n"
                        "Format the response as a JSON object with the following key:\n"
                        '{"cell_type": "Predicted cell type name"}\n\n'
                        "Also, from the marker gene list, select the most important marker genes that strongly support this cell type, and briefly explain why each one is relevant.\n"
                        "Return the explanation as a JSON object where each key is a gene name and the value is a short reason.\n\n"
                        "Final format should be a single JSON object with two keys:\n"
                        "- `cell_type`: predicted cell type (string)\n"
                        "- `marker_explanations`: dictionary of gene: explanation\n\n"
                        "**Return only a single-line JSON object, no extra explanation or code block.**\n"
                        "Format:\n"
                        '{"cell_type": "PredictedCellType", "marker_explanations": {"GENE1": "Reason1", "GENE2": "Reason2"}}'
                    )
                }
            ]
        elif mode == "subset":
            messages = [
                {
                    "role": "system",
                    "content": "You are a domain expert in cell biology and single-cell transcriptomics who identifies subset annotation based on marker genes and explains the reasoning for cell type annotation based on marker genes."
                    },
                {
                    "role": "user",
                    "content": (
                        "Given the tissue type, the parent cell type, and a list of marker genes for a subcluster, identify the most likely subset cell type (subtype/activation state/developmental stage/lineage breanch) for this subcluster.\n"
                        f"Marker genes (ONLY use these; ignore others): {marker_genes_str}\n"
                        f"Tissue type: {tissue_type}\n"
                        f"Parent cell type: {parent_celltype}\n"
                        "   - These markers are characteristic of the parent cell type and can be used as an identity anchor when choosing the subset label."
                        "Instructions;\n"
                        "1. From the provided marker genes, identify the most specific and relevant markers that define a known subset of the parent cell type in this tissue context.\n"
                        "2. Based on these markers, determine the most likely subset label strictly within the parent cell type and tissue context.\n"
                        "3. Provide brief explanations for why each selected marker supports this subset annotation (gene -> reason mapping).\n"
                        "4. If evidence is borderline between two subsets, choose the best single label but reflect ambiguity in the score reason.\n\n"
                        "Final format should be a single JSON object with two keys:\n"
                        "- `cell_type`: predicted cell type (string)\n"
                        "- `marker_explanations`: dictionary of gene: explanation\n\n"
                        "**Return only a single-line JSON object, no extra explanation or code block.**\n"
                        "Format:\n"
                        '{"cell_type": "PredictedCellType", "marker_explanations": {"GENE1": "Reason1", "GENE2": "Reason2"}}'
                    )
                }
            ]


        # --- GPT Call ---
        try:
            completion = client.chat.completions.create(
                model=model,
                messages=messages,
                temperature=0.2, 
                max_tokens=2048
            )
            response_content = completion.choices[0].message.content

            cell_type_value = json.loads(response_content).get("cell_type", "Unknown")
            print(f"Raw response for cluster {cluster}: {cell_type_value}")

            try:
                response_data = json.loads(response_content)
                cell_type_name = response_data.get("cell_type", "Unknown")
                explanations = response_data.get("marker_explanations", {})
            except json.JSONDecodeError:
                print(f"Cluster {cluster}: Failed to parse JSON, using raw text instead")
                cell_type_name = response_content.strip()
                explanations = {"error": response_content.strip()}

        except Exception as e:
            print(f"Failed to annotate cluster {cluster}: {e}")
            cell_type_name = "Unknown"
            explanations = {"error": str(e)}

        # --- Save per cluster ---
        cluster_annotations[str(cluster)] = cell_type_name
        marker_explanations[str(cluster)] = explanations

    # --- Save results to adata ---
    df_annotation = pd.DataFrame({
        "cluster": list(cluster_annotations.keys()),
        "LLM_annotation": list(cluster_annotations.values()),
        "markers": [top_dict.get(c, []) for c in cluster_annotations.keys()],
    })

    adata.uns["GPT_annotation_" + mode] = df_annotation

    # --- Save JSON and CSV ---
    if mode=="subset":
        annotations_payload = {
            str(cluster): {
                "cell_type": cluster_annotations[str(cluster)],
                "markers": top_dict.get(str(cluster), []),
            }
            for cluster in cluster_annotations
        }
    else:
        annotations_payload = cluster_annotations

    output_path = os.path.join(output_dir, f"gpt_annotations_{mode}.json")
    with open(output_path, "w") as f:
        json.dump(annotations_payload, f, indent=4)

    cluster_explanations = {}
    for cluster in cluster_annotations:
        payload = {
            "cell_type": cluster_annotations[str(cluster)],
            "marker_explanations": marker_explanations.get(str(cluster), {}),
        }
        if mode == "subset":
            payload["markers"] = top_dict.get(str(cluster), "")
        cluster_explanations[str(cluster)] = payload

    with open(os.path.join(output_dir, f"gpt_explanations_{mode}.json"), "w") as f:
        json.dump(cluster_explanations, f, indent=4)


    explanation_rows = []
    for cluster, entry in cluster_explanations.items():
        cell_type = entry.get("cell_type", "Unknown")

        for gene, reason in entry.get("marker_explanations", {}).items():
            row = {
                "cluster": cluster,
                "cell_type": cell_type,
                "gene": gene,
                "explanation": reason,
            }
            explanation_rows.append(row)

    pd.DataFrame(explanation_rows).to_csv(
        os.path.join(output_dir, f"gpt_explanations_{mode}.csv"), 
        index=False
    )

    print(f"GPT-based annotation complete (mode={mode}). Results saved to: {output_dir}")
    return adata


def gemini_anno(
        adata: sc.AnnData, 
        tissue_type: str, 
        api_key: str, 
        model_name="gemini-2.5-flash-lite",
        mode="major", # major or subset
        parent_celltype=None,
    ) -> sc.AnnData:
    """
    Gemini-based cell type annotation from marker genes.

    Parameters
    ----------
    adata : AnnData
        Single-cell AnnData object with adata.uns['marker_list'] or adata.uns['marker_list_subset'],
        and adata.obs['cluster'] present.
    tissue_type : str
        Tissue name for context (e.g., "human PBMC").
    api_key : str
        Anthropic API key.
    model : str
        Gemini model name.
    mode : str
        "major" → major cell type annotation
        "subset" → subset/activation/state under a parent cell type (requires parent_celltype)
    """

    #import google.generativeai as genai
    from google import genai
    from tqdm import tqdm
    import pandas as pd
    import json
    import os
    import scanpy as sc
    import matplotlib.pyplot as plt
    import re
    import time
    from random import uniform

    if mode not in ["major", "subset"]:
        raise ValueError("marker_mode must be one of: 'major', 'subset'")

    if mode == "major":
        key_name = "marker_list"
    elif mode == "subset":
        key_name = "marker_list_subset"

    # --- Input checks ---
    if mode == "major":
        if key_name not in adata.uns:
            raise ValueError(f"adata.uns['{key_name}'] not found. Please run marker filtering first.")
    elif mode == "subset":
        if key_name not in adata.uns:
            raise ValueError(f"adata.uns['{key_name}'] not found. Please run marker filtering first.")

    if "cluster" not in adata.obs:
        raise ValueError("adata.obs must contain a 'cluster' column.")

    marker_df = adata.uns[key_name]

    cluster_annotations = {}
    marker_explanations = {}
    top_dict = {}
    neg_top_dic = {}

    unique_clusters = marker_df["cluster"].unique()

    output_dir = "./cellia_output"
    os.makedirs(output_dir, exist_ok=True)

    for cluster in tqdm(unique_clusters, desc="LLM annotation per cluster"):
        marker_list = marker_df[marker_df["cluster"] == cluster]["gene"].tolist()
        marker_genes_str = ", ".join(marker_list)
        top_dict[str(cluster)] = marker_genes_str

        if mode == "major":
            prompt = (
                "You are a helpful assistant that identifies cell types based on marker genes and explains the reasoning for cell type annotation based on marker genes."
                f"Using single-cell RNA-seq data, we performed clustering and filtered the marker genes in each cluster. "
                f"We identified the filtered marker genes in a particular cluster as follows:\n{marker_genes_str}\n\n"
                f"The tissue sample is from {tissue_type}.\n"
                "Given this set of marker genes identified in a cluster, which cell type best matches this expression profile? "
                "Only provide the cell type name.\n"
                "Format the response as a JSON object with the following key:\n"
                '{"cell_type": "Predicted cell type name"}\n\n'
                "Also, from the marker gene list, select the most important marker genes that strongly support this cell type, and briefly explain why each one is relevant.\n"
                "Return the explanation as a JSON object where each key is a gene name and the value is a short reason.\n\n"
                "Justify the score with a short reason referencing the relevant markers.\n\n"
                "Final format should be a single JSON object with two keys:\n"
                "- `cell_type`: predicted cell type (string)\n"
                "- `marker_explanations`: dictionary of gene: explanation\n\n"
                "**Return only a single-line JSON object, no extra explanation or code block.**\n"
                "Format:\n"
                '{"cell_type": "PredictedCellType", "marker_explanations": {"GENE1": "Reason1", "GENE2": "Reason2"}}'
            )

        elif mode == "subset":
            prompt = (
                "You are a domain expert in cell biology and single-cell transcriptomics who identifies subset annotation based on marker genes and explains the reasoning for cell type annotation based on marker genes."
                "Given the tissue type, the parent cell type, and a list of marker genes for a subcluster, identify the most likely subset cell type (subtype/activation state/developmental stage/lineage breanch) for this subcluster.\n"
                f"Tissue type: {tissue_type}\n"
                f"Parent cell type: {parent_celltype}\n"
                "   - These markers are characteristic of the parent cell type and can be used as an identity anchor when choosing the subset label."
                "Instructions:\n"
                "1. From the provided marker genes, identify the most specific and relevant markers that define a known subset of the parent cell type in this tissue context.\n"
                "2. Based on these markers, determine the most likely subset label strictly within the parent cell type and tissue context.\n"
                "3. Provide brief explanations for why each selected marker supports this subset annotation (gene -> reason mapping).\n\n"
                "Final format should be a single JSON object with two keys:\n"
                "- `cell_type`: predicted cell type (string)\n"
                "- `marker_explanations`: dictionary of gene: explanation\n\n"
                "**Return only a single-line JSON object, no extra explanation or code block.**\n"
                "Format:\n"
                '{"cell_type": "PredictedCellType", "marker_explanations": {"GENE1": "Reason1", "GENE2": "Reason2"}}'
            )

        try:
            time.sleep(20)

            client=genai.Client(api_key=api_key)

            response=client.models.generate_content(
                model=model_name,
                contents=prompt,
            )
            response_content = response.text.strip()

            response_content = re.sub(r"^```(?:json)?|```$", "", response_content.strip(), flags=re.IGNORECASE)

            cell_type_value = json.loads(response_content).get("cell_type", "Unknown")
            print(f"Raw response for cluster {cluster}: {cell_type_value}")

            try:
                response_data = json.loads(response_content)
                cell_type_name = response_data.get("cell_type", "Unknown")
                explanations = response_data.get("marker_explanations", {})
            except json.JSONDecodeError:
                print(f"Cluster {cluster}: Failed to parse JSON, using raw text instead")
                cell_type_name = response_content.strip()
                explanations = {"error": response_content.strip()}

        except Exception as e:
            print(f"Failed to annotate cluster {cluster}: {e}")
            cell_type_name = "Unknown"
            explanations = {"error": str(e)}

        # Save
        cluster_annotations[str(cluster)] = cell_type_name
        marker_explanations[str(cluster)] = explanations

    # Save annotations to adata
    df_annotation = pd.DataFrame({
        "cluster": list(cluster_annotations.keys()),
        "LLM_annotation": list(cluster_annotations.values()),
        "markers": [top_dict.get(c, []) for c in cluster_annotations.keys()],
    })

    adata.uns["Gemini_annotation_" + mode] = df_annotation

    # -----------------------------
    # Save JSON (annotation)
    # -----------------------------
    annotations_payload = cluster_annotations

    output_path = os.path.join(output_dir, f"gemini_annotations_{mode}.json")
    with open(output_path, "w") as f:
        json.dump(annotations_payload, f, indent=4)

    # -----------------------------
    # Save JSON (explanations) + CSV
    # -----------------------------
    cluster_explanations = {}
    for cluster in cluster_annotations:
        payload = {
            "cell_type": cluster_annotations[str(cluster)],
            "marker_explanations": marker_explanations.get(str(cluster), {}),
        }
        cluster_explanations[str(cluster)] = payload

    with open(os.path.join(output_dir, f"gemini_explanations_{mode}.json"), "w") as f:
        json.dump(cluster_explanations, f, indent=4)

    # Save explanations as CSV
    explanation_rows = []
    for cluster, entry in cluster_explanations.items():
        cell_type = entry.get("cell_type", "Unknown")

        for gene, reason in entry.get("marker_explanations", {}).items():
            row = {
                "cluster": cluster,
                "cell_type": cell_type,
                "gene": gene,
                "explanation": reason,
                "evidence_score": entry.get("evidence_score", "unknown"),
                "evidence_reason": entry.get("evidence_reason", ""),
            }
            explanation_rows.append(row)

    pd.DataFrame(explanation_rows).to_csv(os.path.join(output_dir, f"gemini_explanations_{mode}.csv"), index=False)

    print(f"Gemini-based annotation complete (mode={mode}). Results saved to: {output_dir}")
    return adata

def claude_anoo(
        adata: sc.AnnData,
        tissue_type: str,
        api_key: str,
        model="claude-sonnet-4-5",
        mode="major",       # "major" or "subset"
        parent_celltype=None,
    ) -> sc.AnnData:
    """
    Claude-based cell type annotation from marker genes.

    Parameters
    ----------
    adata : AnnData
        Single-cell AnnData object with adata.uns['marker_list'] or adata.uns['marker_list_subset'],
        and adata.obs['cluster'] present.
    tissue_type : str
        Tissue name for context (e.g., "human PBMC").
    api_key : str
        Anthropic API key.
    model : str
        Claude model name.
    mode : str
        "major" → major cell type annotation
        "subset" → subset/activation/state under a parent cell type (requires parent_celltype)
    db_path : str
        (kept for API compatibility; not used directly in this function)
    """
    import json
    import os
    import pandas as pd
    from tqdm import tqdm
    import re

    # Anthropic SDK
    try:
        import anthropic
    except ImportError as e:
        raise ImportError(
            "anthropic package is required. Install with `pip install anthropic`"
        ) from e

    # --- Input checks ---
    if mode not in ["major", "subset"]:
        raise ValueError("marker_mode must be one of: 'major', 'subset'")

    key_name = "marker_list" if mode == "major" else "marker_list_subset"

    if key_name not in adata.uns:
        raise ValueError(f"adata.uns['{key_name}'] not found. Please run marker filtering first.")

    if "cluster" not in adata.obs:
        raise ValueError("adata.obs must contain a 'cluster' column.")

    if mode == "subset" and not parent_celltype:
        raise ValueError("parent_celltype is required when mode='subset'")

    marker_df = adata.uns[key_name]

    # --- Anthropic client ---
    client = anthropic.Anthropic(api_key=api_key)

    # --- Init output containers ---
    cluster_annotations = {}
    marker_explanations = {}
    top_dict = {}
    unique_clusters = marker_df["cluster"].unique()

    output_dir = "./cellia_output"
    os.makedirs(output_dir, exist_ok=True)

    # --- Prompt templates ---
    if mode == "major":
        system_prompt = (
            "You are a helpful assistant that identifies cell types based on marker genes "
            "and explains the reasoning for cell type annotation based on marker genes."
            "Return ONLY raw JSON. Do NOT wrap in backticks or a code block."
        )
        def user_prompt(marker_genes_str: str) -> str:
            return (
                "Using single-cell RNA-seq data, we performed clustering and filtered the marker genes in each cluster. "
                f"We identified the filtered marker genes in a particular cluster as follows:\n{marker_genes_str}\n\n"
                f"The tissue sample is from {tissue_type}.\n"
                "Given this set of marker genes identified in a cluster, which cell type best matches this expression profile? "
                "Only provide the cell type name.\n"
                "Format the response as a JSON object with the following key:\n"
                '{"cell_type": "Predicted cell type name"}\n\n'
                "Also, from the marker gene list, select the most important marker genes that strongly support this cell type, "
                "and briefly explain why each one is relevant.\n"
                "Return the explanation as a JSON object where each key is a gene name and the value is a short reason.\n\n"
                "Final format should be a single JSON object with three keys:\n"
                "- `cell_type`: predicted cell type (string)\n"
                "- `marker_explanations`: dictionary of gene: explanation\n"
                "**Return only a single-line JSON object, no extra explanation or code block.**\n"
                'Format:\n'
                '{"cell_type": "PredictedCellType", "marker_explanations": {"GENE1": "Reason1", "GENE2": "Reason2"}}'
            )
    else:  # subset
        system_prompt = (
            "You are a domain expert in cell biology and single-cell transcriptomics who identifies subset annotations "
            "based on marker genes and explains the reasoning."
            "Return ONLY raw JSON. Do NOT wrap in backticks or a code block."
        )
        def user_prompt(marker_genes_str: str) -> str:
            return (
                "Given the tissue type, the parent cell type, and a list of marker genes for a subcluster, identify the most likely subset "
                "(subtype/activation state/developmental stage/lineage branch) for this subcluster.\n"
                f"Marker genes (ONLY use these; ignore others): {marker_genes_str}\n"
                f"Tissue type: {tissue_type}\n"
                f"Parent cell type: {parent_celltype}\n"
                "Instructions:\n"
                "1. From the provided marker genes, identify the most specific and relevant markers that define a known subset of the parent cell type in this tissue context.\n"
                "2. Based on these markers, determine the most likely subset label strictly within the parent cell type and tissue context.\n"
                "3. Provide brief explanations for why each selected marker supports this subset annotation (gene -> reason mapping).\n\n"
                "Final format should be a single JSON object with three keys:\n"
                "- `cell_type`: predicted subset (string)\n"
                "- `marker_explanations`: dictionary of gene: explanation\n"
                "**Return only a single-line JSON object, no extra explanation or code block.**\n"
                'Format:\n'
                '{"cell_type": "PredictedSubset", "marker_explanations": {"GENE1": "Reason1", "GENE2": "Reason2"}}'
            )

    def _anthropic_text_from_message(msg):
        """Safely flatten Anthropic message.content (list of blocks) to text."""
        # msg.content is a list of content blocks; pick text blocks and join.
        try:
            blocks = getattr(msg, "content", None) or []
            texts = []
            for b in blocks:
                if getattr(b, "type", None) == "text":
                    texts.append(getattr(b, "text", ""))
                elif isinstance(b, dict) and b.get("type") == "text":
                    texts.append(b.get("text", ""))
            return " ".join(t for t in texts if t).strip()
        except Exception:
            return str(msg)

    # --- Loop over clusters ---
    for cluster in tqdm(unique_clusters, desc="LLM annotation per cluster"):
        marker_list = marker_df[marker_df["cluster"] == cluster]["gene"].tolist()
        marker_genes_str = ", ".join(marker_list)
        top_dict[str(cluster)] = marker_genes_str

        content = user_prompt(marker_genes_str)

        messages = [{"role": "user", "content": content}]


        try:
            message = client.messages.create(
                model=model,
                max_tokens=2048,
                temperature=0.2,
                system=system_prompt,
                messages=messages,
                stop_sequences=["```"],
            )
            response_content = _anthropic_text_from_message(message)

            cell_type_value = json.loads(response_content).get("cell_type", "Unknown")
            print(f"Raw response for cluster {cluster}: {cell_type_value}")

            # Try to parse JSON (expects single-line JSON per prompt)
            try:
                response_data = json.loads(response_content)
                cell_type_name = response_data.get("cell_type", "Unknown")
                explanations = response_data.get("marker_explanations", {}) or {}
            except json.JSONDecodeError:
                # Fallback: store raw text for debugging
                cell_type_name = response_content.strip() or "Unknown"
                explanations = {"_parse_error": response_content.strip()}

        except Exception as e:
            print(f"Failed to annotate cluster {cluster}: {e}")
            cell_type_name = "Unknown"
            explanations = {"_error": str(e)}

        # Save per cluster
        cid = str(cluster)
        cluster_annotations[cid] = cell_type_name
        marker_explanations[cid] = explanations

    # --- Save results to adata ---
    df_annotation = pd.DataFrame({
        "cluster": list(cluster_annotations.keys()),
        "LLM_annotation": list(cluster_annotations.values()),
        "markers": [top_dict.get(c, []) for c in cluster_annotations.keys()],
    })

    adata.uns["Claude_annotation_" + mode] = df_annotation  # keep key name for backward-compat

    # --- Save JSON and CSV ---
    # annotations json: subset일 때 negative_markers 포함해서 저장
    annotations_payload = cluster_annotations

    output_path = os.path.join(output_dir, f"claude_annotations_{mode}.json")
    with open(output_path, "w") as f:
        json.dump(annotations_payload, f, indent=4)

    # explanations json: subset일 때 negative_markers 포함
    cluster_explanations = {}
    for cid in cluster_annotations.keys():
        payload = {
            "cell_type": cluster_annotations[cid],
            "marker_explanations": marker_explanations.get(cid, {}),
        }
        cluster_explanations[cid] = payload
        
    with open(os.path.join(output_dir, f"claude_explanations_{mode}.json"), "w") as f:
        json.dump(cluster_explanations, f, indent=4)

    # explanations csv: subset일 때 negative_markers 컬럼 포함
    explanation_rows = []
    for cid, entry in cluster_explanations.items():
        cell_type = entry.get("cell_type", "Unknown")
        for gene, reason in (entry.get("marker_explanations") or {}).items():
            row = {
                "cluster": cid,
                "cell_type": cell_type,
                "gene": gene,
                "explanation": reason,
            }
            explanation_rows.append(row)

    import pandas as pd
    pd.DataFrame(explanation_rows).to_csv(
        os.path.join(output_dir, f"claude_explanations_{mode}.csv"),
        index=False
    )

    print(f"Claude-based annotation complete (mode={mode}). Results saved to: {output_dir}")
    return adata


def cellia_run(
    adata: sc.AnnData,
    tissue_db: str,
    tissue_type: str,
    api_key: str,
    n_top_markers: int = 15,

    groupby: str = "cluster",

    deg_mode: str = "major",          # "major" or "subset" 
    db_mode: str = "db",              # "db" or "subset_db" 
    db_path: str = "./database/Marker_DB.csv",

    llm_provider: str = "gpt",        # "gpt" or "gemini" or "claude"
    llm_mode: str = "major",          # "major" or "subset" 
    model: Optional[str] = None,      
    parent_celltype: Optional[str] = None,
    subset_db: Optional[str] = None,  # required if db_mode == "subset_db"
    output_dir: str = "./cellia_output",
) -> sc.AnnData:

    if db_mode not in ["db", "subset_db"]:
        raise ValueError("db_mode must be 'db' or 'subset_db'")

    if llm_provider not in ["gpt", "gemini", "claude"]:
        raise ValueError("llm_provider must be one of: 'gpt', 'gemini', 'claude'")

    if llm_mode not in ["db", "subset"]:
        raise ValueError("llm_mode must be one of: 'db', 'subset'")
    
    if db_mode == "subset_db" and not subset_db:
        raise ValueError("subset_db must be provided when db_mode='subset_db'")
    
    if llm_mode == "subset" and not parent_celltype:
        raise ValueError("parent_celltype must be provided when llm_mode='subset'")

    adata = find_markers(
        adata=adata,
        groupby=groupby
    )

    adata = filter_markers(
        adata=adata,
        k=n_top_markers,
        mode=db_mode,
        tissue_db=tissue_db,
        subset_db=subset_db,
        db_path=db_path,
        deg_mode=deg_mode,
    )


    if llm_provider == "gpt":
        _model = model or "gpt-4.1-2025-04-14"
        adata = gpt_anno(
            adata=adata,
            tissue_type=tissue_type,
            api_key=api_key,
            model=_model,
            mode=llm_mode,
            parent_celltype=parent_celltype,
            db_path=db_path,
        )

    elif llm_provider == "gemini":
        _model = model or "gemini-2.5-flash"
        adata = gemini_anno(
            adata=adata,
            tissue_type=tissue_type,
            api_key=api_key,
            model_name=_model,
            mode=llm_mode,
            parent_celltype=parent_celltype,
        )

    elif llm_provider == "claude":
        _model = model or "claude-sonnet-4-5"
        adata = claude_anoo(
            adata=adata,
            tissue_type=tissue_type,
            api_key=api_key,
            model=_model,
            mode=llm_mode,
            parent_celltype=parent_celltype,
        )

    return adata

__all__ = ["cellia_run", "find_markers", "filter_markers", "gpt_anno"]

