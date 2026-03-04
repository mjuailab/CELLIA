from dash import Dash, dcc, html, Input, Output, State, no_update
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
import numpy as np
from typing import Optional, Dict, Any
import os
import json



# GeneCards
def gene_url(gene: str) -> str:
    gene = str(gene).strip()
    return f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene}"


cap_template = go.layout.Template(
    layout=go.Layout(
        font=dict(family="Inter, 'Helvetica Neue', Arial, sans-serif", size=14, color="#0b1220"),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="#FFFFFF",
        hoverlabel=dict(
            bgcolor="#FFFFFF",
            font_size=13,
            font_family="Inter, Arial",
            bordercolor="rgba(15,23,42,0.12)",
            font=dict(color="#0b1220"),
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1,
            bgcolor="rgba(255,255,255,0.85)",
            bordercolor="rgba(15,23,42,0.10)",
            borderwidth=1,
        ),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, title=None),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, title=None),
        coloraxis_colorbar=dict(outlinewidth=0, ticks="outside", ticklen=4),
        margin=dict(l=40, r=40, t=40, b=40),
    )
)
pio.templates["cap_light_pro"] = cap_template
pio.templates.default = "cap_light_pro"


PASTEL = (
    px.colors.qualitative.Pastel
    + px.colors.qualitative.Set3
    + px.colors.qualitative.D3
    + px.colors.qualitative.T10
)


def _discrete_map(series):
    cats = list(pd.Series(series).astype(str).unique())
    return {c: PASTEL[i % len(PASTEL)] for i, c in enumerate(sorted(cats))}


def _parse_marker_expl(value: str) -> dict:
    if not isinstance(value, str) or not value.strip():
        return {}
    try:
        obj = json.loads(value)
        if isinstance(obj, dict):
            return {str(k): str(v) for k, v in obj.items()}
    except Exception:
        pass
    expl = {}
    for p in value.replace(";", ",").split(","):
        if ":" in p:
            g, r = p.split(":", 1)
            if g.strip():
                expl[g.strip()] = r.strip()
    return expl



def build_marker_info_from_uns(
    adata,
    use_uns_markers: str = "marker_list",
    use_uns_llm_annotation: str = "GPT_annotation_db",
    cluster_key="cluster",
    num_top_k: int = 15,
    rationale_json: Optional[Dict[str, Any]] = None,
):
    df_m = pd.DataFrame(adata.uns[use_uns_markers]).copy()

    df_g = pd.DataFrame(
        adata.uns.get(
            use_uns_llm_annotation,
            pd.DataFrame(columns=["cluster", "LLM_annotation", "markers", "marker_explanations"]),
        )
    ).copy()

    if "gene" in df_m.columns:
        df_m["gene"] = df_m["gene"].astype(str).str.strip()

    df_m["cluster"] = df_m["cluster"].astype(str)
    if not df_g.empty and "cluster" in df_g.columns:
        df_g["cluster"] = df_g["cluster"].astype(str)

    def _split_markers(s):
        if pd.isna(s):
            return []
        return [x.strip() for x in str(s).replace(";", ",").split(",") if x.strip()]

    if "markers" in df_g.columns:
        df_g["markers_list"] = df_g["markers"].map(_split_markers)
    else:
        df_g["markers_list"] = [[] for _ in range(len(df_g))]

    db_by_cluster = (
        df_m.sort_values(["cluster", "avg_log2FC"], ascending=[True, False])
        .groupby("cluster")
        .apply(
            lambda sub: {
                "db_genes": sub["gene"].tolist(),
                "top15_db": sub["gene"].head(num_top_k).tolist(),
                "stats": sub[["gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"]].to_dict("records"),
            }
        )
        .to_dict()
    )

    def process_llm_group(sub):
        first_row = sub.iloc[0]
        explanations = {}
        if "marker_explanations" in first_row and pd.notna(first_row["marker_explanations"]):
            explanations = _parse_marker_expl(first_row["marker_explanations"])
        return pd.Series(
            {
                "cell_type": first_row.get("LLM_annotation", "Unknown"),
                "llm_markers": first_row.get("markers_list", []),
                "marker_explanations": explanations,
            }
        )

    llm_by_cluster = {}
    if not df_g.empty and "cluster" in df_g.columns:
        llm_by_cluster = df_g.groupby("cluster").apply(process_llm_group).to_dict("index")

    rationale_map = {str(k): v for k, v in rationale_json.items()} if rationale_json else {}

    clusters = sorted(
        set(df_m["cluster"])
        | (set(df_g["cluster"]) if ("cluster" in df_g.columns and not df_g.empty) else set())
        | set(adata.obs[cluster_key].astype(str))
    )

    marker_info = {}
    for cid in clusters:
        llm = llm_by_cluster.get(cid, {})
        db = db_by_cluster.get(cid, {})
        rj = rationale_map.get(cid, {})

        cell_type = rj.get("cell_type", llm.get("cell_type", "Unknown"))

        expl = {}
        llm_explained_genes = set()

        if rj.get("marker_explanations"):
            expl.update(rj["marker_explanations"])
            llm_explained_genes.update(rj["marker_explanations"].keys())
        elif llm.get("marker_explanations"):
            expl.update(llm["marker_explanations"])
            llm_explained_genes.update(llm["marker_explanations"].keys())
            for gene in db.get("top15_db", []):
                expl.setdefault(gene, "from DB (top 15 by avg_log2FC)")
        else:
            llm_genes = llm.get("llm_markers", [])
            llm_explained_genes.update(llm_genes)
            for gene in llm_genes:
                expl[gene] = "LLM-selected marker"
            for gene in db.get("top15_db", []):
                expl.setdefault(gene, "from DB (top 15 by avg_log2FC)")

        marker_info[cid] = {
            "cell_type": cell_type,
            "marker_explanations": expl,
            "db_stats": db.get("stats", []),
            "llm_markers": llm.get("llm_markers", []),
            "llm_explained_genes": list(llm_explained_genes),
            "db_genes": db.get("db_genes", []),
            "top15_db": db.get("top15_db", []),
        }

    return marker_info


# UMAP 
def _pick_umap_key(adata, preferred="X_umap") -> str:
    if preferred in adata.obsm:
        return preferred
    for k in adata.obsm.keys():
        if "umap" in str(k).lower():
            return k
    raise ValueError("No UMAP embedding found in adata.obsm (expected 'X_umap' or a key containing 'umap').")


def make_umap_df(adata, marker_info, cluster_key="cluster"):
    umap_key = _pick_umap_key(adata, preferred="X_umap")
    XY = np.asarray(adata.obsm[umap_key])
    if XY.ndim != 2 or XY.shape[1] < 2:
        raise ValueError(f"UMAP embedding '{umap_key}' must be 2D with at least 2 columns, got shape={XY.shape}")

    cell_ids = pd.Index(adata.obs_names).astype(str)

    df = pd.DataFrame(
        {
            "UMAP1": XY[:, 0],
            "UMAP2": XY[:, 1],
            "_cell_id": cell_ids,
            "barcodes": cell_ids,  # alias for compatibility
        }
    )

    df["cluster"] = adata.obs[cluster_key].astype(str).values
    df["cell_type"] = df["cluster"].map(lambda c: marker_info.get(c, {}).get("cell_type", "Unknown"))
    df["annotation"] = df["cluster"] + " (" + df["cell_type"] + ")"

    def get_markers_str(cluster):
        genes = list(marker_info.get(cluster, {}).get("marker_explanations", {}).keys())[:10]
        return "<br>".join(genes) if genes else "N/A"

    df["top_markers_str"] = df["cluster"].apply(get_markers_str)
    return df


def precompute_cluster_gene_stats_dense(adata, marker_info, cluster_key="cluster"):
    genes_set = {g for v in marker_info.values() for g in v.get("llm_markers", []) + v.get("top15_db", [])}
    genes = [g for g in genes_set if g in adata.var_names]
    if not genes:
        return {"genes": [], "avg": [], "pct": []}

    Xsub = adata[:, genes].X
    X = Xsub.toarray() if hasattr(Xsub, "toarray") else np.asarray(Xsub)

    df_expr = pd.DataFrame(X, columns=genes, index=adata.obs_names)
    df_expr[cluster_key] = adata.obs[cluster_key].astype(str).values

    avg = df_expr.groupby(cluster_key)[genes].mean().reset_index().to_dict("records")
    pct = (df_expr[genes] > 0).groupby(df_expr[cluster_key]).mean().reset_index().to_dict("records")
    return {"genes": genes, "avg": avg, "pct": pct}



def make_main_umap(umap_df, color_by="annotation", highlight_cluster_id=None):
    id_col = "_cell_id" if "_cell_id" in umap_df.columns else ("barcodes" if "barcodes" in umap_df.columns else None)
    if id_col is None:
        umap_df = umap_df.copy()
        umap_df["_cell_id"] = np.arange(len(umap_df)).astype(str)
        id_col = "_cell_id"

    # base gray (all cells)
    base = px.scatter(
        umap_df,
        x="UMAP1",
        y="UMAP2",
        render_mode="auto",
        hover_name=id_col,
        custom_data=["cluster", "cell_type", "top_markers_str"],
    )
    base.update_traces(
        marker=dict(size=6, opacity=0.18, color="#cbd5e1", line=dict(width=0)),
        hovertemplate=(
            "<b>Cell</b>: %{hovertext}<br>"
            "<b>Cluster</b>: %{customdata[0]}<br>"
            "<b>Cell type</b>: %{customdata[1]}<br>"
            "<br><b>Top markers</b><br>%{customdata[2]}"
            "<extra></extra>"
        ),
        showlegend=False,
    )
    fig = go.Figure(base.data)

    # highlight selected cluster
    if highlight_cluster_id is not None:
        sub = umap_df[umap_df["cluster"] == str(highlight_cluster_id)]
        if not sub.empty:
            is_discrete = (sub[color_by].dtype == "object") or str(sub[color_by].dtype).startswith("category")
            color_map = _discrete_map(umap_df[color_by]) if is_discrete else None

            hi = px.scatter(
                sub,
                x="UMAP1",
                y="UMAP2",
                color=color_by,
                render_mode="auto",
                hover_name=id_col,
                custom_data=["cluster", "cell_type", "top_markers_str"],
                color_discrete_map=color_map,
            )
            hi.update_traces(
                marker=dict(
                    size=8,
                    opacity=0.96,
                    line=dict(width=3, color="rgba(91,124,250,0.25)"),  # halo
                ),
                hovertemplate=(
                    "<b>Cell</b>: %{hovertext}<br>"
                    "<b>Cluster</b>: %{customdata[0]}<br>"
                    "<b>Cell type</b>: %{customdata[1]}<br>"
                    "<br><b>Top markers</b><br>%{customdata[2]}"
                    "<extra></extra>"
                ),
            )

            for tr in hi.data:
                tr.showlegend = True
                fig.add_trace(tr)

    fig.update_layout(
        title=None,
        uirevision="main-umap",
        legend=dict(orientation="h", yanchor="bottom", y=1.01, xanchor="right", x=1),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, title=None),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, title=None),
    )
    return fig


def make_mini_umap(umap_df, height=280):
    mini = px.scatter(
        umap_df,
        x="UMAP1",
        y="UMAP2",
        color="cluster",
        custom_data=["cluster"],
        render_mode="auto",
        color_discrete_map=_discrete_map(umap_df["cluster"]),
    )
    mini.update_layout(
        height=height,
        margin=dict(l=0, r=0, t=0, b=0),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        xaxis=dict(showticklabels=False, showgrid=False, zeroline=False, title=None),
        yaxis=dict(showticklabels=False, showgrid=False, zeroline=False, title=None),
        showlegend=False,
        uirevision="mini-umap",
    )
    mini.update_traces(marker=dict(size=5, opacity=0.80), hovertemplate="Cluster %{customdata[0]}<extra></extra>")
    return mini



# Card / dotplot
def make_explanation_card(cluster_id, marker_info, num_top_k: int = 15):
    info = marker_info.get(cluster_id, {})
    ct = info.get("cell_type", "Unknown")
    expl = info.get("marker_explanations", {})
    top15 = info.get("top15_db", [])
    llm_highlight = set(info.get("llm_explained_genes", []))

    strip = html.Div(
        style={
            "height": "4px",
            "borderRadius": "999px",
            "marginBottom": "14px",
            "opacity": "0.80",
        }
    )
    #strip = html.Div(style={"display": "none"}) / "background": "linear-gradient(90deg,#5b7cfa,#22c55e,#f59e0b)",

    chips = []
    for gene in top15[:num_top_k]:
        gene_str = str(gene)
        icon = html.I(className="fas fa-brain", title="LLM-highlighted marker") if gene_str in llm_highlight else None

        chips.append(
            html.A(
                [gene_str, icon],
                href=gene_url(gene_str),
                target="_blank",
                className="chip chip-link",
            )
        )

    items = [
        html.Li(
            [
                html.A(str(g), href=gene_url(g), target="_blank", className="gene-ref"),
                html.Span(f": {note}"),
            ],
            className="expl-item",
        )
        for g, note in expl.items()
    ]

    return html.Div(
        [
            strip,
            html.Div(
                [
                    html.Div("Selected cluster", className="eyebrow"),
                    html.H3(f"Cluster {cluster_id} → {ct}", className="card-title"),
                ]
            ),
            html.Hr(),
            html.Div(
                [
                    html.Div(f"Top {num_top_k} markers", className="eyebrow"),
                    html.Div(chips if chips else html.Span("No markers found.", className="muted"), className="chip-row"),
                    html.Div("Click gene → open GeneCards reference.", className="hint"),
                ],
                className="section",
            ),
            html.Hr(),
            html.Div(
                [
                    html.Div("Marker explanations", className="eyebrow"),
                    (html.Ul(items, className="expl-list") if items else html.Span("No explanations provided.", className="muted")),
                ],
                className="section",
            ),
        ],
        className="card card--explain",
    )


def make_dotplot_from_cache(cluster_id, marker_info, cache, topk=15):
    info = marker_info.get(cluster_id, {})
    llm, dbk = info.get("llm_markers", []), info.get("top15_db", [])
    all_genes_for_cluster = llm + [g for g in dbk if g not in llm]
    genes_to_plot = [g for g in all_genes_for_cluster if g in cache.get("genes", [])][:topk]

    if not genes_to_plot:
        return go.Figure().update_layout(
            annotations=[dict(text="No marker data available.", showarrow=False, font=dict(size=14, color="#64748b"))],
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            plot_bgcolor="#ffffff",
        )

    avg, pct = pd.DataFrame(cache["avg"]), pd.DataFrame(cache["pct"])
    avg_sub = avg.set_index("cluster")[genes_to_plot]
    pct_sub = pct.set_index("cluster")[genes_to_plot]

    def wrap_celltype(text, width=20):
        import textwrap
        return "<br>".join(textwrap.wrap(str(text), width=width))

    y_labels = []
    for cid in avg_sub.index:
        ct = marker_info.get(str(cid), {}).get("cell_type", "Unknown")
        ct_wrapped = wrap_celltype(ct, width=20)
        y_labels.append(f"Cluster {cid}<br><span style='font-size:17px;color:#64748b'>{ct_wrapped}</span>")

    avg_sub.index = pct_sub.index = y_labels

    df_long = (
        avg_sub.reset_index()
        .melt("index", var_name="Gene", value_name="AvgExpr")
        .merge(pct_sub.reset_index().melt("index", var_name="Gene", value_name="PctExpr"), on=["index", "Gene"])
    )

    orrd_scale = px.colors.sequential.OrRd
    custom_colorscale = [[0.0, "#ffffff"], [1e-9, orrd_scale[0]]]
    custom_colorscale.extend([[i / (len(orrd_scale) - 1), color] for i, color in enumerate(orrd_scale)])

    fig = px.scatter(
        df_long,
        x="Gene",
        y="index",
        size="PctExpr",
        color="AvgExpr",
        color_continuous_scale=custom_colorscale,
        size_max=16,
        custom_data=["AvgExpr", "PctExpr"],
        title=None,
    )

    fig.update_traces(
        hovertemplate=(
            "<b>Gene</b>: %{x}<br>"
            "<b>Cluster</b>: %{y}<br><br>"
            "Avg. Expression: %{customdata[0]:.3f}<br>"
            "Pct. Expressed: %{customdata[1]:.1%}"
            "<extra></extra>"
        ),
        marker=dict(line=dict(width=1, color="#7f1d1d")),
    )

    fig.update_layout(
        template="cap_light_pro",
        xaxis_title="Gene",
        yaxis_title=None,
        margin=dict(l=185, r=0, t=24, b=0),
        xaxis=dict(
            showticklabels=True,
            tickangle=-45,
            automargin=True,
            ticks="outside",
            ticklen=4,
            showgrid=False,
            zeroline=False,
            tickfont=dict(size=16),
        ),
        yaxis=dict(
            showticklabels=True,
            automargin=False,
            ticks="outside",
            ticklen=4,
            showgrid=False,
            zeroline=False,
            tickfont=dict(size=14, family="Inter, Arial"),
        ),
        coloraxis_colorbar=dict(title="Avg.<br>Expr.", ticks="outside", ticklen=4, outlinewidth=0),
        uirevision="dot-plot",
    )

    return fig

# Rationale loader
def _load_rationale_file(path: str) -> Optional[dict]:
    if not path or not os.path.exists(path):
        return None
    ext = os.path.splitext(path)[1].lower()
    try:
        if ext in (".json", ".jsonl"):
            with open(path, "r", encoding="utf-8") as f:
                data = json.load(f)
            if isinstance(data, dict):
                return {str(k): v for k, v in data.items()}
            elif isinstance(data, list):
                df = pd.DataFrame(data)
            else:
                return None
        else:
            df = pd.read_csv(path, sep="\t" if ext in (".tsv", ".tab") else ",", encoding="utf-8")
    except Exception as e:
        print(f"Error loading rationale file {path}: {e}")
        return None

    if "cluster" not in df.columns:
        return None

    df["cluster"] = df["cluster"].astype(str)
    out = {}
    for _, row in df.iterrows():
        cid = row["cluster"]
        expl = _parse_marker_expl(row["marker_explanations"]) if ("marker_explanations" in df.columns and pd.notna(row.get("marker_explanations"))) else {}
        out[cid] = {"cell_type": row.get("cell_type"), "marker_explanations": expl}
    return out


# App launcher
def launch_cap_style_app(
    adata,
    use_uns_markers: str = "marker_list",
    use_uns_llm_annotation: str = "GPT_annotation_db",
    port=8051,
    debug=True,
    num_top_k: int = 15,
    rationale_json_path: Optional[str] = None,
    cluster_key: str = "cluster",
):
    rj = _load_rationale_file(rationale_json_path)

    marker_info = build_marker_info_from_uns(
        adata,
        use_uns_markers,
        use_uns_llm_annotation,
        cluster_key=cluster_key,
        num_top_k=num_top_k,
        rationale_json=rj,
    )

    umap_df = make_umap_df(adata, marker_info, cluster_key=cluster_key)
    cache_stats = precompute_cluster_gene_stats_dense(adata, marker_info, cluster_key=cluster_key)

    app = Dash(__name__, external_stylesheets=["https://use.fontawesome.com/releases/v5.15.4/css/all.css"])
    app.title = "Cell Annotation Explorer"

    MINI_CONTAINER_H = 300
    MINI_FIG_H = 280

    app.index_string = f"""
<!DOCTYPE html>
<html>
<head>
  {{%metas%}}
  <title>Cell Annotation Explorer</title>
  {{%favicon%}}
  {{%css%}}
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">

  <style>
    :root {{
      --bg: #f5f7fb;
      --panel: rgba(255,255,255,0.82);
      --ink: #0b1220;
      --muted: #64748b;

      --accent: #5b7cfa;
      --accent-soft: rgba(91,124,250,0.15);
      --accent-glow: rgba(91,124,250,0.35);

      --border: rgba(15,23,42,0.10);
      --border-strong: rgba(15,23,42,0.18);

      --shadow-xs: 0 1px 1px rgba(15,23,42,0.04);
      --shadow-sm: 0 2px 8px rgba(15,23,42,0.06);
      --shadow-md: 0 16px 40px rgba(15,23,42,0.12);

      --radius-sm: 10px;
      --radius-md: 16px;
      --radius-lg: 22px;

      --font: 'Inter', ui-sans-serif, system-ui;
    }}

    html, body {{
      height: 100vh;
      margin: 0;
      background: var(--bg);
      color: var(--ink);
      font-family: var(--font);
      overflow: hidden;
      -webkit-font-smoothing: antialiased;
      -moz-osx-font-smoothing: grayscale;
    }}

    body::before {{
      content: "";
      position: fixed;
      inset: -140px;
      z-index: -1;
      background:
        radial-gradient(800px 420px at 10% 8%, rgba(91,124,250,0.18), transparent 60%),
        radial-gradient(700px 360px at 85% 10%, rgba(34,197,94,0.14), transparent 55%),
        radial-gradient(620px 320px at 80% 90%, rgba(251,191,36,0.16), transparent 60%);
      filter: blur(6px);
    }}

    .app-container {{
      display: grid;
      grid-template-columns: 360px 1fr;
      height: 100vh;
    }}

    .sidebar {{
      padding: 20px;
      border-right: 1px solid var(--border);
      background: linear-gradient(180deg, rgba(255,255,255,0.88), rgba(255,255,255,0.72));
      backdrop-filter: blur(12px);
      overflow-y: auto;
    }}

    .brand {{
      display: flex;
      gap: 14px;
      align-items: center;
      padding: 16px;
      border-radius: var(--radius-lg);
      background: linear-gradient(135deg, rgba(255,255,255,0.95), rgba(255,255,255,0.75));
      border: 1px solid var(--border);
      box-shadow: var(--shadow-sm);
    }}

    .brand-badge {{
      width: 44px;
      height: 44px;
      border-radius: 14px;
      background: radial-gradient(circle at 30% 30%, #ffffff, var(--accent-soft)),
                  linear-gradient(135deg, var(--accent), #8ea0ff);
      color: white;
      font-weight: 900;
      display: flex;
      align-items: center;
      justify-content: center;
      box-shadow: inset 0 0 0 1px rgba(255,255,255,0.5),
                  0 6px 16px var(--accent-soft);
    }}

    .brand h2 {{
      font-size: 1.05rem;
      margin: 0;
      letter-spacing: -0.02em;
    }}

    .brand p {{
      margin: 6px 0 0;
      font-size: 0.88rem;
      color: var(--muted);
      line-height: 1.4;
    }}

    .section-title {{
      margin: 18px 6px 10px;
      font-size: 0.72rem;
      text-transform: uppercase;
      letter-spacing: .14em;
      color: var(--muted);
      font-weight: 800;
    }}

    .card {{
      background: var(--panel);
      border-radius: var(--radius-lg);
      border: 1px solid var(--border);
      box-shadow: var(--shadow-sm);
      position: relative;
      transition: transform .22s ease, box-shadow .22s ease, border-color .22s ease;
    }}

    .card::after {{
      content: "";
      position: absolute;
      inset: 0;
      border-radius: inherit;
      pointer-events: none;
      box-shadow: inset 0 1px 0 rgba(255,255,255,0.65);
    }}

    .card:hover {{
      box-shadow: var(--shadow-md);
      transform: translateY(-1px);
      border-color: var(--border-strong);
    }}

    .pad {{ padding: 14px; }}

    .control-group {{
      display: flex;
      flex-direction: column;
      gap: 14px;
    }}

    .hint {{
      color: var(--muted);
      font-size: 12px;
      margin-top: 8px;
    }}

    #mini-umap-container {{
      height: {MINI_CONTAINER_H}px;
      padding: 10px;
      border-radius: var(--radius-lg);
      background: linear-gradient(180deg, rgba(255,255,255,0.92), rgba(255,255,255,0.70));
      border: 1px solid var(--border);
      box-shadow: var(--shadow-sm);
    }}
    #mini-umap-container .dash-graph {{ height: 100%; }}

    .main-content {{
      padding: 22px 26px;
      overflow-y: auto;
    }}

    .topbar {{
      display: flex;
      justify-content: space-between;
      align-items: flex-end;
      margin-bottom: 16px;
      gap: 12px;
    }}

    .topbar h1 {{
      font-size: 1.18rem;
      margin: 0;
      letter-spacing: -0.02em;
    }}

    .topbar .subtitle {{
      font-size: 0.90rem;
      color: var(--muted);
      margin-top: 6px;
    }}

    .badge {{
      padding: 8px 14px;
      border-radius: 999px;
      background: linear-gradient(135deg, #ffffff, rgba(255,255,255,0.65));
      border: 1px solid var(--border);
      font-size: 12px;
      font-weight: 800;
      box-shadow: var(--shadow-xs);
      white-space: nowrap;
    }}

    .main-grid {{
      display: grid;
      grid-template-columns: 1.5fr 3.5fr;
      gap: 26px;
      align-items: start;
    }}

    .graph-card {{ padding: 16px; }}

    .graph-title {{
      font-size: 0.95rem;
      font-weight: 900;
      margin: 0 0 10px 6px;
      letter-spacing: -0.01em;
    }}

    .card--explain {{
      margin-top: 16px;
      padding: 18px;
    }}

    .eyebrow {{
      font-size: 0.72rem;
      color: var(--muted);
      text-transform: uppercase;
      letter-spacing: .12em;
      font-weight: 900;
    }}

    .card-title {{
      margin: 8px 0 0;
      font-size: 1.05rem;
      font-weight: 900;
      letter-spacing: -0.02em;
    }}

    hr {{
      border: none;
      height: 1px;
      background: rgba(15,23,42,0.10);
      margin: 14px 0;
    }}

    .chip-row {{
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      margin-top: 10px;
    }}

    .chip {{
      display: inline-flex;
      align-items: center;
      gap: 8px;
      padding: 6px 12px;
      border-radius: 999px;
      background: linear-gradient(135deg, rgba(255,255,255,0.90), rgba(255,255,255,0.65));
      border: 1px solid var(--border);
      font-size: 0.80rem;
      font-weight: 650;
      box-shadow: var(--shadow-xs);
      transition: transform .14s ease, box-shadow .14s ease, border-color .14s ease;
      color: inherit;
    }}

    .chip:hover {{
      transform: translateY(-1px);
      box-shadow: var(--shadow-sm);
      border-color: rgba(91,124,250,0.25);
    }}

    /* clickable gene chip */
    .chip-link {{
      text-decoration: none;
      cursor: pointer;
    }}
    .chip-link:hover {{
      background: rgba(91,124,250,0.18);
      box-shadow: 0 0 0 1px rgba(91,124,250,0.25);
    }}

    /* gene link inside explanations */
    .gene-ref {{
      color: var(--ink);
      text-decoration: none;
      font-weight: 800;
    }}
    .gene-ref:hover {{
      text-decoration: underline;
      color: var(--accent);
    }}

    .chip .fa-brain {{
      color: var(--accent);
      text-shadow: 0 0 8px rgba(91,124,250,0.22);
      opacity: 0.95;
    }}

    .expl-list {{
      margin: 10px 0 0;
      padding: 0 0 0 18px;
      list-style-type: '— ';
    }}

    .expl-item {{
      margin: 6px 0;
      line-height: 1.55;
      padding-left: 8px;
    }}

    .muted {{ color: var(--muted); }}

    .js-plotly-plot .modebar {{
      opacity: 0.08;
      transition: opacity .2s ease;
    }}
    .js-plotly-plot:hover .modebar {{
      opacity: 0.45;
    }}

    @media (max-width: 1200px) {{
      .app-container {{ grid-template-columns: 1fr; }}
      .sidebar {{ border-right: none; border-bottom: 1px solid var(--border); }}
      .main-grid {{ grid-template-columns: 1fr; }}
    }}
  </style>
</head>
<body>
  {{%app_entry%}}
  <footer>{{%config%}}{{%scripts%}}{{%renderer%}}</footer>
</body>
</html>
"""

    init_main_fig = make_main_umap(umap_df, color_by="annotation", highlight_cluster_id=None)
    init_mini_fig = make_mini_umap(umap_df, height=MINI_FIG_H)

    sidebar = html.Div(
        [
            html.Div(
                [
                    html.Div("CA", className="brand-badge"),
                    html.Div(
                        [
                            html.H2("Annotation Explorer"),
                            html.P("Select clusters from the mini UMAP to update the main UMAP and marker panels."),
                        ]
                    ),
                ],
                className="brand",
            ),
            html.Div("Display options", className="section-title"),
            html.Div(
                [
                    html.Div(
                        className="control-group",
                        children=[
                            dcc.Dropdown(
                                id="color-by",
                                options=[{"label": x.replace("_", " ").title(), "value": x} for x in ["annotation", "cluster", "cell_type"]],
                                value="annotation",
                                clearable=False,
                            ),
                            html.Div(
                                [
                                    html.Div("Top K genes in dot plot", style={"fontSize": "12px", "fontWeight": "800", "color": "#334155"}),
                                    dcc.Slider(
                                        id="topk",
                                        min=5, max=15, step=1, value=15,
                                        marks={5: "5", 10: "10", 15: "15"},
                                        tooltip={"placement": "bottom"},
                                    ),
                                    html.Div("Tip: Use mini UMAP for precise cluster selection.", className="hint"),
                                ]
                            ),
                        ],
                    ),
                ],
                className="card pad",
            ),
            html.Div("Navigation", className="section-title"),
            html.Div(
                dcc.Graph(
                    id="mini-umap",
                    figure=init_mini_fig,
                    style={"height": "100%", "width": "100%"},
                    config={"displayModeBar": False, "responsive": True},
                    clear_on_unhover=True,
                ),
                id="mini-umap-container",
            ),
        ],
        className="sidebar",
    )

    main_content = html.Div(
        [
            html.Div(
                [
                    html.Div(
                        [
                            html.H1("Cell annotation explorer"),
                            html.Div("Selected cluster is emphasized; non-selected cells remain light gray.", className="subtitle"),
                        ]
                    ),
                    html.Div(id="selected-badge", className="badge"),
                ],
                className="topbar",
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.H3("UMAP", className="graph-title"),
                            html.Div(
                                dcc.Graph(
                                    id="main-umap",
                                    figure=init_main_fig,
                                    style={"height": "100%", "width": "100%"},
                                    config={"displayModeBar": False, "responsive": True},
                                    clear_on_unhover=True,
                                ),
                                style={"height": "500px"},
                            ),
                        ],
                        className="card graph-card",
                    ),
                    html.Div(
                        [
                            html.H3("Marker gene expression", className="graph-title"),
                            html.Div(dcc.Graph(id="dot-plot"), style={"height": "500px"}),
                        ],
                        className="card graph-card",
                    ),
                ],
                className="main-grid",
            ),
            html.Div(id="marker-explanation-box", className="card card--explain"),
        ],
        className="main-content",
    )

    app.layout = html.Div(
        [
            dcc.Store(id="selected_cluster", data=None),
            sidebar,
            main_content,
        ],
        className="app-container",
    )

    # selection from secondary umap
    @app.callback(
        Output("selected_cluster", "data"),
        Input("mini-umap", "clickData"),
        State("selected_cluster", "data"),
        prevent_initial_call=True,
    )
    def set_selected_cluster(mini_click, current):
        if not mini_click:
            return no_update
        try:
            return str(mini_click["points"][0]["customdata"][0])
        except Exception:
            return no_update

    @app.callback(Output("selected-badge", "children"), Input("selected_cluster", "data"))
    def update_badge(cid):
        return f"Selected: Cluster {cid}" if cid else "Selected: —"

    @app.callback(
        Output("main-umap", "figure"),
        [Input("color-by", "value"), Input("selected_cluster", "data")],
    )
    def update_main_umap(color_by, selected_cluster):
        try:
            return make_main_umap(umap_df, color_by=color_by, highlight_cluster_id=selected_cluster)
        except Exception as e:
            print("[update_main_umap ERROR]", repr(e))
            return make_main_umap(umap_df, color_by="annotation", highlight_cluster_id=selected_cluster)

    @app.callback(
        [Output("marker-explanation-box", "children"), Output("dot-plot", "figure")],
        [Input("selected_cluster", "data"), Input("topk", "value")],
    )
    def update_details(selected_cluster, topk):
        if not selected_cluster:
            initial_card = html.Div(
                [
                    html.Div("Getting started", className="eyebrow"),
                    html.P(
                        "Click a cluster in the navigation mini UMAP (left).",
                        style={"marginTop": "10px", "lineHeight": "1.7", "color": "#334155"},
                    ),
                    html.Div("Brain icon indicates LLM-highlighted markers.", className="hint", style={"marginTop": "10px"}),
                ]
            )
            fig = go.Figure().update_layout(
                plot_bgcolor="#ffffff",
                annotations=[dict(text="Select a cluster to view dot plot", showarrow=False, font=dict(size=14, color="#64748b"))],
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
            )
            return initial_card, fig

        cid = str(selected_cluster)
        card = make_explanation_card(cid, marker_info, num_top_k=num_top_k)
        fig = make_dotplot_from_cache(cid, marker_info, cache_stats, topk=int(topk))
        return card, fig

    app.run(debug=debug, port=port)