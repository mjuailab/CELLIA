import argparse
from pathlib import Path
import scanpy as sc

def build_parser():
    p = argparse.ArgumentParser(
        description="Run CELLIA workflow"
    )

    p.add_argument("--input", "-i", required=True, help="Input AnnData (.h5ad)")
    p.add_argument("--output", "-o", required=True, help="Output AnnData (.h5ad)")

    p.add_argument("--tissue_db", required=True, help="DB tissue string for Marker_DB filtering (e.g., PBMC)")
    p.add_argument("--tissue_type", required=True, help="Prompt tissue context (e.g., human PBMC)")
    p.add_argument("--api_key", required=True, help="API key for the selected LLM provider")

    p.add_argument("--groupby", default="cluster", help="adata.obs column for clustering (default: cluster)")
    p.set_defaults(tie_correct=True)

    p.add_argument("--n_top_markers", type=int, default=15, help="Top-k markers per cluster (default: 15)")

    p.add_argument("--deg_mode", choices=["major", "subset"], default="major", help="DEG threshold preset")
    p.add_argument("--db_mode", choices=["db", "subset_db"], default="db", help="DB filtering preset")
    p.add_argument("--subset_db", default=None, help="Subset cell type string used when db-mode=subset_db")
    p.add_argument("--db_path", default="./database/Marker_DB.csv", help="Marker DB csv path")

    p.add_argument("--llm_provider", choices=["gpt", "gemini", "claude"], default="gpt", help="LLM provider")
    p.add_argument("--llm_mode", choices=["major", "subset"], default="major", help="Prompt mode")
    p.add_argument("--model", default=None, help="Provider model name")
    p.add_argument("--parent_celltype", default=None, help="Required if llm-mode==subset")
    p.add_argument

    return p


def main():
    args = build_parser().parse_args()

    adata = sc.read_h5ad(args.input)

    if args.llm_mode == "subset" and not args.parent_celltype:
        raise SystemExit("--parent-celltype is required when --llm-mode subset")

    from cellia import cellia_run 

    adata = cellia_run(
        adata=adata,
        tissue_db=args.tissue_db,
        tissue_type=args.tissue_type,
        api_key=args.api_key,
        n_top_markers=args.n_top_markers,
        groupby=args.groupby,
        deg_mode=args.deg_mode,
        db_mode=args.db_mode,
        db_path=args.db_path,
        llm_provider=args.llm_provider,
        llm_mode=args.llm_mode,
        model=args.model,
        parent_celltype=args.parent_celltype,
        subset_db=args.subset_db,
    )

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path)

    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()