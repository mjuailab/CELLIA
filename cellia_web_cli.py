import argparse
from pathlib import Path
import scanpy as sc
from cellia_web import launch_cap_style_app


def build_parser():
    p = argparse.ArgumentParser(description="Launch CELLIA Annotation Explorer (Dash)")

    p.add_argument("-i", "--input", required=True, help="Input AnnData (.h5ad)")
    p.add_argument("--port", type=int, default=8051, help="Dash server port (default: 8051)")
    p.add_argument("--debug", action="store_true", help="Enable debug mode")
    p.add_argument("--no_debug", dest="debug", action="store_false")
    p.set_defaults(debug=True)

    p.add_argument("--cluster_key", default="cluster", help="adata.obs key for cluster (default: cluster)")
    p.add_argument("--markers_uns", default="marker_list", help="adata.uns key for marker table (default: marker_list)")
    p.add_argument(
        "--llm_uns",
        default="GPT_annotation_db",
        help="adata.uns key for LLM annotation table (default: GPT_annotation_db)",
    )

    p.add_argument("--topk", type=int, default=15, help="Top K genes displayed in marker dot plot (default: 15)")
    p.add_argument(
        "--rationale",
        type=str,
        default=None,
        help="Path to JSON file containing LLM explanations/rationales.",
    )
    return p


def main():
    args = build_parser().parse_args()

    in_path = Path(args.input)
    if not in_path.exists():
        raise SystemExit(f"[ERROR] Input file not found: {in_path}")

    adata = sc.read_h5ad(in_path)

    if args.cluster_key not in adata.obs:
        raise SystemExit(f"[ERROR] cluster_key '{args.cluster_key}' not found in adata.obs")

    if args.markers_uns not in adata.uns:
        raise SystemExit(
            f"[ERROR] markers_uns '{args.markers_uns}' not found in adata.uns. "
            f"Available uns keys: {list(adata.uns.keys())[:30]}"
        )

    print("[INFO] Launching app with:")
    print(f"  input       = {in_path}")
    print(f"  port        = {args.port}")
    print(f"  debug       = {args.debug}")
    print(f"  cluster_key = {args.cluster_key}")
    print(f"  markers_uns = {args.markers_uns}")
    print(f"  llm_uns     = {args.llm_uns}")
    print(f"  topk        = {args.topk}")
    print(f"  rationale   = {args.rationale}")

    launch_cap_style_app(
        adata=adata,
        use_uns_markers=args.markers_uns,
        use_uns_llm_annotation=args.llm_uns,
        port=args.port,
        debug=args.debug,
        num_top_k=args.topk,
        rationale_json_path=args.rationale,
        cluster_key=args.cluster_key,
    )


if __name__ == "__main__":
    main()