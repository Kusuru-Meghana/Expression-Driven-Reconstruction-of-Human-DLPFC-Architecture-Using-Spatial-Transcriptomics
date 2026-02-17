import scanpy as sc
from pathlib import Path


def run_marker_analysis(adata, sample_id, results_dir, n_top=20):
    """
    Marker gene analysis per Leiden cluster.

    Steps:
    1. Rank genes per cluster
    2. Save marker table
    3. Save dotplot + heatmap
    """

    print(f"\n🔹 Running marker analysis for {sample_id}...")

    # ----------------------------
    # 1️⃣ Differential expression
    # ----------------------------
    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method="wilcoxon",
        use_raw=True
    )
    marker_df = sc.get.rank_genes_groups_df(adata, None)

    # ----------------------------
    # Create results directories
    # ----------------------------
    marker_dir = results_dir / "markers"
    fig_dir = results_dir / "figures" / "markers"

    marker_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    sc.settings.figdir = fig_dir

    # ----------------------------
    # 2️⃣ Save marker table
    # ----------------------------
    marker_df = sc.get.rank_genes_groups_df(adata, None)

    output_csv = marker_dir / f"{sample_id}_markers.csv"
    marker_df.to_csv(output_csv, index=False)

    print(f"   Marker table saved: {output_csv}")

    # ----------------------------
    # 3️⃣ Save Dotplot
    # ----------------------------
    sc.pl.rank_genes_groups_dotplot(
        adata,
        n_genes=n_top,
        show=False,
        save=f"_{sample_id}_dotplot.png"
    )

    # ----------------------------
    # 4️⃣ Save Heatmap
    # ----------------------------
    sc.pl.rank_genes_groups_heatmap(
        adata,
        n_genes=n_top,
        show=False,
        save=f"_{sample_id}_heatmap.png"
    )

    print(f"   Marker plots saved for {sample_id}")

    return adata
