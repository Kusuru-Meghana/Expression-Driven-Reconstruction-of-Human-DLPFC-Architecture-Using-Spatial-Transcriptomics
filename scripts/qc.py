# ============================================
# QC Module for Visium Spatial Transcriptomics
# ============================================

import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path


def run_qc(adata, sample_id, figures_dir):
    """
    Run full QC pipeline:
    - Calculate QC metrics
    - Pre-filter summary + plots
    - Apply filtering thresholds
    - Post-filter summary + plots
    - Return filtered AnnData
    """

    print("\n====================================")
    print(f"Running QC for Sample: {sample_id}")
    print("====================================")

    # ----------------------------------------
    # 1️⃣ Identify mitochondrial genes
    # ----------------------------------------
    adata.var["mt"] = adata.var_names.str.startswith("MT-")

    # ----------------------------------------
    # 2️⃣ Calculate QC metrics
    # ----------------------------------------
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        inplace=True
    )

    metrics = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]

    figures_dir = Path(figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)

    # ========================================
    # 🔵 PRE-FILTER SUMMARY
    # ========================================
    print("\n--- PRE-FILTER QC SUMMARY ---")

    for metric in metrics:
        print(f"\n{metric}")
        print(f"   Mean : {adata.obs[metric].mean():.2f}")
        print(f"   Min  : {adata.obs[metric].min():.2f}")
        print(f"   Max  : {adata.obs[metric].max():.2f}")

    # ========================================
    # 🔵 PRE-FILTER PLOTS
    # ========================================
    sc.pl.violin(adata, metrics, jitter=0.4, multi_panel=True, show=False)
    plt.suptitle(f"PRE-FILTER QC - Sample {sample_id}")
    plt.savefig(figures_dir / f"{sample_id}_PRE_QC_violin.png")
    plt.close()

    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", show=False)
    plt.title(f"PRE Total Counts vs Genes - {sample_id}")
    plt.savefig(figures_dir / f"{sample_id}_PRE_counts_vs_genes.png")
    plt.close()

    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", show=False)
    plt.title(f"PRE Total Counts vs Mito - {sample_id}")
    plt.savefig(figures_dir / f"{sample_id}_PRE_counts_vs_mito.png")
    plt.close()

    # ========================================
    # 🟡 FILTERING
    # ========================================
    print("\nApplying filtering thresholds...")

    before = adata.n_obs

    adata = adata[
        (adata.obs["n_genes_by_counts"] > 200) &
        (adata.obs["n_genes_by_counts"] < 4500) &
        (adata.obs["pct_counts_mt"] < 30)
    ].copy()

    after = adata.n_obs

    print(f"Spots before filtering: {before}")
    print(f"Spots after filtering : {after}")
    print(f"Spots removed         : {before - after}")

    # ========================================
    # 🔴 POST-FILTER SUMMARY
    # ========================================
    print("\n--- POST-FILTER QC SUMMARY ---")

    for metric in metrics:
        print(f"\n{metric}")
        print(f"   Mean : {adata.obs[metric].mean():.2f}")
        print(f"   Min  : {adata.obs[metric].min():.2f}")
        print(f"   Max  : {adata.obs[metric].max():.2f}")

    # ========================================
    # 🔴 POST-FILTER PLOTS
    # ========================================
    sc.pl.violin(adata, metrics, jitter=0.4, multi_panel=True, show=False)
    plt.suptitle(f"POST-FILTER QC - Sample {sample_id}")
    plt.savefig(figures_dir / f"{sample_id}_POST_QC_violin.png")
    plt.close()

    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", show=False)
    plt.title(f"POST Total Counts vs Genes - {sample_id}")
    plt.savefig(figures_dir / f"{sample_id}_POST_counts_vs_genes.png")
    plt.close()

    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", show=False)
    plt.title(f"POST Total Counts vs Mito - {sample_id}")
    plt.savefig(figures_dir / f"{sample_id}_POST_counts_vs_mito.png")
    plt.close()

    print("\nQC completed for sample:", sample_id)
    print("====================================\n")

    return adata
