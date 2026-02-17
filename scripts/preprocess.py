import scanpy as sc
from pathlib import Path


def run_preprocessing(
    adata,
    sample_id,
    results_dir,
    n_top_genes=3000,
    n_pcs=20
):
    """
    Preprocessing stage for spatial transcriptomics data.

    Steps:
    1. Normalize total counts
    2. Log-transform (ln(1+x))
    3. Save log-normalized data to .raw
    4. Identify highly variable genes
    5. Subset to HVGs
    6. Scale data
    7. Run PCA
    8. Save plots and processed AnnData object
    """

    print(f"\n🔹 Preprocessing sample {sample_id}...")

    # ----------------------------
    # 1️⃣ Normalize total counts
    # ----------------------------
    sc.pp.normalize_total(adata, target_sum=1e4)

    # ----------------------------
    # 2️⃣ Log transform (ln(1+x))
    # ----------------------------
    sc.pp.log1p(adata)

    # ----------------------------
    # ✅ Save log-normalized data
    # (used later for marker analysis)
    # ----------------------------
    adata.raw = adata.copy()

    # ----------------------------
    # 3️⃣ Highly Variable Genes
    # ----------------------------
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat",
        n_top_genes=n_top_genes
    )

    print(f"   Number of HVGs selected: {adata.var.highly_variable.sum()}")

    # ----------------------------
    # Create figure directory
    # ----------------------------
    figures_path = results_dir / "figures" / "preprocessing"
    figures_path.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = figures_path

    # ----------------------------
    # Save HVG plot
    # ----------------------------
    sc.pl.highly_variable_genes(
        adata,
        save=f"_{sample_id}_hvg.png",
        show=False
    )

    # ----------------------------
    # 4️⃣ Keep only HVGs
    # ----------------------------
    adata = adata[:, adata.var.highly_variable].copy()

    # ----------------------------
    # 5️⃣ Scale data
    # ----------------------------
    sc.pp.scale(adata, max_value=10)

    # ----------------------------
    # 6️⃣ PCA
    # ----------------------------
    sc.tl.pca(
        adata,
        svd_solver="arpack",
        n_comps=n_pcs
    )

    # ----------------------------
    # Save PCA variance plot
    # ----------------------------
    sc.pl.pca_variance_ratio(
        adata,
        log=True,
        save=f"_{sample_id}_pca_variance.png",
        show=False
    )

    # ----------------------------
    # Save processed AnnData object
    # ----------------------------
    processed_path = results_dir / "processed"
    processed_path.mkdir(parents=True, exist_ok=True)

    output_file = processed_path / f"{sample_id}_processed.h5ad"
    adata.write(output_file)

    print(f"   Preprocessing complete for {sample_id}")
    print(f"   Saved: {output_file}")

    return adata
