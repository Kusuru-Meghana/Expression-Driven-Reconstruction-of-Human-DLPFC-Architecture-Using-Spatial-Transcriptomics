import scanpy as sc
from pathlib import Path


def integrate_samples(
    results_dir: Path,
    sample_ids: list,
    n_top_genes: int = 3000,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    resolution: float = 0.5,
    random_state: int = 0,
):
    """
    Multi-sample integration of Visium spatial transcriptomics data.

    Loads *_filtered.h5ad files from results/,
    integrates them using batch-aware HVGs + ComBat,
    then performs PCA, Leiden clustering, and UMAP.

    Outputs:
    - results/integrated/integrated.h5ad
    - results/figures/integration/umap_integrated.png
    - results/figures/integration/spatial_by_sample/*
    - results/markers/integrated_markers.csv
    """

    results_dir = Path(results_dir)

    # --------------------------------------------
    # Output directories
    # --------------------------------------------
    integrated_dir = results_dir / "integrated"
    fig_dir = results_dir / "figures" / "integration"
    spatial_dir = fig_dir / "spatial_by_sample"
    marker_dir = results_dir / "markers"

    integrated_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    spatial_dir.mkdir(parents=True, exist_ok=True)
    marker_dir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------
    # 1️⃣ Load filtered objects
    # --------------------------------------------
    adatas = []

    for sid in sample_ids:
        file_path = results_dir / f"{sid}_filtered.h5ad"

        if not file_path.exists():
            raise FileNotFoundError(
                f"Missing file: {file_path}\n"
                f"Expected filtered files inside results/"
            )

        adata = sc.read_h5ad(file_path)
        adata.obs["sample_id"] = sid
        adatas.append(adata)

    # --------------------------------------------
    # 2️⃣ Concatenate samples
    # --------------------------------------------
    adata = sc.concat(
        adatas,
        label="sample_id",
        keys=sample_ids,
        join="outer",
        index_unique="-",
        merge="unique",
    )

    print("Integrated shape:", adata.shape)

    # --------------------------------------------
    # 3️⃣ Normalize + log transform (joint)
    # --------------------------------------------
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Save full log-normalized data for markers later
    adata.raw = adata.copy()

    # --------------------------------------------
    # 4️⃣ Batch-aware HVG selection
    # --------------------------------------------
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat",
        n_top_genes=n_top_genes,
        batch_key="sample_id",
    )

    print("Joint HVGs selected:", adata.var["highly_variable"].sum())

    # Subset to HVGs
    adata = adata[:, adata.var["highly_variable"]].copy()

    # --------------------------------------------
    # 5️⃣ Batch correction (ComBat)
    # --------------------------------------------
    sc.pp.combat(adata, key="sample_id")

    # --------------------------------------------
    # 6️⃣ Scale + PCA
    # --------------------------------------------
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(
        adata,
        n_comps=n_pcs,
        svd_solver="arpack",
        random_state=random_state,
    )

    # --------------------------------------------
    # 7️⃣ Neighbors + Leiden + UMAP
    # --------------------------------------------
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        random_state=random_state,
    )

    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added="integrated_leiden",
        flavor="igraph",
        n_iterations=2,
        directed=False,
        random_state=random_state,
    )

    sc.tl.umap(adata, random_state=random_state)

    print("Integrated clusters found:",
          adata.obs["integrated_leiden"].nunique())

    # --------------------------------------------
    # 8️⃣ Save UMAP
    # --------------------------------------------
    sc.settings.figdir = fig_dir

    sc.pl.umap(
        adata,
        color=["integrated_leiden", "sample_id"],
        wspace=0.4,
        show=False,
        save="_integrated.png",
    )

    # --------------------------------------------
    # 9️⃣ Spatial plots per sample
    # --------------------------------------------
    for sid in sample_ids:
        sub = adata[adata.obs["sample_id"] == sid].copy()
        sc.settings.figdir = spatial_dir

        sc.pl.spatial(
            sub,
            color="integrated_leiden",
            spot_size=1.2,
            show=False,
            save=f"_{sid}_integrated_clusters.png",
        )

    # --------------------------------------------
    # 🔟 Marker analysis (integrated clusters)
    # --------------------------------------------
    sc.tl.rank_genes_groups(
        adata,
        groupby="integrated_leiden",
        method="wilcoxon",
        use_raw=False,
    )

    marker_df = sc.get.rank_genes_groups_df(adata, None)
    marker_df.to_csv(marker_dir / "integrated_markers.csv", index=False)

    sc.settings.figdir = fig_dir
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby="integrated_leiden",
        n_genes=10,
        show=False,
        save="_integrated_markers.png",
    )

    # --------------------------------------------
    # 1️⃣1️⃣ Save integrated object
    # --------------------------------------------
    output_file = integrated_dir / "integrated.h5ad"
    adata.write(output_file)

    print("Saved integrated object:", output_file)

    return adata
