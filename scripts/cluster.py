import scanpy as sc
from pathlib import Path


def run_clustering(
    adata,
    sample_id,
    results_dir,
    n_neighbors=15,
    n_pcs=20,
    resolution=0.5
):

    print(f"\n🔹 Clustering sample {sample_id}...")

    # ----------------------------
    # 1️⃣ Compute neighbors (PCA space)
    # ----------------------------
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        random_state=0   # ✅ ADD HERE
    )

    # ----------------------------
    # 2️⃣ Leiden clustering
    # ----------------------------
    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added="leiden",
        flavor="igraph",       # future-proof
        n_iterations=2,
        directed=False,
        random_state=0         # ✅ ADD HERE
    )

    print(f"   Number of clusters found: {adata.obs['leiden'].nunique()}")

    # ----------------------------
    # 3️⃣ UMAP
    # ----------------------------
    sc.tl.umap(
        adata,
        random_state=0    # ✅ ADD HERE
    )

    # ----------------------------
    # Create figure directories
    # ----------------------------
    umap_dir = results_dir / "figures" / "clustering" / "umap"
    spatial_dir = results_dir / "figures" / "clustering" / "spatial"

    umap_dir.mkdir(parents=True, exist_ok=True)
    spatial_dir.mkdir(parents=True, exist_ok=True)

    # ----------------------------
    # Save UMAP plot
    # ----------------------------
    sc.settings.figdir = umap_dir

    sc.pl.umap(
        adata,
        color="leiden",
        save=f"_{sample_id}_umap.png",
        show=False
    )

    # ----------------------------
    # Save Spatial cluster plot
    # ----------------------------
    sc.settings.figdir = spatial_dir

    sc.pl.spatial(
        adata,
        color="leiden",
        spot_size=1.2,
        save=f"_{sample_id}_spatial_clusters.png",
        show=False
    )

    # ----------------------------
    # Save clustered object
    # ----------------------------
    clustered_path = results_dir / "processed"
    output_file = clustered_path / f"{sample_id}_clustered.h5ad"

    adata.write(output_file)

    print(f"   Clustering complete for {sample_id}")
    print(f"   Saved: {output_file}")

    return adata
