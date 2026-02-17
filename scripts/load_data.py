import scanpy as sc
from pathlib import Path


def load_sample(sample_id, data_dir):

    sample_path = data_dir / sample_id
    print(f"Loading data from: {sample_path}")

    adata = sc.read_visium(sample_path)

    adata.var_names_make_unique()

    return adata
