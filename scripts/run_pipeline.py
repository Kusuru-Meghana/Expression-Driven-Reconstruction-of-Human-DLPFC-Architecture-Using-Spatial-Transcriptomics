# ============================================
# Main Spatial Transcriptomics Pipeline
# ============================================

from pathlib import Path

from load_data import load_sample
from qc import run_qc
from preprocess import run_preprocessing
from cluster import run_clustering
from markers import run_marker_analysis


# --------------------------------------------
# 1️⃣ Define project paths
# --------------------------------------------

# This assumes run_pipeline.py is inside /scripts
PROJECT_DIR = Path(__file__).resolve().parents[1]

DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "results"

QC_FIG_DIR = RESULTS_DIR / "figures" / "qc"

# Create results folder if it doesn't exist
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


# --------------------------------------------
# 2️⃣ Sample IDs (ALL SAMPLES)
# --------------------------------------------

SAMPLES = [
    "151507",
    "151508",
    "151509",
    "151510",
    "151669",
    "151670"
]


print("====================================")
print("Starting Spatial Transcriptomics Pipeline")
print("====================================")


# --------------------------------------------
# 3️⃣ Run pipeline for each sample
# --------------------------------------------

for sample_id in SAMPLES:

    print("\n------------------------------------")
    print(f"Processing sample: {sample_id}")
    print("------------------------------------")

    # ----------------------------
    # Load Data
    # ----------------------------
    adata = load_sample(sample_id, DATA_DIR)

    # ----------------------------
    # Quality Control
    # ----------------------------
    adata = run_qc(adata, sample_id, QC_FIG_DIR)

    # ----------------------------
    # Preprocessing
    # ----------------------------
    adata = run_preprocessing(adata, sample_id, RESULTS_DIR)

    # ----------------------------
    # Clustering
    # ----------------------------
    adata = run_clustering(adata, sample_id, RESULTS_DIR)

    # ----------------------------
    # Marker Analysis
    # ----------------------------
    adata = run_marker_analysis(adata, sample_id, RESULTS_DIR)

    print(f"Sample {sample_id} completed successfully.")


print("\n====================================")
print("Pipeline finished successfully.")
print("====================================")
