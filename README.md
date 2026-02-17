# Expression-Driven-Reconstruction-of-Human-DLPFC-Architecture-Using-Spatial-Transcriptomics

## Project Overview

This project investigates whether transcriptional profiles alone can reconstruct anatomical organization in human dorsolateral prefrontal cortex (DLPFC) tissue.

Using six independent 10x Genomics Visium sections from the spatialLIBD dataset, clustering was performed exclusively on gene expression features, intentionally excluding spatial coordinates during model construction.

Spatial information was used only for downstream validation.

## Research Question

Can intrinsic transcriptional similarity recover spatially coherent cortical domains without spatial priors?

## Dataset

- Platform: 10x Genomics Visium

- Tissue: Human DLPFC

- Sections: 151507, 151508, 151509, 151510, 151669, 151670

- Thousands of spatially indexed spots per section

- Each spot contains:

     - High-dimensional gene expression vector

     - Spatial coordinates (used for validation only)

## Computational Pipeline

The analysis was implemented as a modular Scanpy-based workflow.

## Per-Section Processing

For each section:

1. Quality control

   - Filtering based on total counts, gene counts, mitochondrial fraction

2. Library-size normalization

3. Log(1 + x) transformation

4. Highly variable gene selection (top 3000 genes)

5. Principal Component Analysis (PCA)

6. K-nearest neighbor graph construction

7. Leiden community detection

Clustering was performed purely in PCA space derived from transcriptional features.

Clusters were subsequently mapped to spatial coordinates to evaluate anatomical alignment.

## Cross-Sample Integration

To assess reproducibility across tissue sections:

- Batch-aware highly variable gene selection

- ComBat batch correction

- Integrated PCA embedding

- Re-clustering of merged dataset

Post-correction embeddings demonstrated strong sample mixing, indicating that recovered domains were not driven by batch effects.

## Spatial Quantification

Spatial structure was evaluated using Moran’s I to measure global spatial autocorrelation.

Observed maximum Moran’s I ≈ 0.07.

In layered cortical tissue, moderate positive autocorrelation reflects biologically realistic gradient-based organization rather than sharply bounded compartments.

This quantitative validation supports the hypothesis that expression-defined clusters correspond to anatomical structure.

## Key Outcomes

- Expression-only clustering recovers spatially coherent domains.

- Transcriptional domains are reproducible across independent sections.

- Batch correction preserves biological structure while removing technical variation.

- Spatial autocorrelation confirms non-random organization.

## Limitations & Extensions

- Clustering excluded spatial priors by design.

- Moran’s I provides global, not local, spatial insight.

- Parameter sensitivity (resolution, neighbors) could be systematically profiled.

- Spatially aware graph models represent a natural extension.

## Technical Stack

- Python

- Scanpy / AnnData

- scikit-learn

- ComBat (batch correction)

- Spatial statistics (Moran’s I)

## Repository Structure

```
scripts/        Modular pipeline components
notebooks/      Per-section and integration workflows
results/        Processed outputs and visualizations
```
Raw Visium data excluded due to size constraints.
