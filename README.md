# AI-Omics-for-Cancer-Regulatory-Gene-Identification

[![R](https://img.shields.io/badge/Language-R-blue)](#)
[![Workflow](https://img.shields.io/badge/Analysis-Transcriptomics%20%7C%20ML%20%7C%20PPI-green)](#)
[![Status](https://img.shields.io/badge/Project-Reproduction-orange)](#)

A reproduction-oriented workflow for identifying candidate regulatory genes in **aflatoxin B1 (AFB1)-induced hepatocellular carcinoma (HCC)** using transcriptomic analysis, target-database integration, co-expression network analysis, machine learning, and functional interpretation.

---

## Overview

This repository reproduces a published study on **AFB1-induced HCC** by integrating:

- multi-dataset transcriptomic preprocessing
- candidate target collection from public databases
- batch correction and differential expression analysis
- WGCNA module detection
- PPI network and enrichment analysis
- machine-learning-based feature selection and benchmarking
- final signature-model validation and interpretation

<p align="center">
  <img src="images/workflow_overview.png" alt="Workflow overview" width="900"><br>
  <em>Figure. Overall workflow of candidate-gene identification and final signature validation.</em>
</p>

> Replace the image path above with your own workflow figure after uploading it to an `images/` folder.

---

## Repository structure

```text
scripts/
├── 01_HCCdata_preparation.R
├── 02_target_database_integration.R
├── 03_batch_correction_and_DEG_analysis.R
├── 04_WGCNA_module_analysis_and_DEG_integration.R
├── 05_PPI_and_enrichment_analysis.R
├── 06_data_preparation_and_LassoRF.R
├── 07_machine_learning_model_comparison.R
└── 08_final_signature_model_validation_and_interpretation.R
