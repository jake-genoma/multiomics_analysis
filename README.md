# Multi-Omics Integrative Analysis Platform (OMNIAP)

![Shiny](https://img.shields.io/badge/Shiny-2.0-blue?logo=r&logoColor=white)
![License](https://img.shields.io/github/license/yourusername/multiomics-analysis-app)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXX)

An interactive Shiny application for systems biology analysis of transcriptomics, proteomics, and metabolomics datasets with AI-powered insights.

<p align="center">
  <img src="docs/screenshots/app-dashboard.png" width="800" alt="App Dashboard">
</p>

## Key Features

### 🔬 Multi-Omics Integration
- **Data Upload**: Supports CSV/TSV files for all omics layers
- **ID Conversion**: Automatic gene/protein/metabolite identifier mapping
- **Sample Matching**: Visual validation of sample consistency

### 📊 Advanced Analytics
- **Stratified Correlation**: Pearson/Spearman with metadata grouping
- **Differential Analysis**: limma (microarray) and DESeq2/edgeR (RNA-seq) pipelines
- **Pathway Enrichment**: KEGG, Reactome, and GO with interactive visualizations
- **Machine Learning**: Random Forest, SVM, XGBoost with cross-validation

### 🤖 AI-Powered Insights
- **Biomarker Discovery**: Integrated feature importance analysis
- **Hypothesis Generation**: LLM-driven mechanistic interpretations
- **Quality Control**: Automated data quality assessment

## Installation

### Prerequisites
- R (≥ 4.2.0)
- RStudio (recommended)

### Method 1: From GitHub
```r
# Install dependencies
install.packages("remotes")
remotes::install_github("yourusername/multiomics-analysis-app", dependencies = TRUE)

# Run the app
shiny::runGitHub("yourusername/multiomics-analysis-app")
