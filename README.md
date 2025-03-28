# FlashRMSD
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15097621.svg)](https://doi.org/10.5281/zenodo.15097621)

**FlashRMSD** is a high-performance, symmetry-corrected RMSD calculation tool designed for speed, accuracy, and robustness. It leverages optimized atom featurization and matching algorithms to outperform existing solutions, making it ideal for molecular modeling, conformer analysis, and automated pipeline integration.

### 🔧 Key Features
- Ultra-fast execution with optimized backtracking and pruning
- Accurate symmetry handling for reliable RMSD calculations
- Lightweight and integrable with minimal dependencies

FlashRMSD is especially suited for high-throughput and precision-demanding applications where both performance and correctness are critical.

---

## 🚀 Getting Started

### 1. Build the Tool

To compile FlashRMSD:

```bash
make FlashRMSD
```

The compiled binary will be located in the `tools/` directory.

---

### 2. Usage

```bash
tools/FlashRMSD <query_file_path> [<template_file_path>] 
            [-n]   # Run naive version of algorithm
            [-x]   # Compute cross-RMSD across all conformations in query file
            [-h]   # Include hydrogen atoms
            [-b]   # Use bond types
            [-v]   # Verbose output (runtimes, comments)
            [-a]   # Print atom assignments
```

---

### 3. Examples

**Simple RMSD Calculation**

```bash
$ tools/FlashRMSD data/examples/60C_1.sdf data/examples/60C_2.sdf
0.225653
```

**Verbose Output**

```bash
$ tools/FlashRMSD data/examples/60C_1.sdf data/examples/60C_2.sdf -v
Conformer 0: 
RMSD: 0.225653
Time taken: 0.000519
```

**Cross Format Comparison**

```bash
$ tools/FlashRMSD data/examples/60C_1.sdf data/examples/60C_2.mol2
0.225653
```

**Multi-Conformer Query File**

```bash
$ tools/FlashRMSD data/examples/60C.sdf data/examples/60C_1.sdf
0.000000
0.225653
0.051647
0.223659
0.186011
0.183432
0.225263
1.971918
0.193740
```

**Cross-RMSD for All Conformations**

```bash
$ tools/FlashRMSD data/examples/60C.sdf -x
0.000000 0.225653 0.051647 0.223659 0.186011 0.183432 0.225263 1.971918 0.193740 
0.225653 0.000000 0.196025 0.115033 0.134188 0.109236 0.125729 1.995769 0.103963 
0.051647 0.196025 0.000000 0.189846 0.167663 0.160704 0.191554 1.970697 0.166611 
0.223659 0.115033 0.189846 0.000000 0.086300 0.078457 0.020948 1.949745 0.059418 
0.186011 0.134188 0.167663 0.086300 0.000000 0.062541 0.080091 1.962163 0.067065 
0.183432 0.109236 0.160704 0.078457 0.062541 0.000000 0.080833 1.990632 0.021347 
0.225263 0.125729 0.191554 0.020948 0.080091 0.080833 0.000000 1.948782 0.063350 
1.971918 1.995769 1.970697 1.949745 1.962163 1.990632 1.948782 0.000000 1.984511 
0.193740 0.103963 0.166611 0.059418 0.067065 0.021347 0.063350 1.984511 0.000000
```

---

## 🧪 Benchmarking (Optional)

This repository includes the source code of the [DockRMSD](https://doi.org/10.1186/s13321-019-0362-7) tool, along with a modified version (**DockRMSDExt**) used for benchmarking comparisons.

To build the benchmarking tools:

```bash
make DockRMSD
make DockRMSDExt
```

Or build all tools at once:

```bash
make all
```

> ⚠️ Note: DockRMSD's source was included for reproducibility. Please refer to the original project for license and reuse information.

---

## 📂 Benchmark Dataset

This repository is accompanied by a comprehensive benchmark dataset, designed to evaluate RMSD calculation tools across diverse and structurally challenging molecules.

🔗 [Benchmark Dataset on Zenodo](https://doi.org/10.5281/zenodo.15097621)  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15097621.svg)](https://doi.org/10.5281/zenodo.15097621)

The dataset includes:
- 45,000+ molecules from CCD and BIRD repositories
- Multi-conformer files (SDF and MOL2)
- Per-pose files for all-to-all benchmarking
- A `results.csv` file with ground truth cross-RMSD values (blank where values could not be computed)

> **Citation:**  
> If you use the dataset, please cite:  
> **DOI: [10.5281/zenodo.15097621](https://doi.org/10.5281/zenodo.15097621)**
