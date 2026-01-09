# TP-Finder-PMD

**TP-Finder-PMD** is a Python-based workflow for identifying **Transformation Products (TPs)** in non-targeted metabolomics data. 

It automates the process of Paired Mass Difference (PMD) analysis, MS/MS spectrum similarity calculation, and statistical filtering (RSD, Spearman correlation, FDR correction) to discover significant precursor-TP pairs.

## Key Features

* **IS Normalization**: Automatically normalizes sample data using Internal Standards.
* **PMD Analysis**: Searches for specific mass shifts (e.g., +O, +2H, -CH2) indicative of biotransformation.
* **MS/MS Similarity**: Calculates cosine similarity between precursor and potential TPs to reduce false positives.
* **Statistical Filtering**:
    * **Detection Rate**: Filters out low-abundance features.
    * **RSD Filtering**: Removes unstable features based on Relative Standard Deviation.
    * **Correlation Analysis**: Spearman correlation with False Discovery Rate (FDR) control.

## Prerequisites

* Python 3.8 or higher
* Dependencies listed in `requirements.txt`

## Installation

1.  Clone this repository:
    ```bash
    git clone [https://github.com/EMELY-C/PMD-Analysis-Workflow.git]
    cd TP-Finder-PMD
    ```

2.  Install the required packages:
    ```bash
    pip install -r requirements.txt
    ```

## Data Preparation

Please refer to the provided **`demo_input.xlsx`** file included in this repository. 

This file demonstrates the required structure (including Sheet 1 for data and Sheet 2 for internal standards) and column formatting needed to run the script.

## Usage

1.  Place your data file in the project directory.
2.  Open the script `PMD-find-TPs.py` and ensure the `input_excel_file` variable matches your filename.
3.  Run the script:
    ```bash
    python PMD-find-TPs.py
    ```

## Output Files

The workflow generates the following files in order:

1.  **1_IS_Corrected_Imputed.xlsx**: Data after IS normalization, detection rate filtering, and missing value imputation.
2.  **2_Final_TP_Analysis_By_Ion.xlsx**: List of candidate TP pairs, separated by mass difference type (Sheets named by mass shift).
3.  **3_Matrix_Reshaped_Raw.xlsx**: Reshaped intensity matrix for paired features.
4.  **4_Matrix_RSD_Filtered.xlsx**: Data matrix after filtering by Relative Standard Deviation (RSD).
5.  **5_Final_Correlation_Result.xlsx**: **Final Report**. Contains TP pairs that passed all statistical thresholds.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
