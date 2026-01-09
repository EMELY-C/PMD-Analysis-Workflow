import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
import warnings
import os
import sys

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# ================= 1. Global Parameters =================

# !!! NOTE: Please update the file path below to match your local environment !!!
# Suggestion: Place the input file in the same directory as this script.
input_excel_file = "MSDIAL-PM-NEG.xlsx"

# Define output filenames
output_file_1 = '1_IS_Corrected_Imputed.xlsx'        # Preprocessed data
output_file_2 = '2_Final_TP_Analysis_By_Ion.xlsx'    # Paired results (Sheet by Mass Diff)
output_file_3 = '3_Matrix_Reshaped_Raw.xlsx'         # Intermediate: Reshaped matrix (Raw)
output_file_4 = '4_Matrix_RSD_Filtered.xlsx'         # Intermediate: RSD filtered matrix
output_file_5 = '5_Final_Correlation_Result.xlsx'    # Final Result: Correlation Report

# Ionization Mode
polar = 'NEG'  # Options: 'POS' or 'NEG'

# Mass Difference Settings (Transformation Mass Shifts)
if polar == 'POS':
    ions = [44.9851, 15.9949, 31.98982, 14.0157, 28.0313]
else:
    # Common transformations in NEG mode
    ions = [44.9851, 15.9949, 31.98982, 14.0157, 28.0313, 30.01057, 12.0000, 2.01565]

# Threshold Settings
MSMS_SIMILARITY_THRESHOLD = 0.4   # MS/MS cosine similarity > 0.4
DETECTION_RATE_THRESHOLD = 0.7    # Detection rate > 70%
RSD_THRESHOLD = 500               # RSD < 500% (Loosely filtered, can be adjusted to 30%)
CORRELATION_R_THRESHOLD = 0.6     # Spearman correlation r > 0.6
FDR_P_THRESHOLD = 0.05            # FDR adjusted p-value < 0.05


# ================= 2. Core Functions =================

def preprocess_msms(str_val, parent_mz, error=0.01):
    """
    Parse and clean MS/MS spectrum string.
    Removes precursor ion and potential noise.
    """
    if not isinstance(str_val, str): return np.array([])
    try:
        spectrum = []
        # Format expectation: "mz:intensity mz:intensity"
        for p in str_val.split(' '):
            if ':' in p:
                mz, inten = map(float, p.split(':'))
                spectrum.append([mz, inten])
        spectrum = np.array(spectrum)
    except:
        return np.array([])

    if len(spectrum) == 0: return spectrum

    # Filter 1: Remove ions larger than precursor or too close to precursor
    kept = []
    rejected_mzs = []
    for item in spectrum:
        mz = item[0]
        if (mz < (parent_mz - 10)) or (abs(mz - parent_mz) < error):
            kept.append(item)
        else:
            rejected_mzs.append(mz)

    kept = np.array(kept)
    if len(kept) == 0: return np.array([])

    # Filter 2: Noise removal based on specific neutral losses (e.g., CO2) if needed
    final_spectrum = []
    rejected_mzs = np.array(rejected_mzs)
    if len(rejected_mzs) > 0:
        for item in kept:
            mz = item[0]
            diffs = rejected_mzs - mz
            # Checking for common background noise diffs (approx 44 or 80 Da)
            is_noise = np.any((np.abs(diffs - 43.98983) < error) | (np.abs(diffs - 79.95682) < error))
            if not is_noise:
                final_spectrum.append(item)
    else:
        final_spectrum = kept

    final_spectrum = np.array(final_spectrum)
    if len(final_spectrum) == 0: return np.array([])

    # Normalize intensity and filter low abundance fragments
    max_inten = np.max(final_spectrum[:, 1])
    if max_inten > 0:
        final_spectrum[:, 1] /= max_inten
        # Keep fragments with relative intensity >= 5%
        final_spectrum = final_spectrum[final_spectrum[:, 1] >= 0.05]

    return final_spectrum


def calculate_similarity(str1, str2, mz1, mz2, error=0.01):
    """
    Calculate Cosine Similarity between two MS/MS spectra.
    """
    s1 = preprocess_msms(str1, mz1, error)
    s2 = preprocess_msms(str2, mz2, error)
    if len(s1) == 0 or len(s2) == 0: return 0.0

    # Match fragments within error tolerance
    count1 = sum(1 for item in s1 if np.min(np.abs(s2[:, 0] - item[0])) <= error)
    count2 = sum(1 for item in s2 if np.min(np.abs(s1[:, 0] - item[0])) <= error)

    # Calculate score
    return np.sqrt(count1 / len(s1)) * np.sqrt(count2 / len(s2))


def calculate_rsd(series):
    """
    Calculate Relative Standard Deviation (RSD%)
    Formula: (std / mean) * 100
    """
    mean = series.mean()
    if mean == 0: return 0
    return (series.std() / mean) * 100


def impute_row(row, cols):
    """
    Impute missing values (zeros) with min/sqrt(2) of the row.
    """
    vals = row[cols]
    non_zeros = vals[vals > 0]
    if len(non_zeros) > 0:
        return vals.replace(0, non_zeros.min() / np.sqrt(2))
    return vals


# ================= 3. Step 1: Data Preprocessing =================

print(f"========== Step 1: Data Loading & Preprocessing ==========")

# Check if input file exists
if not os.path.exists(input_excel_file):
    print(f"Error: Input file '{input_excel_file}' not found.")
    print("Please place the file in the directory or update the 'input_excel_file' path.")
    sys.exit()

try:
    # Sheet 0: Sample Data, Sheet 1: Internal Standard (IS) Data
    df_data = pd.read_excel(input_excel_file, sheet_name=0)
    df_is = pd.read_excel(input_excel_file, sheet_name=1)
except Exception as e:
    print(f"Error loading file: {e}")
    sys.exit()

# Identify Sample Columns
try:
    # Assuming columns after 'MS/MS spectrum' are sample columns
    msms_idx = df_data.columns.get_loc('MS/MS spectrum')
    sample_cols = df_data.columns[msms_idx + 1:]
    print(f"Detected {len(sample_cols)} sample columns.")
except:
    print("Error: Column 'MS/MS spectrum' not found in input data.")
    sys.exit()

# Internal Standard (IS) Normalization
print("Performing Internal Standard (IS) normalization...")
try:
    is_row = df_is.iloc[0][sample_cols].astype(float).replace(0, np.nan)
    df_corrected = df_data.copy()
    df_corrected[sample_cols] = df_corrected[sample_cols].div(is_row).fillna(0)
except Exception as e:
    print(f"Error during IS normalization: {e}")
    sys.exit()

# Filter by Detection Rate
print(f"Filtering features (Detection Rate > {DETECTION_RATE_THRESHOLD * 100}%)...")
df_corrected['Detection_Rate'] = (df_corrected[sample_cols] > 0).mean(axis=1)
df_filtered = df_corrected[df_corrected['Detection_Rate'] > DETECTION_RATE_THRESHOLD].copy()
print(f"Features remaining: {len(df_filtered)} / {len(df_data)}")

if df_filtered.empty:
    print("Error: No data remaining after filtration.")
    sys.exit()

# Missing Value Imputation
print("Imputing missing values...")
df_filtered[sample_cols] = df_filtered.apply(lambda row: impute_row(row, sample_cols), axis=1)

# Save Preprocessed Data
df_filtered.to_excel(output_file_1, index=False)
print(f"-> Preprocessed file saved: {output_file_1}")


# ================= 4. Step 2: TP Identification (Pairing) =================

print(f"\n========== Step 2: Transformation Product (TP) Search ==========")
mzs = df_filtered['Average Mz'].values
msms_list = df_filtered['MS/MS spectrum'].values

with pd.ExcelWriter(output_file_2) as writer:
    has_tp_data = False

    for ion in ions:
        current_ion_rows = []
        pair_counter = 1

        # Iterate through features to find pairs matching the mass difference
        for i in tqdm(range(len(mzs)), desc=f"Scanning Ion {ion:.4f}"):
            for j in range(i + 1, len(mzs)):
                diff = abs(mzs[i] - mzs[j])

                # Tolerance: 0.0005 Da (Strict PMD matching)
                if abs(diff - ion) < 0.0005:
                    if mzs[i] < mzs[j]:
                        idx_1, idx_2 = i, j
                    else:
                        idx_1, idx_2 = j, i

                    # Calculate MS/MS Similarity
                    sim = calculate_similarity(str(msms_list[idx_1]), str(msms_list[idx_2]), mzs[idx_1], mzs[idx_2])

                    if sim > MSMS_SIMILARITY_THRESHOLD:
                        row_1 = df_filtered.iloc[idx_1].copy()
                        row_2 = df_filtered.iloc[idx_2].copy()

                        # Annotate Pair Info
                        row_1['Pair_ID'] = f"pair{pair_counter}"
                        row_2['Pair_ID'] = f"pair{pair_counter}.1"
                        pmd_info = f"PMD ({ion}) Da"
                        row_1['PMD_Info'] = pmd_info
                        row_2['PMD_Info'] = pmd_info
                        row_1['Frag_Ratio_Score'] = sim
                        row_2['Frag_Ratio_Score'] = sim

                        current_ion_rows.append(row_1)
                        current_ion_rows.append(row_2)
                        pair_counter += 1

        if len(current_ion_rows) > 0:
            df_sheet = pd.DataFrame(current_ion_rows)
            # Reorder columns to put new info at the end
            base_cols = df_filtered.columns.tolist()
            new_cols = ['Pair_ID', 'PMD_Info', 'Frag_Ratio_Score']
            final_cols = [c for c in base_cols if c not in new_cols] + new_cols

            df_sheet = df_sheet[final_cols]
            # Save pairs to a specific sheet named after the mass difference
            df_sheet.to_excel(writer, sheet_name=str(ion), index=False)
            has_tp_data = True

    if not has_tp_data:
        print("No matching pairs found. Process terminated.")
        sys.exit()
    else:
        print(f"-> Paired results saved: {output_file_2}")


# ================= 5. Step 3: Statistical Analysis =================

print(f"\n========== Step 3: Matrix Reshaping, RSD Filtering & Correlation ==========")

try:
    all_sheets = pd.read_excel(output_file_2, sheet_name=None)
except Exception as e:
    print(f"Error loading Step 2 results: {e}")
    sys.exit()

# Initialize Excel Writers
writer_raw = pd.ExcelWriter(output_file_3)
writer_filtered = pd.ExcelWriter(output_file_4)
writer_final = pd.ExcelWriter(output_file_5)

has_final_res = False

for sheet_name, df in all_sheets.items():
    print(f"\nProcessing Sheet: {sheet_name} ...")

    # --- 1. Matrix Reshaping ---
    try:
        available_samples = [c for c in sample_cols if c in df.columns]
        # Extract Pair_ID and Sample Data
        df_matrix = df[['Pair_ID'] + available_samples].copy()
        # Transpose: Rows = Samples, Columns = Features (Pair_IDs)
        df_matrix = df_matrix.set_index('Pair_ID').T
        df_matrix.index.name = 'Sample_Name'
        df_matrix = df_matrix.apply(pd.to_numeric, errors='coerce')
    except Exception as e:
        print(f"  - Reshaping failed: {e}")
        continue

    # Save Raw Matrix
    df_matrix.to_excel(writer_raw, sheet_name=sheet_name)

    # --- 2. RSD Calculation & Filtering ---
    # Calculate RSD for each feature column
    rsd_series = df_matrix.apply(calculate_rsd, axis=0)

    # Keep features with RSD < Threshold
    valid_cols = rsd_series[rsd_series < RSD_THRESHOLD].index.tolist()
    df_matrix_filtered = df_matrix[valid_cols]

    # Save Filtered Matrix
    df_matrix_filtered.to_excel(writer_filtered, sheet_name=sheet_name)
    print(f"  - RSD Filter (<{RSD_THRESHOLD}%): Kept {len(valid_cols)} / {len(df_matrix.columns)} features")

    # --- 3. Spearman Correlation & FDR Correction ---
    unique_pairs = df['Pair_ID'].str.split('.').str[0].unique()
    results = []
    p_values = []

    valid_cols_set = set(df_matrix_filtered.columns)

    for base_pair in unique_pairs:
        col_1 = base_pair          # e.g., pair1
        col_2 = f"{base_pair}.1"   # e.g., pair1.1

        # Both features must pass RSD filter
        if (col_1 in valid_cols_set) and (col_2 in valid_cols_set):
            vec_1 = df_matrix_filtered[col_1]
            vec_2 = df_matrix_filtered[col_2]

            # Calculate Spearman Correlation
            if vec_1.std() == 0 or vec_2.std() == 0:
                rho, p = 0, 1.0
            else:
                rho, p = stats.spearmanr(vec_1, vec_2)

            # Retrieve PMD Info for reporting
            pmd_info = df[df['Pair_ID'] == col_1]['PMD_Info'].iloc[0] if 'PMD_Info' in df.columns else ''

            results.append({
                'Pair_Group': base_pair,
                'Feature_1': col_1,
                'Feature_2': col_2,
                'Spearman_Rho': rho,
                'P_Value': p,
                'PMD_Info': pmd_info
            })
            p_values.append(p)

    # FDR Correction (Benjamini-Hochberg)
    if results:
        # Handle potential NaNs in p-values
        p_values = [1.0 if np.isnan(x) else x for x in p_values]
        reject, fdr_p, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

        final_rows = []
        for i, res in enumerate(results):
            res['FDR_P_Value'] = fdr_p[i]

            # Final Filter Criteria: r > threshold AND FDR p < threshold
            if res['Spearman_Rho'] > CORRELATION_R_THRESHOLD and res['FDR_P_Value'] < FDR_P_THRESHOLD:
                final_rows.append(res)

        if final_rows:
            df_res = pd.DataFrame(final_rows)
            cols_order = ['Pair_Group', 'Spearman_Rho', 'FDR_P_Value', 'P_Value', 'PMD_Info', 'Feature_1', 'Feature_2']
            df_res = df_res[cols_order]

            df_res.to_excel(writer_final, sheet_name=sheet_name, index=False)
            has_final_res = True
            print(f"  - Final Result: Found {len(df_res)} significant pairs.")
        else:
            print("  - Final Result: No pairs met the correlation criteria.")

# Close and Save Writers
writer_raw.close()
writer_filtered.close()
writer_final.close()

print("\n========== Process Completed ==========")
print(f"1. Preprocessed Data:   {output_file_1}")
print(f"2. Paired List:         {output_file_2}")
print(f"3. Raw Matrix:          {output_file_3}")
print(f"4. RSD Filtered Matrix: {output_file_4}")
print(f"5. Final Report:        {output_file_5}")

if not has_final_res:
    print("\nNote: The final report (5) might be empty if no pairs passed all statistical filters.")
