#!/bin/bash
#SBATCH -J sample_processing_%a
#SBATCH -o sample_processing_%a.out.txt
#SBATCH -e sample_processing_%a.err.txt
#SBATCH --mem=200G
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zm77@duke.edu
#SBATCH --array=270

# Load required modules
module load Anaconda3
source ~/.bashrc
conda activate /hpc/group/adrc/zm77/software/conda_envs/zmsc
module load HDF5/1.14.3-rhel9

# Define directories and files
H5SEURAT_DIR="/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/Cell type objects"
MAPPING_FILE="/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/424sample.cell.csv"
BASE_DIR="/work/zm77/ROSMAP.MDD.2025312/scSampler"
SAMPLE_DIR="${BASE_DIR}/full_samples"
DOWNSAMPLED_DIR="${BASE_DIR}/downsampled_samples"
TEMP_DIR="${BASE_DIR}/temp"
SCRIPT_DIR="/hpc/group/adrc/zm77/LUTZ.MDD/ROSMAP.BIG.10032024/new.analysis.03112025/scripts"
SAMPLE_IDS_FILE="${BASE_DIR}/sample_ids.txt"

# Create directories
# mkdir -p "${SAMPLE_DIR}"
# mkdir -p "${DOWNSAMPLED_DIR}"
# mkdir -p "${TEMP_DIR}"

# Generate the sample IDs file if it doesn't exist
# if [ ! -f "${SAMPLE_IDS_FILE}" ]; then
  # echo "Creating sample IDs file..."
  # Rscript "${SCRIPT_DIR}/create_sample_ids.R" "${MAPPING_FILE}" "${SAMPLE_IDS_FILE}"
# fi

# Get the sample ID for this array job
SAMPLE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_IDS_FILE}")

echo "Processing sample: ${SAMPLE_ID}"

# Process this sample
Rscript "${SCRIPT_DIR}/process_single_sample.R" "${SAMPLE_ID}" "${H5SEURAT_DIR}" "${MAPPING_FILE}" "${SAMPLE_DIR}"

# Check if the sample was successfully processed
if [ -f "${SAMPLE_DIR}/${SAMPLE_ID}.rds" ]; then
  echo "Sample ${SAMPLE_ID} successfully processed"
  
  # Now run scSampler on this sample
  #echo "Running scSampler on sample ${SAMPLE_ID}"
  
  # Create a temporary directory for this sample
  # SAMPLE_TEMP_DIR="${TEMP_DIR}/${SAMPLE_ID}"
  # mkdir -p "${SAMPLE_TEMP_DIR}"
  
  # Run Python script for downsampling
  # /hpc/group/adrc/zm77/software/conda_envs/zmsc/bin/python "${SCRIPT_DIR}/subsample_single_seurat.py" \
    # --input_file "${SAMPLE_DIR}/${SAMPLE_ID}.rds" \
    # --output_dir "${DOWNSAMPLED_DIR}" \
    # --temp_dir "${SAMPLE_TEMP_DIR}" \
    # --fraction 0.05 \
    # --n_batch 16
    
  # Clean up temporary files
  #rm -rf "${SAMPLE_TEMP_DIR}"
  
  echo "Completed processing sample ${SAMPLE_ID}"
else
  echo "Error: Failed to process sample ${SAMPLE_ID}"
  exit 1
fi