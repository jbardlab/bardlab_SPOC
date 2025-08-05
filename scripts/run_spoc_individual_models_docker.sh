#!/bin/bash
# Individual model SPOC analysis script using Docker
# Usage: ./run_spoc_individual_models.sh <prediction_folder> [output_folder]

set -e

# Script configuration
CONTAINER_IMAGE="ghcr.io/jbardlab/spoc:v0.1"
LOCAL_CONTAINER="spoc:v1"
SCRIPT_DIR="$(dirname "$0")"
REPO_DIR="$(dirname "$0")/.."

# Default values
OUTPUT_FOLDER=""
VERBOSE=false
DRY_RUN=false

# Function to display usage
usage() {
    echo "Usage: $0 <prediction_folder> [output_folder] [options]"
    echo ""
    echo "Arguments:"
    echo "  prediction_folder   Single AlphaFold2/ColabFold prediction folder"
    echo "  output_folder       Directory for output CSV files (default: prediction_folder/individual_models)"
    echo ""
    echo "Options:"
    echo "  --verbose          Show detailed output from SPOC analysis"
    echo "  --dry-run          Show what would be processed without running analysis"
    echo "  --help             Show this help message"
    echo ""
    echo "Description:"
    echo "  This script analyzes each model individually in single-complex mode,"
    echo "  giving you separate SPOC scores for each model in the prediction folder."
    echo ""
    echo "Examples:"
    echo "  $0 /path/to/single_prediction"
    echo "  $0 /path/to/single_prediction /path/to/results --verbose"
    exit 1
}

# Function to check if Docker is available and pull/build image
setup_docker() {
    if ! command -v docker &> /dev/null; then
        echo "Error: Docker is not installed or not in PATH"
        exit 1
    fi
    
    echo "Setting up Docker image..."
    if docker pull "$CONTAINER_IMAGE" 2>/dev/null; then
        docker tag "$CONTAINER_IMAGE" "$LOCAL_CONTAINER"
        echo "✓ Docker image pulled successfully"
    else
        echo "⚠ Failed to pull Docker image, building locally..."
        if docker build --platform linux/amd64 -t "$LOCAL_CONTAINER" "$REPO_DIR"; then
            echo "✓ Docker image built successfully"
        else
            echo "✗ Failed to build Docker image"
            exit 1
        fi
    fi
}

# Function to extract model info from filename
get_model_info() {
    local pdb_file="$1"
    local filename=$(basename "$pdb_file")
    
    # Extract model number (look for patterns like model_1, model_2, etc.)
    if [[ $filename =~ model_([0-9]+) ]]; then
        echo "${BASH_REMATCH[1]}"
    elif [[ $filename =~ rank_([0-9]+) ]]; then
        # Map rank to model number based on typical AlphaFold naming
        case "${BASH_REMATCH[1]}" in
            "001") echo "2" ;;
            "002") echo "4" ;;
            "003") echo "1" ;;
            *) echo "${BASH_REMATCH[1]}" ;;
        esac
    else
        echo "unknown"
    fi
}

# Function to extract seed info from filename
get_seed_info() {
    local pdb_file="$1"
    local filename=$(basename "$pdb_file")
    
    # Extract seed number (look for patterns like seed_000, seed_001, etc.)
    if [[ $filename =~ seed_([0-9]+) ]]; then
        echo "${BASH_REMATCH[1]}"
    else
        echo "unknown"
    fi
}

# Function to run SPOC analysis on a single model
run_single_model_analysis() {
    local prediction_folder="$1"
    local output_folder="$2"
    local pdb_file="$3"
    
    local filename=$(basename "$pdb_file")
    local model_num=$(get_model_info "$pdb_file")
    local seed_num=$(get_seed_info "$pdb_file")
    local folder_name=$(basename "$prediction_folder")
    
    # Create temporary directory for this model only
    local temp_model_dir="/tmp/spoc_individual_${folder_name}_model_${model_num}_seed_${seed_num}"
    mkdir -p "$temp_model_dir"
    
    echo "  Processing model $model_num (seed $seed_num): $filename"
    
    # Copy the PDB file
    cp "$pdb_file" "$temp_model_dir/"
    
    # Extract model and seed pattern from PDB filename for matching JSON files
    if [[ $filename =~ (model_[0-9]+_seed_[0-9]+) ]]; then
        local model_seed_pattern="${BASH_REMATCH[1]}"
        echo "    → Looking for files matching pattern: *${model_seed_pattern}*"
        
        # Copy corresponding JSON files with same model_X_seed_Y pattern
        cp "$prediction_folder"/*"$model_seed_pattern"*.json "$temp_model_dir/" 2>/dev/null || true
    else
        echo "    → Warning: Could not extract model_seed pattern from $filename"
        # Fallback: try to match based on the full base name
        local base_name=$(echo "$filename" | sed 's/\.pdb$//')
        cp "$prediction_folder"/*"$base_name"*.json "$temp_model_dir/" 2>/dev/null || true
    fi
    
    # Copy shared files (PAE, template domain names, config) - these are needed for all models
    cp "$prediction_folder"/*predicted_aligned_error*.json "$temp_model_dir/" 2>/dev/null || true
    cp "$prediction_folder"/*template_domain_names*.json "$temp_model_dir/" 2>/dev/null || true
    cp "$prediction_folder"/config.json "$temp_model_dir/" 2>/dev/null || true
    
    # Convert to absolute paths for Docker
    local abs_temp_dir=$(realpath "$temp_model_dir")
    local abs_output_folder=$(realpath "$output_folder")
    local abs_repo_dir=$(realpath "$REPO_DIR")
    
    local output_file="$output_folder/${folder_name}_model_${model_num}_seed_${seed_num}_spoc_results.csv"
    
    # Check if output already exists
    if [[ -f "$output_file" ]]; then
        echo "    ⚠ Output file already exists: $(basename "$output_file")"
        echo "    → Skipping (delete file to reprocess)"
        rm -rf "$temp_model_dir"
        return 0
    fi
    
    # Build Docker command
    local docker_cmd="docker run --platform linux/amd64 --rm"
    docker_cmd+=" -v \"$abs_temp_dir\":/input"
    docker_cmd+=" -v \"$abs_repo_dir\":/repo"
    docker_cmd+=" -v \"$abs_output_folder\":/output"
    docker_cmd+=" \"$LOCAL_CONTAINER\""
    docker_cmd+=" python /repo/scripts/run_custom_nobio_v2.py /input"
    docker_cmd+=" --rf_params /repo/models/rf_afm_no_bio.joblib"
    docker_cmd+=" --output \"/output/${folder_name}_model_${model_num}_seed_${seed_num}_spoc_results.csv\""
    docker_cmd+=" --ipsae_script /repo/scripts/ipsae.py"
    docker_cmd+=" --single_complex"
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "    → Would run: $docker_cmd"
        rm -rf "$temp_model_dir"
        return 0
    fi
    
    # Run the analysis
    local start_time=$(date +%s)
    if [[ "$VERBOSE" == "true" ]]; then
        echo "    → Running SPOC analysis (verbose mode)..."
        eval "$docker_cmd"
    else
        echo "    → Running SPOC analysis..."
        eval "$docker_cmd" >/dev/null 2>&1
    fi
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    # Cleanup temp directory
    rm -rf "$temp_model_dir"
    
    # Check if analysis was successful
    if [[ -f "$output_file" ]]; then
        # Extract SPOC score from output
        local spoc_score="N/A"
        if [[ -s "$output_file" ]]; then
            spoc_score=$(tail -n 1 "$output_file" | cut -d',' -f2 2>/dev/null || echo "N/A")
        fi
        echo "    ✓ Model $model_num (seed $seed_num) completed in ${duration}s, SPOC score: $spoc_score"
        return 0
    else
        echo "    ✗ Model $model_num (seed $seed_num) analysis failed - no output file created"
        return 1
    fi
}

# Parse command line arguments
if [[ $# -eq 0 ]]; then
    usage
fi

# Check for help first
if [[ "$1" == "--help" || "$1" == "-h" ]]; then
    usage
fi

PREDICTION_FOLDER="$1"
shift

# Check if second argument is a directory (output folder) or an option
if [[ $# -gt 0 && ! "$1" =~ ^-- ]]; then
    OUTPUT_FOLDER="$1"
    shift
fi

# Parse options
while [[ $# -gt 0 ]]; do
    case $1 in
        --verbose)
            VERBOSE=true
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --help|-h)
            usage
            ;;
        *)
            echo "Error: Unknown option $1"
            usage
            ;;
    esac
done

# Validate input folder
if [[ ! -d "$PREDICTION_FOLDER" ]]; then
    echo "Error: Prediction folder does not exist: $PREDICTION_FOLDER"
    exit 1
fi

# Set default output folder if not provided
if [[ -z "$OUTPUT_FOLDER" ]]; then
    OUTPUT_FOLDER="$PREDICTION_FOLDER/individual_models"
fi

# Create output folder
mkdir -p "$OUTPUT_FOLDER"

# Find PDB files
pdb_files=($(find "$PREDICTION_FOLDER" -maxdepth 1 -name "*.pdb"))

if [[ ${#pdb_files[@]} -eq 0 ]]; then
    echo "Error: No PDB files found in $PREDICTION_FOLDER"
    exit 1
fi

echo "=== Individual Model SPOC Analysis ==="
echo "Prediction folder: $PREDICTION_FOLDER"
echo "Output folder: $OUTPUT_FOLDER"
echo "Found ${#pdb_files[@]} model(s) to analyze individually"
echo ""

# Setup Docker if not in dry-run mode
if [[ "$DRY_RUN" != "true" ]]; then
    setup_docker
fi

# Process each model individually
success_count=0
failed_count=0
start_time=$(date +%s)

for pdb_file in "${pdb_files[@]}"; do
    if run_single_model_analysis "$PREDICTION_FOLDER" "$OUTPUT_FOLDER" "$pdb_file"; then
        ((success_count++))
    else
        ((failed_count++))
    fi
done

end_time=$(date +%s)
total_duration=$((end_time - start_time))

# Summary
echo ""
echo "=== Analysis Summary ==="
echo "Total models processed: ${#pdb_files[@]}"
echo "Successful analyses: $success_count"
echo "Failed analyses: $failed_count"
echo "Total time: ${total_duration}s"

if [[ "$DRY_RUN" != "true" ]]; then
    echo ""
    echo "Results saved to: $OUTPUT_FOLDER"
    echo "Individual model results:"
    
    # Show summary table of results
    echo ""
    echo "Model | Seed | SPOC Score | File"
    echo "------|------|------------|-----"
    
    for csv_file in "$OUTPUT_FOLDER"/*.csv; do
        if [[ -f "$csv_file" ]]; then
            filename=$(basename "$csv_file")
            if [[ $filename =~ _model_([0-9]+)_seed_([0-9]+)_spoc_results ]]; then
                model_num="${BASH_REMATCH[1]}"
                seed_num="${BASH_REMATCH[2]}"
            else
                model_num="?"
                seed_num="?"
            fi
            spoc_score=$(tail -n 1 "$csv_file" | cut -d',' -f2 2>/dev/null || echo "N/A")
            printf "%-5s | %-4s | %-10s | %s\n" "$model_num" "$seed_num" "$spoc_score" "$filename"
        fi
    done
    
    echo ""
    echo "Creating combined results file..."
    
    # Create combined CSV file
    combined_file="$OUTPUT_FOLDER/combined_individual_models_results.csv"
    
    # First, create header by adding new columns to the original header
    first_csv=$(find "$OUTPUT_FOLDER" -name "*_model_*_seed_*_spoc_results.csv" | head -1)
    if [[ -f "$first_csv" ]]; then
        # Get original header and add new columns at the beginning
        original_header=$(head -n 1 "$first_csv")
        echo "model_number,seed_number,input_pdb_filename,$original_header" > "$combined_file"
        
        # Process each individual results file
        for csv_file in "$OUTPUT_FOLDER"/*_model_*_seed_*_spoc_results.csv; do
            if [[ -f "$csv_file" ]]; then
                filename=$(basename "$csv_file")
                
                # Extract model and seed numbers
                if [[ $filename =~ _model_([0-9]+)_seed_([0-9]+)_spoc_results ]]; then
                    model_num="${BASH_REMATCH[1]}"
                    seed_num="${BASH_REMATCH[2]}"
                else
                    model_num="unknown"
                    seed_num="unknown"
                fi
                
                # Find corresponding PDB filename by looking for the original PDB file
                pdb_filename="unknown"
                for pdb_file in "$PREDICTION_FOLDER"/*.pdb; do
                    if [[ -f "$pdb_file" ]]; then
                        pdb_basename=$(basename "$pdb_file")
                        if [[ $pdb_basename =~ model_${model_num}_seed_${seed_num} ]]; then
                            pdb_filename="$pdb_basename"
                            break
                        fi
                    fi
                done
                
                # Get the data row (skip header)
                data_row=$(tail -n 1 "$csv_file")
                
                # Combine new columns with original data
                echo "$model_num,$seed_num,$pdb_filename,$data_row" >> "$combined_file"
            fi
        done
        
        echo "✓ Combined results saved to: $combined_file"
        echo ""
        echo "Combined file contains $(tail -n +2 "$combined_file" | wc -l) model results with columns:"
        echo "  - model_number: AlphaFold model number"
        echo "  - seed_number: Random seed used"
        echo "  - input_pdb_filename: Original PDB file name"
        echo "  - All original SPOC analysis columns..."
    else
        echo "⚠ No individual model results found to combine"
    fi
fi

if [[ $failed_count -gt 0 ]]; then
    exit 1
else
    echo ""
    echo "✓ All individual model analyses completed successfully!"
fi
