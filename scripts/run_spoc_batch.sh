#!/bin/bash
# Simple batch SPOC analysis script using Docker
# Usage: ./run_spoc_batch.sh <input_folder> [output_folder] [options]

set -e

# Script configuration
CONTAINER_IMAGE="ghcr.io/jbardlab/spoc:v0.1"
LOCAL_CONTAINER="spoc:v1"
SCRIPT_DIR="$(dirname "$0")"
REPO_DIR="$(dirname "$0")/.."

# Default values
OUTPUT_FOLDER=""
SINGLE_COMPLEX_MODE=false
VERBOSE=false
DRY_RUN=false

# Function to display usage
usage() {
    echo "Usage: $0 <input_folder> [output_folder] [options]"
    echo ""
    echo "Arguments:"
    echo "  input_folder    Directory containing AlphaFold2/ColabFold prediction folders"
    echo "  output_folder   Directory for output CSV files (default: input_folder/spoc_results)"
    echo ""
    echo "Options:"
    echo "  --single-complex    Enable single complex mode (analyze single models by repeating 3x)"
    echo "  --verbose          Show detailed output from SPOC analysis"
    echo "  --dry-run          Show what would be processed without running analysis"
    echo "  --help             Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 /path/to/predictions"
    echo "  $0 /path/to/predictions /path/to/results --single-complex"
    echo "  $0 /path/to/predictions --verbose --dry-run"
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

# Function to find prediction folders
find_prediction_folders() {
    local input_dir="$1"
    local folders=()
    
    # Look for folders containing PDB files
    while IFS= read -r -d '' folder; do
        if find "$folder" -maxdepth 1 -name "*.pdb" | grep -q .; then
            folders+=("$folder")
        fi
    done < <(find "$input_dir" -type d -print0)
    
    printf '%s\n' "${folders[@]}"
}

# Function to run SPOC analysis on a single folder
run_spoc_analysis() {
    local prediction_folder="$1"
    local output_folder="$2"
    local folder_name=$(basename "$prediction_folder")
    local output_file="$output_folder/${folder_name}_spoc_results.csv"
    
    # Convert to absolute paths for Docker
    local abs_prediction_folder=$(realpath "$prediction_folder")
    local abs_output_folder=$(realpath "$output_folder")
    local abs_repo_dir=$(realpath "$REPO_DIR")
    
    echo "  Processing: $folder_name"
    
    # Check if output already exists
    if [[ -f "$output_file" ]]; then
        echo "    ⚠ Output file already exists: $output_file"
        echo "    → Skipping (delete file to reprocess)"
        return 0
    fi
    
    # Build Docker command
    local docker_cmd="docker run --platform linux/amd64 --rm"
    docker_cmd+=" -v \"$abs_prediction_folder\":/input"
    docker_cmd+=" -v \"$abs_repo_dir\":/repo"
    docker_cmd+=" -v \"$abs_output_folder\":/output"
    docker_cmd+=" \"$LOCAL_CONTAINER\""
    docker_cmd+=" python /repo/scripts/run_custom_nobio_v2.py /input"
    docker_cmd+=" --rf_params /repo/models/rf_afm_no_bio.joblib"
    docker_cmd+=" --output \"/output/${folder_name}_spoc_results.csv\""
    docker_cmd+=" --ipsae_script /repo/scripts/ipsae.py"
    
    if [[ "$SINGLE_COMPLEX_MODE" == "true" ]]; then
        docker_cmd+=" --single_complex"
    fi
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "    → Would run: $docker_cmd"
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
    
    # Check if analysis was successful
    if [[ -f "$output_file" ]]; then
        # Extract SPOC score from output
        local spoc_score="N/A"
        if [[ -s "$output_file" ]]; then
            spoc_score=$(tail -n 1 "$output_file" | cut -d',' -f2 2>/dev/null || echo "N/A")
        fi
        echo "    ✓ Analysis completed in ${duration}s, SPOC score: $spoc_score"
    else
        echo "    ✗ Analysis failed - no output file created"
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

INPUT_FOLDER="$1"
shift

# Check if second argument is a directory (output folder) or an option
if [[ $# -gt 0 && ! "$1" =~ ^-- ]]; then
    OUTPUT_FOLDER="$1"
    shift
fi

# Parse options
while [[ $# -gt 0 ]]; do
    case $1 in
        --single-complex)
            SINGLE_COMPLEX_MODE=true
            shift
            ;;
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
if [[ ! -d "$INPUT_FOLDER" ]]; then
    echo "Error: Input folder does not exist: $INPUT_FOLDER"
    exit 1
fi

# Set default output folder if not provided
if [[ -z "$OUTPUT_FOLDER" ]]; then
    OUTPUT_FOLDER="$INPUT_FOLDER/spoc_results"
fi

# Create output folder
mkdir -p "$OUTPUT_FOLDER"

echo "=== SPOC Batch Analysis ==="
echo "Input folder: $INPUT_FOLDER"
echo "Output folder: $OUTPUT_FOLDER"
if [[ "$SINGLE_COMPLEX_MODE" == "true" ]]; then
    echo "Mode: Single complex (models repeated 3x)"
else
    echo "Mode: Standard (multiple models if available)"
fi
echo ""

# Setup Docker if not in dry-run mode
if [[ "$DRY_RUN" != "true" ]]; then
    setup_docker
fi

# Find all prediction folders
echo "Scanning for prediction folders..."
prediction_folders=($(find_prediction_folders "$INPUT_FOLDER"))

if [[ ${#prediction_folders[@]} -eq 0 ]]; then
    echo "Error: No prediction folders containing PDB files found in $INPUT_FOLDER"
    exit 1
fi

echo "Found ${#prediction_folders[@]} prediction folder(s) to process:"
for folder in "${prediction_folders[@]}"; do
    echo "  - $(basename "$folder")"
done
echo ""

# Process each prediction folder
success_count=0
failed_count=0
start_time=$(date +%s)

for folder in "${prediction_folders[@]}"; do
    if run_spoc_analysis "$folder" "$OUTPUT_FOLDER"; then
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
echo "Total folders processed: ${#prediction_folders[@]}"
echo "Successful analyses: $success_count"
echo "Failed analyses: $failed_count"
echo "Total time: ${total_duration}s"

if [[ "$DRY_RUN" != "true" ]]; then
    echo ""
    echo "Results saved to: $OUTPUT_FOLDER"
    echo "Output files:"
    ls -la "$OUTPUT_FOLDER"/*.csv 2>/dev/null | sed 's/^/  /'
fi

if [[ $failed_count -gt 0 ]]; then
    exit 1
else
    echo "✓ All analyses completed successfully!"
fi
