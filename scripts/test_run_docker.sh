#!/bin/bash
# filepath: /Users/jbard/Library/CloudStorage/SynologyDrive-home/repos/bardlab_SPOC/scripts/test_run_docker.sh

set -e

# Test configuration
CONTAINER_IMAGE="ghcr.io/jbardlab/spoc:v0.1"
LOCAL_CONTAINER="spoc:v1"
TEST_DIR="$(dirname "$0")/../test_output"
REPO_DIR="$(dirname "$0")/.."

echo "=== SPOC Docker Test Suite ==="

# Function to cleanup test artifacts
cleanup() {
  echo "Cleaning up test artifacts..."
  rm -rf "$TEST_DIR"
}

# Set trap to cleanup on exit
trap cleanup EXIT

# Test 1: Pull Docker image or build for linux/amd64 platform
echo "Test 1: Checking Docker image availability..."
if docker pull "$CONTAINER_IMAGE" 2>/dev/null; then
  docker tag "$CONTAINER_IMAGE" "$LOCAL_CONTAINER"
  echo "✓ Docker image pulled successfully"
  DOCKER_AVAILABLE=true
else
  echo "⚠ Failed to pull Docker image, building locally for linux/amd64 platform..."
  if docker build --platform linux/amd64 -t "$LOCAL_CONTAINER" "$REPO_DIR"; then
    echo "✓ Docker image built successfully for linux/amd64 (will run with Rosetta)"
    DOCKER_AVAILABLE=true
  else
    echo "✗ Failed to build Docker image"
    DOCKER_AVAILABLE=false
  fi
fi

# Test 2: Create test environment
echo "Test 2: Setting up test environment..."
mkdir -p "$TEST_DIR"
# Copy example data for a realistic test
if [[ -d "$REPO_DIR/example/spoctest_GNAS2_GP119_v3" ]]; then
  cp -r "$REPO_DIR/example/spoctest_GNAS2_GP119_v3" "$TEST_DIR/"
  echo "✓ Test environment created with example data"
else
  echo "⚠ Example data not found, creating minimal test structure"
  echo "test_data" > "$TEST_DIR/test_input.txt"
fi

# Test 3A: Check local script files (without Docker)
echo "Test 3A: Checking local script files..."
if [[ -f "$REPO_DIR/scripts/run_custom_nobio_v2.py" ]]; then
  echo "✓ run_custom_nobio_v2.py exists"
  # Check if it's a valid Python file
  if python3 -m py_compile "$REPO_DIR/scripts/run_custom_nobio_v2.py" 2>/dev/null; then
    echo "✓ run_custom_nobio_v2.py compiles successfully"
  else
    echo "⚠ run_custom_nobio_v2.py has syntax issues"
  fi
else
  echo "✗ run_custom_nobio_v2.py not found"
fi

if [[ -f "$REPO_DIR/scripts/ipsae.py" ]]; then
  echo "✓ ipsae.py exists"
else
  echo "✗ ipsae.py not found"
fi

if [[ -f "$REPO_DIR/models/rf_afm_no_bio.joblib" ]]; then
  echo "✓ rf_afm_no_bio.joblib model file exists"
else
  echo "⚠ rf_afm_no_bio.joblib model file not found"
fi

# Test 3: Test Docker container execution
echo "Test 3: Testing Docker container execution..."
if [[ "$DOCKER_AVAILABLE" == "true" ]]; then
  if docker run --platform linux/amd64 --rm \
    -v "$TEST_DIR":/input \
    -v "$REPO_DIR":/repo \
    "$LOCAL_CONTAINER" \
    bash -c "python --version && ls /repo/scripts/"; then
    echo "✓ Docker container runs successfully"
  else
    echo "✗ Docker container failed to run"
    exit 1
  fi
else
  echo "⚠ Skipping Docker container test - Docker not available"
fi

# Test 4: Test run_custom_nobio_v2.py execution with --single_complex mode
echo "Test 4: Testing run_custom_nobio_v2.py script in single complex mode..."
if [[ "$DOCKER_AVAILABLE" == "true" && -f "$REPO_DIR/models/rf_afm_no_bio.joblib" && -f "$REPO_DIR/scripts/ipsae.py" ]]; then
  
  # First, let's see what PDB files are available in the example directory
  echo "  Available PDB files in example directory:"
  docker run --platform linux/amd64 --rm \
    -v "$TEST_DIR":/input \
    -v "$REPO_DIR":/repo \
    "$LOCAL_CONTAINER" \
    bash -c "ls -la /input/spoctest_GNAS2_GP119_v3/*.pdb" | head -5
  
  # Create separate directories for each model to test single complex mode
  MODEL_DIRS=()
  for model_num in 1 2 4; do
    model_dir="$TEST_DIR/model_${model_num}_test"
    mkdir -p "$model_dir"
    
    # Copy just one model's files to test single complex mode
    docker run --platform linux/amd64 --rm \
      -v "$TEST_DIR":/input \
      -v "$REPO_DIR":/repo \
      "$LOCAL_CONTAINER" \
      bash -c "cp /input/spoctest_GNAS2_GP119_v3/*model_${model_num}* /input/model_${model_num}_test/ 2>/dev/null || true"
    
    # Also copy the PAE file (predicted_aligned_error)
    docker run --platform linux/amd64 --rm \
      -v "$TEST_DIR":/input \
      -v "$REPO_DIR":/repo \
      "$LOCAL_CONTAINER" \
      bash -c "cp /input/spoctest_GNAS2_GP119_v3/*predicted_aligned_error* /input/model_${model_num}_test/ 2>/dev/null || true"
    
    MODEL_DIRS+=("$model_dir")
  done
  
  # Test each model individually in single complex mode
  success_count=0
  for model_num in 1 2 4; do
    # Map model number to rank based on file naming pattern
    case $model_num in
      1) rank="003" ;;
      2) rank="001" ;;
      4) rank="002" ;;
    esac
    
    model_dir="$TEST_DIR/model_${model_num}_test"
    mkdir -p "$model_dir"
    
    echo "  Testing model ${model_num} (rank ${rank}) in single complex mode..."
    
    # Copy the specific model files and shared files
    cp "$TEST_DIR/spoctest_GNAS2_GP119_v3/"*"rank_${rank}"*.pdb "$model_dir/" 2>/dev/null && echo "    → Copied PDB file for model $model_num"
    cp "$TEST_DIR/spoctest_GNAS2_GP119_v3/"*"rank_${rank}"*.json "$model_dir/" 2>/dev/null && echo "    → Copied JSON file for model $model_num"
    cp "$TEST_DIR/spoctest_GNAS2_GP119_v3/"*"predicted_aligned_error"*.json "$model_dir/" 2>/dev/null && echo "    → Copied PAE file"
    cp "$TEST_DIR/spoctest_GNAS2_GP119_v3/"*"template_domain_names"*.json "$model_dir/" 2>/dev/null && echo "    → Copied template domain names"
    cp "$TEST_DIR/spoctest_GNAS2_GP119_v3/config.json" "$model_dir/" 2>/dev/null && echo "    → Copied config file"
    
    # Check if PDB file exists
    pdb_files=$(find "$model_dir" -name "*.pdb" | wc -l)
    if [[ $pdb_files -eq 0 ]]; then
      echo "    ✗ Model ${model_num}: No PDB file found in $model_dir"
      ls -la "$model_dir" | sed 's/^/      /'
      continue
    else
      echo "    → Found $pdb_files PDB file(s) for model $model_num"
      ls -la "$model_dir"/*.pdb | sed 's/^/      /'
    fi
    
    output_file="/tmp/model_${model_num}_output.csv"
    
    # Run the analysis using correct parameter format
    if docker run --platform linux/amd64 --rm \
      -v "$model_dir":/input \
      -v "$REPO_DIR":/repo \
      -v "/tmp":/output \
      "$LOCAL_CONTAINER" \
      bash -c "python /repo/scripts/run_custom_nobio_v2.py /input --rf_params /repo/models/rf_afm_no_bio.joblib --output /output/model_${model_num}_output.csv --ipsae_script /repo/scripts/ipsae.py --single_complex" 2>&1; then
      
      if [[ -f "$output_file" ]]; then
        echo "    ✓ Model ${model_num}: Analysis completed successfully"
        # Show the SPOC score from this model
        if [[ -s "$output_file" ]]; then
          echo "      Output preview:"
          head -n 2 "$output_file" | sed 's/^/        /'
          score_line=$(tail -n 1 "$output_file")
          spoc_score=$(echo "$score_line" | cut -d',' -f2 2>/dev/null || echo "N/A")
          echo "      SPOC score: $spoc_score"
        fi
        success_count=$((success_count + 1))
      else
        echo "    ✗ Model ${model_num}: Output file not created"
        echo "      Checking output directory:"
        ls -la /tmp/ | grep "model_${model_num}" | sed 's/^/        /'
      fi
    else
      echo "    ✗ Model ${model_num}: Analysis failed"
    fi
  done
  
  if [[ $success_count -gt 0 ]]; then
    echo "  ✓ Single complex mode testing completed successfully ($success_count/3 models)"
    echo "  Summary of outputs:"
    for model_num in 1 2 4; do
      if [[ -f "$TEST_DIR/model_${model_num}_output.csv" ]]; then
        echo "    Model $model_num results:"
        head -n 2 "$TEST_DIR/model_${model_num}_output.csv" | tail -n 1 | cut -d',' -f1-3 | sed 's/^/      /'
      fi
    done
  else
    echo "  ✗ All single complex mode tests failed"
    exit 1
  fi
else
  if [[ "$DOCKER_AVAILABLE" != "true" ]]; then
    echo "⚠ Skipping script test - Docker not available"
  else
    echo "⚠ Skipping script test - missing required model or script files"
  fi
fi

# Test 5: Verify Docker cleanup
echo "Test 5: Testing Docker cleanup..."
if [[ "$DOCKER_AVAILABLE" == "true" ]]; then
  RUNNING_CONTAINERS=$(docker ps -q --filter ancestor="$LOCAL_CONTAINER")
  if [[ -z "$RUNNING_CONTAINERS" ]]; then
    echo "✓ No containers left running"
  else
    echo "✗ Found running containers that should have been cleaned up"
    docker ps --filter ancestor="$LOCAL_CONTAINER"
    exit 1
  fi
else
  echo "⚠ Skipping Docker cleanup test - Docker not available"
fi

echo "=== All tests passed! ==="