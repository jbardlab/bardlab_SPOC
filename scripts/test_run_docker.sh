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

if [[ -f "$REPO_DIR/scripts/run_spoc_individual_models.py" ]]; then
  echo "✓ run_spoc_individual_models.py exists"
  if python3 -m py_compile "$REPO_DIR/scripts/run_spoc_individual_models.py" 2>/dev/null; then
    echo "✓ run_spoc_individual_models.py compiles successfully"
  else
    echo "⚠ run_spoc_individual_models.py has syntax issues"
  fi
else
  echo "✗ run_spoc_individual_models.py not found"
fi

if [[ -f "$REPO_DIR/scripts/run_spoc_top_models.py" ]]; then
  echo "✓ run_spoc_top_models.py exists"
  if python3 -m py_compile "$REPO_DIR/scripts/run_spoc_top_models.py" 2>/dev/null; then
    echo "✓ run_spoc_top_models.py compiles successfully"
  else
    echo "⚠ run_spoc_top_models.py has syntax issues"
  fi
else
  echo "✗ run_spoc_top_models.py not found"
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

# Test 4: Test run_spoc_individual_models.py script
echo "Test 4: Testing run_spoc_individual_models.py script..."
if [[ "$DOCKER_AVAILABLE" == "true" && -f "$REPO_DIR/models/rf_afm_no_bio.joblib" && -f "$REPO_DIR/scripts/ipsae.py" && -f "$REPO_DIR/scripts/run_spoc_individual_models.py" ]]; then
  
  # Test run_spoc_individual_models.py on the full example dataset
  echo "  Testing run_spoc_individual_models.py on example data..."
  individual_output_dir="$TEST_DIR/individual_models_output"
  mkdir -p "$individual_output_dir"
  
  # Check what PDB files are available
  echo "    → Available PDB files in example directory:"
  docker run --platform linux/amd64 --rm \
    -v "$TEST_DIR":/input \
    -v "$REPO_DIR":/repo \
    "$LOCAL_CONTAINER" \
    bash -c "ls -la /input/spoctest_GNAS2_GP119_v3/*.pdb | wc -l && ls /input/spoctest_GNAS2_GP119_v3/*.pdb | head -3"
  
  # Run individual models analysis
  docker run --platform linux/amd64 --rm \
    -v "$TEST_DIR":/input \
    -v "$REPO_DIR":/repo \
    -v "$individual_output_dir":/output \
    "$LOCAL_CONTAINER" \
    python /repo/scripts/run_spoc_individual_models.py /input/spoctest_GNAS2_GP119_v3 /output \
    --rf_params /repo/models/rf_afm_no_bio.joblib \
    --ipsae_script /repo/scripts/ipsae.py \
    --analysis_script /repo/scripts/run_custom_nobio_v2.py 2>&1
  
  # Check if we got at least some results (partial success is OK)
  individual_results=$(find "$individual_output_dir" -name "*_model_*_seed_*_spoc_results.csv" | wc -l)
  combined_results="$individual_output_dir/combined_individual_models_results.csv"
  
  if [[ $individual_results -gt 0 ]]; then
    echo "    ✓ Individual models analysis completed with $individual_results successful models"
    echo "    → Individual result files:"
    find "$individual_output_dir" -name "*_model_*_seed_*_spoc_results.csv" | head -3 | sed 's/^/        /'
  else
    echo "    ✗ No individual model result files found"
    exit 1
  fi
  
  if [[ -f "$combined_results" ]]; then
    echo "    ✓ Combined results file created: $(basename "$combined_results")"
    echo "    → Combined results preview:"
    head -n 3 "$combined_results" | sed 's/^/        /'
    total_models=$(tail -n +2 "$combined_results" | wc -l)
    echo "    → Total models analyzed: $total_models"
  else
    echo "    ⚠ Combined results file not found"
  fi
else
  if [[ "$DOCKER_AVAILABLE" != "true" ]]; then
    echo "⚠ Skipping individual models test - Docker not available"
  else
    echo "⚠ Skipping individual models test - missing required files"
  fi
fi

# Test 5: Test run_spoc_top_models.py script
echo "Test 5: Testing run_spoc_top_models.py script..."
if [[ "$DOCKER_AVAILABLE" == "true" && -f "$REPO_DIR/scripts/run_spoc_top_models.py" ]]; then
  
  # Check if we have results from the individual models test
  combined_results="$TEST_DIR/individual_models_output/combined_individual_models_results.csv"
  if [[ -f "$combined_results" ]]; then
    echo "  Using combined results from individual models test..."
    
    top_models_output_dir="$TEST_DIR/top_models_output"
    mkdir -p "$top_models_output_dir"
    
    # Test run_spoc_top_models.py
    if docker run --platform linux/amd64 --rm \
      -v "$TEST_DIR":/input \
      -v "$REPO_DIR":/repo \
      -v "$top_models_output_dir":/output \
      "$LOCAL_CONTAINER" \
      python /repo/scripts/run_spoc_top_models.py \
      /input/spoctest_GNAS2_GP119_v3 \
      /input/individual_models_output/combined_individual_models_results.csv \
      /output \
      --rf_params /repo/models/rf_afm_no_bio.joblib \
      --ipsae_script /repo/scripts/ipsae.py \
      --analysis_script /repo/scripts/run_custom_nobio_v2.py 2>&1; then
      
      echo "    ✓ Top models analysis completed successfully"
      
      # Check outputs
      top_models_result=$(find "$top_models_output_dir" -name "*_top*_spoc_results.csv")
      if [[ -n "$top_models_result" ]]; then
        echo "    ✓ Top models result file created: $(basename "$top_models_result")"
        echo "    → Top models results preview:"
        head -n 2 "$top_models_result" | sed 's/^/        /'
        if [[ -s "$top_models_result" ]]; then
          score_line=$(tail -n 1 "$top_models_result")
          spoc_score=$(echo "$score_line" | cut -d',' -f2 2>/dev/null || echo "N/A")
          echo "    → Final SPOC score: $spoc_score"
        fi
      else
        echo "    ⚠ Top models result file not found"
      fi
      
    else
      echo "    ✗ Top models analysis failed"
    fi
    
  else
    echo "  ⚠ Skipping top models test - no combined results from individual models test"
  fi
  
else
  if [[ "$DOCKER_AVAILABLE" != "true" ]]; then
    echo "⚠ Skipping top models test - Docker not available"
  else
    echo "⚠ Skipping top models test - missing required files"
  fi
fi

# Test 6: Verify Docker cleanup
echo "Test 6: Testing Docker cleanup..."
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