#!/usr/bin/env python3
"""
Reruns SPOC analysis using the top 3 models from individual analysis.
Assumes running inside the SPOC Docker container.
"""

import os
import sys
import argparse
import re
import glob
import shutil
import tempfile
import subprocess
import time
import pandas as pd
from pathlib import Path
from typing import List, Tuple, Optional


def copy_and_rename_model_files(prediction_folder: str, pdb_file: str, temp_dir: str, new_model_index: int) -> None:
    """Copy and rename relevant files for a single model to a temporary directory, preserving the filename prefix."""
    original_filename = os.path.basename(pdb_file)
    
    # --- Find the part of the filename to replace ---
    model_match = re.search(r'model_\d+', original_filename)
    if not model_match:
        print(f"    → Warning: Could not find 'model_X' pattern in {original_filename}. Skipping rename for this file.")
        # Fallback to simple copy if pattern not found
        shutil.copy2(os.path.join(prediction_folder, original_filename), temp_dir)
        return

    # --- Rename PDB file ---
    new_pdb_filename = re.sub(r'model_\d+', f'model_{new_model_index}', original_filename)
    shutil.copy2(os.path.join(prediction_folder, original_filename), os.path.join(temp_dir, new_pdb_filename))
    print(f"    → Copied and renamed PDB: {original_filename} → {new_pdb_filename}")

    # --- Find and rename corresponding JSON files ---
    # Extract original model and seed pattern from PDB filename for matching
    model_seed_match = re.search(r'(model_\d+_seed_\d+)', original_filename)
    if model_seed_match:
        model_seed_pattern = model_seed_match.group(1)
        
        # Find corresponding JSON files
        json_pattern = os.path.join(prediction_folder, f"*{model_seed_pattern}*.json")
        json_files = glob.glob(json_pattern)
        
        for json_file in json_files:
            original_json_filename = os.path.basename(json_file)
            # Create new JSON filename by replacing the model number part
            new_json_filename = re.sub(r'model_\d+', f'model_{new_model_index}', original_json_filename)
            shutil.copy2(json_file, os.path.join(temp_dir, new_json_filename))
            print(f"    → Copied and renamed JSON: {original_json_filename} → {new_json_filename}")
            
    else:
        # Fallback for filenames that don't match the expected pattern
        print(f"    → Warning: Could not extract model_seed pattern from {original_filename}. JSON files may not be copied correctly.")
        base_name = os.path.splitext(original_filename)[0]
        json_pattern = os.path.join(prediction_folder, f"*{base_name}*.json")
        for json_file in glob.glob(json_pattern):
            # Simple rename, might not be perfect
            new_json_filename = f"model_{new_model_index}_{os.path.basename(json_file)}"
            shutil.copy2(json_file, os.path.join(temp_dir, new_json_filename))
            print(f"    → Copied and renamed JSON (fallback): {os.path.basename(json_file)} → {new_json_filename}")


def run_spoc_analysis(temp_dir: str, output_file: str, rf_params: str, ipsae_script: str, verbose: bool = False) -> bool:
    """Run SPOC analysis on the temporary directory."""
    # Use bash to activate conda environment and run the script
    cmd = [
        "bash", "-c", 
        f"source /opt/conda/etc/profile.d/conda.sh && conda activate spoc_venv && python /repo/scripts/run_custom_nobio_v2.py {temp_dir} --rf_params {rf_params} --output {output_file} --ipsae_script {ipsae_script}"
    ]
    
    try:
        if verbose:
            print("    → Running SPOC analysis (verbose mode)...")
            subprocess.run(cmd, check=True, capture_output=False)
        else:
            print("    → Running SPOC analysis...")
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        return True
    except subprocess.CalledProcessError as e:
        print(f"    ✗ SPOC analysis failed with return code {e.returncode}")
        if verbose and e.stdout:
            print("STDOUT:", e.stdout)
        if verbose and e.stderr:
            print("STDERR:", e.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Rerun SPOC with the top 3 models from individual analysis."
    )
    
    parser.add_argument(
        "prediction_folder",
        help="Single AlphaFold2/ColabFold prediction folder containing PDBs and JSONs."
    )
    
    parser.add_argument(
        "results_csv",
        help="Path to the 'combined_individual_models_results.csv' file."
    )
    
    parser.add_argument(
        "output_folder",
        nargs='?',
        help="Directory for output CSV file (default: prediction_folder/top_models_analysis)"
    )
    
    parser.add_argument(
        "--rf_params",
        default="/repo/models/rf_afm_no_bio.joblib",
        help="Path to Random Forest parameters file."
    )
    
    parser.add_argument(
        "--ipsae_script",
        default="/repo/scripts/ipsae.py",
        help="Path to IPSAE script."
    )
    
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Show detailed output from SPOC analysis."
    )
    
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be processed without running analysis."
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.isdir(args.prediction_folder):
        print(f"Error: Prediction folder does not exist: {args.prediction_folder}")
        sys.exit(1)
        
    if not os.path.exists(args.results_csv):
        print(f"Error: Results CSV file not found: {args.results_csv}")
        sys.exit(1)
    
    # Set default output folder if not provided
    if not args.output_folder:
        args.output_folder = os.path.join(args.prediction_folder, "top_models_analysis")
    
    os.makedirs(args.output_folder, exist_ok=True)
    
    print("=== Top 3 Models SPOC Rerun ===")
    print(f"Prediction folder: {args.prediction_folder}")
    print(f"Results CSV: {args.results_csv}")
    print(f"Output folder: {args.output_folder}")
    print()

    # Read CSV and find top 3 models
    try:
        df = pd.read_csv(args.results_csv)
        if 'SPOC_nobio_score' not in df.columns or 'input_pdb_filename' not in df.columns:
            print("Error: CSV must contain 'SPOC_nobio_score' and 'input_pdb_filename' columns.")
            sys.exit(1)
            
        top_models = df.sort_values(by='SPOC_nobio_score', ascending=False).head(3)
        
        if len(top_models) < 3:
            print(f"Warning: Found only {len(top_models)} models in the results file.")
        
        if len(top_models) == 0:
            print("Error: No models found in the results file.")
            sys.exit(1)

        top_pdb_files = top_models['input_pdb_filename'].tolist()

    except Exception as e:
        print(f"Error reading or processing CSV file: {e}")
        sys.exit(1)

    print("Found top models to analyze:")
    for i, row in top_models.iterrows():
        print(f"  - PDB: {row['input_pdb_filename']}, SPOC Score: {row['SPOC_nobio_score']:.4f}")
    print()

    if args.dry_run:
        print("DRY RUN MODE - would copy the following files to a temporary directory and run SPOC:")
        for pdb_file in top_pdb_files:
            print(f"  - {pdb_file} and associated JSON files.")
        print()
        print("✓ Dry run completed!")
        return

    # Create temporary directory
    num_models = len(top_pdb_files)
    folder_name = os.path.basename(os.path.normpath(args.prediction_folder))
    with tempfile.TemporaryDirectory(prefix=f"spoc_top{num_models}_{folder_name}_") as temp_dir:
        print(f"Created temporary directory for {num_models} models: {temp_dir}")
        
        # Copy and rename top models' files
        for i, pdb_file in enumerate(top_pdb_files):
            copy_and_rename_model_files(args.prediction_folder, pdb_file, temp_dir, i + 1)
        
        # Run analysis
        output_file = os.path.join(args.output_folder, f"{folder_name}_top{num_models}_spoc_results.csv")
        
        print("\nStarting SPOC analysis on top models...")
        start_time = time.time()
        success = run_spoc_analysis(temp_dir, output_file, args.rf_params, args.ipsae_script, args.verbose)
        end_time = time.time()
        duration = int(end_time - start_time)

        # Summary
        print()
        print("=== Analysis Summary ===")
        if success and os.path.exists(output_file):
            print(f"✓ Top models analysis completed successfully in {duration}s.")
            print(f"Results saved to: {output_file}")
        else:
            print("✗ Top models analysis failed.")
            sys.exit(1)

if __name__ == "__main__":
    main()
