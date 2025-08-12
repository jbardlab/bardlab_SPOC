#!/usr/bin/env python3
"""
Individual model SPOC analysis script (Python version)
Analyzes each AlphaFold2/ColabFold model individually in single-complex mode.
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


def get_model_info(pdb_file: str) -> str:
    """Extract model number from PDB filename."""
    filename = os.path.basename(pdb_file)
    
    # Extract model number (look for patterns like model_1, model_2, etc.)
    model_match = re.search(r'model_(\d+)', filename)
    if model_match:
        return model_match.group(1)
    
    # Map rank to model number based on typical AlphaFold naming
    rank_match = re.search(r'rank_(\d+)', filename)
    if rank_match:
        rank_map = {"001": "2", "002": "4", "003": "1"}
        return rank_map.get(rank_match.group(1), rank_match.group(1))
    
    return "unknown"


def get_seed_info(pdb_file: str) -> str:
    """Extract seed number from PDB filename."""
    filename = os.path.basename(pdb_file)
    
    # Extract seed number (look for patterns like seed_000, seed_001, etc.)
    seed_match = re.search(r'seed_(\d+)', filename)
    if seed_match:
        return seed_match.group(1)
    
    return "unknown"


def find_pdb_files(prediction_folder: str) -> List[str]:
    """Find all PDB files in the prediction folder."""
    pdb_pattern = os.path.join(prediction_folder, "*.pdb")
    return glob.glob(pdb_pattern)


def copy_model_files(prediction_folder: str, pdb_file: str, temp_dir: str) -> None:
    """Copy relevant files for a single model to temporary directory."""
    filename = os.path.basename(pdb_file)
    
    # Copy the PDB file
    shutil.copy2(pdb_file, temp_dir)
    
    # Extract model and seed pattern from PDB filename for matching JSON files
    model_seed_match = re.search(r'(model_\d+_seed_\d+)', filename)
    if model_seed_match:
        model_seed_pattern = model_seed_match.group(1)
        print(f"    → Looking for files matching pattern: *{model_seed_pattern}*")
        
        # Copy corresponding JSON files with same model_X_seed_Y pattern
        pattern = os.path.join(prediction_folder, f"*{model_seed_pattern}*.json")
        for json_file in glob.glob(pattern):
            print(f"    → Copying JSON file: {json_file}")
            shutil.copy2(json_file, temp_dir)
    else:
        print(f"    → Warning: Could not extract model_seed pattern from {filename}")
        # Fallback: try to match based on the full base name
        base_name = os.path.splitext(filename)[0]
        pattern = os.path.join(prediction_folder, f"*{base_name}*.json")
        for json_file in glob.glob(pattern):
            shutil.copy2(json_file, temp_dir)

def run_spoc_analysis(temp_dir: str, output_file: str, rf_params: str, ipsae_script: str, analysis_script: str, verbose: bool = False) -> bool:
    """Run SPOC analysis on the temporary directory."""
    # Use bash to activate conda environment and run the script
    cmd = [
        "bash", "-c", 
        f"python {analysis_script} {temp_dir} --rf_params {rf_params} --output {output_file} --ipsae_script {ipsae_script} --single_complex"
    ]
    
    try:
        if verbose:
            print("    → Running SPOC analysis (verbose mode)...")
            result = subprocess.run(cmd, check=True, capture_output=False)
        else:
            print("    → Running SPOC analysis...")
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        return True
    except subprocess.CalledProcessError as e:
        print(f"    ✗ SPOC analysis failed with return code {e.returncode}")
        if verbose and e.stdout:
            print("STDOUT:", e.stdout)
        if verbose and e.stderr:
            print("STDERR:", e.stderr)
        return False


def get_spoc_score(csv_file: str) -> str:
    """Extract SPOC score from output CSV file."""
    try:
        if os.path.exists(csv_file) and os.path.getsize(csv_file) > 0:
            df = pd.read_csv(csv_file)
            if len(df) > 0 and 'SPOC_nobio_score' in df.columns:
                return str(df['SPOC_nobio_score'].iloc[0])
    except Exception:
        pass
    return "N/A"


def analyze_single_model(prediction_folder: str, output_folder: str, pdb_file: str, 
                        rf_params: str, ipsae_script: str, analysis_script: str, verbose: bool = False) -> bool:
    """Analyze a single model in single-complex mode."""
    filename = os.path.basename(pdb_file)
    model_num = get_model_info(pdb_file)
    seed_num = get_seed_info(pdb_file)
    folder_name = os.path.basename(prediction_folder)
    
    print(f"  Processing model {model_num} (seed {seed_num}): {filename}")
    
    # Create output file path
    output_file = os.path.join(output_folder, f"{folder_name}_model_{model_num}_seed_{seed_num}_spoc_results.csv")
    
    # Check if output already exists
    if os.path.exists(output_file):
        print(f"    ⚠ Output file already exists: {os.path.basename(output_file)}")
        print("    → Skipping (delete file to reprocess)")
        return True
    
    # Create temporary directory for this model only
    with tempfile.TemporaryDirectory(prefix=f"spoc_individual_{folder_name}_model_{model_num}_seed_{seed_num}_") as temp_dir:
        try:
            # Copy relevant files to temp directory
            copy_model_files(prediction_folder, pdb_file, temp_dir)
            
            # Run SPOC analysis
            start_time = time.time()
            success = run_spoc_analysis(temp_dir, output_file, rf_params, ipsae_script, analysis_script, verbose)
            end_time = time.time()
            duration = int(end_time - start_time)
            
            if success and os.path.exists(output_file):
                spoc_score = get_spoc_score(output_file)
                print(f"    ✓ Model {model_num} (seed {seed_num}) completed in {duration}s, SPOC score: {spoc_score}")
                return True
            else:
                print(f"    ✗ Model {model_num} (seed {seed_num}) analysis failed - no output file created")
                return False
                
        except Exception as e:
            print(f"    ✗ Model {model_num} (seed {seed_num}) analysis failed with error: {e}")
            return False


def create_combined_csv(output_folder: str, prediction_folder: str) -> str:
    """Create combined CSV file with all individual model results."""
    print("Creating combined results file...")
    
    combined_file = os.path.join(output_folder, "combined_individual_models_results.csv")
    
    # Find all individual result files
    pattern = os.path.join(output_folder, "*_model_*_seed_*_spoc_results.csv")
    result_files = glob.glob(pattern)
    
    if not result_files:
        print("⚠ No individual model results found to combine")
        return ""
    
    # Get header from first file and add new columns
    first_file = result_files[0]
    df_first = pd.read_csv(first_file)
    original_columns = df_first.columns.tolist()
    new_columns = ['model_number', 'seed_number', 'input_pdb_filename'] + original_columns
    
    # Create combined dataframe
    combined_data = []
    
    for csv_file in result_files:
        filename = os.path.basename(csv_file)
        
        # Extract model and seed numbers from filename
        model_match = re.search(r'_model_(\d+)_seed_(\d+)_spoc_results', filename)
        if model_match:
            model_num = model_match.group(1)
            seed_num = model_match.group(2)
        else:
            model_num = "unknown"
            seed_num = "unknown"
        
        # Find corresponding PDB filename
        pdb_filename = "unknown"
        pdb_pattern = os.path.join(prediction_folder, f"*model_{model_num}_seed_{seed_num}*.pdb")
        pdb_files = glob.glob(pdb_pattern)
        if pdb_files:
            pdb_filename = os.path.basename(pdb_files[0])
        
        # Read the data
        try:
            df = pd.read_csv(csv_file)
            if len(df) > 0:
                # Add new columns to the data
                row_data = [model_num, seed_num, pdb_filename] + df.iloc[0].tolist()
                combined_data.append(row_data)
        except Exception as e:
            print(f"    ⚠ Error reading {csv_file}: {e}")
    
    # Create and save combined dataframe
    if combined_data:
        combined_df = pd.DataFrame(combined_data, columns=new_columns)
        combined_df.to_csv(combined_file, index=False)
        
        print(f"✓ Combined results saved to: {combined_file}")
        print(f"Combined file contains {len(combined_df)} model results with columns:")
        print("  - model_number: AlphaFold model number")
        print("  - seed_number: Random seed used")
        print("  - input_pdb_filename: Original PDB file name")
        print("  - All original SPOC analysis columns...")
        
        return combined_file
    else:
        print("⚠ No valid data found to combine")
        return ""


def show_summary_table(output_folder: str) -> None:
    """Show summary table of individual model results."""
    print("Individual model results:")
    print()
    print("Model | Seed | SPOC Score | File")
    print("------|------|------------|-----")
    
    pattern = os.path.join(output_folder, "*_model_*_seed_*_spoc_results.csv")
    result_files = sorted(glob.glob(pattern))
    
    for csv_file in result_files:
        filename = os.path.basename(csv_file)
        
        # Extract model and seed
        model_match = re.search(r'_model_(\d+)_seed_(\d+)_spoc_results', filename)
        if model_match:
            model_num = model_match.group(1)
            seed_num = model_match.group(2)
        else:
            model_num = "?"
            seed_num = "?"
        
        spoc_score = get_spoc_score(csv_file)
        print(f"{model_num:<5} | {seed_num:<4} | {spoc_score:<10} | {filename}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze each AlphaFold2/ColabFold model individually in single-complex mode"
    )
    
    parser.add_argument(
        "prediction_folder",
        help="Single AlphaFold2/ColabFold prediction folder"
    )
    
    parser.add_argument(
        "output_folder",
        nargs='?',
        help="Directory for output CSV files (default: prediction_folder/individual_models)"
    )
    
    parser.add_argument(
        "--rf_params",
        default="/repo/models/rf_afm_no_bio.joblib",
        help="Path to Random Forest parameters file"
    )
    
    parser.add_argument(
        "--ipsae_script",
        default="/repo/scripts/ipsae.py",
        help="Path to IPSAE script"
    )
    
    parser.add_argument(
        "--analysis_script",
        default="/repo/scripts/run_custom_nobio_v2.py",
        help="Path to SPOC analysis script (run_custom_nobio_v2.py)"
    )
    
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Show detailed output from SPOC analysis"
    )
    
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be processed without running analysis"
    )
    
    args = parser.parse_args()
    
    # Validate input folder
    if not os.path.isdir(args.prediction_folder):
        print(f"Error: Prediction folder does not exist: {args.prediction_folder}")
        sys.exit(1)
    
    # Set default output folder if not provided
    if not args.output_folder:
        args.output_folder = os.path.join(args.prediction_folder, "individual_models")
    
    # Create output folder
    os.makedirs(args.output_folder, exist_ok=True)
    
    # Find PDB files
    pdb_files = find_pdb_files(args.prediction_folder)
    
    if not pdb_files:
        print(f"Error: No PDB files found in {args.prediction_folder}")
        sys.exit(1)
    
    print("=== Individual Model SPOC Analysis ===")
    print(f"Prediction folder: {args.prediction_folder}")
    print(f"Output folder: {args.output_folder}")
    print(f"Found {len(pdb_files)} model(s) to analyze individually")
    print()
    
    if args.dry_run:
        print("DRY RUN MODE - showing what would be processed:")
        for pdb_file in pdb_files:
            filename = os.path.basename(pdb_file)
            model_num = get_model_info(pdb_file)
            seed_num = get_seed_info(pdb_file)
            print(f"  Would process model {model_num} (seed {seed_num}): {filename}")
        print()
        print("=== Analysis Summary ===")
        print(f"Total models that would be processed: {len(pdb_files)}")
        print("✓ Dry run completed!")
        return
    
    # Process each model individually
    success_count = 0
    failed_count = 0
    start_time = time.time()
    
    for pdb_file in pdb_files:
        if analyze_single_model(args.prediction_folder, args.output_folder, pdb_file, 
                               args.rf_params, args.ipsae_script, args.analysis_script, args.verbose):
            success_count += 1
        else:
            failed_count += 1
    
    end_time = time.time()
    total_duration = int(end_time - start_time)
    
    # Summary
    print()
    print("=== Analysis Summary ===")
    print(f"Total models processed: {len(pdb_files)}")
    print(f"Successful analyses: {success_count}")
    print(f"Failed analyses: {failed_count}")
    print(f"Total time: {total_duration}s")
    
    if success_count > 0:
        print()
        print(f"Results saved to: {args.output_folder}")
        
        # Show summary table
        show_summary_table(args.output_folder)
        
        print()
        
        # Create combined CSV
        combined_file = create_combined_csv(args.output_folder, args.prediction_folder)
        
        if combined_file:
            print()
    
    if failed_count > 0:
        if success_count > 0:
            print(f"⚠ Completed with {failed_count} failed analyses out of {len(pdb_files)} total")
            print("✓ Partial success - at least some individual model analyses completed!")
        else:
            print("✗ All analyses failed!")
            sys.exit(1)
    else:
        print("✓ All individual model analyses completed successfully!")


if __name__ == "__main__":
    main()
