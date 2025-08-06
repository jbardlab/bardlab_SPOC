# SPOC Batch Analysis Script

## Overview

The `run_spoc_batch.sh` script provides an easy way to run SPOC analysis on multiple AlphaFold2/ColabFold prediction folders using Docker. It automatically discovers prediction folders containing PDB files and processes them individually.

## Usage

```bash
./scripts/run_spoc_batch.sh <input_folder> [output_folder] [options]
```

### Arguments

- `input_folder`: Directory containing AlphaFold2/ColabFold prediction folders
- `output_folder`: Directory for output CSV files (default: `input_folder/spoc_results`)

### Options

- `--single-complex`: Enable single complex mode (analyze single models by repeating 3x)
- `--verbose`: Show detailed output from SPOC analysis
- `--dry-run`: Show what would be processed without running analysis
- `--help`: Show help message

## Examples

### Basic usage
```bash
./scripts/run_spoc_batch.sh /path/to/predictions
```

### Specify custom output directory
```bash
./scripts/run_spoc_batch.sh /path/to/predictions /path/to/results
```

### Single complex mode with verbose output
```bash
./scripts/run_spoc_batch.sh /path/to/predictions --single-complex --verbose
```

### Dry run to see what would be processed
```bash
./scripts/run_spoc_batch.sh /path/to/predictions --dry-run
```

## Input Structure

The script expects your input directory to contain subdirectories with AlphaFold2/ColabFold prediction results. Each subdirectory should contain:

- `.pdb` files (structure predictions)
- `.json` files (confidence scores, PAE data, etc.)

Example structure:
```
predictions/
├── complex1/
│   ├── complex1_rank_001_model_2.pdb
│   ├── complex1_rank_002_model_4.pdb
│   ├── complex1_rank_003_model_1.pdb
│   ├── complex1_predicted_aligned_error_v1.json
│   └── complex1_scores_rank_001_model_2.json
├── complex2/
│   ├── complex2_rank_001_model_2.pdb
│   └── complex2_predicted_aligned_error_v1.json
└── complex3/
    └── ...
```

## Output

The script generates CSV files with comprehensive SPOC analysis results for each prediction folder:

- File naming: `{folder_name}_spoc_results.csv`
- Location: Specified output directory or `{input_folder}/spoc_results/`

### Output columns include:
- `complex_name`: Name of the analyzed complex
- `SPOC_nobio_score`: Main SPOC confidence score (0-1, higher is better)
- `single_complex_mode`: Whether single complex mode was used
- `avg_n_models`, `max_n_models`: Model statistics
- `iptm_mean`, `iptm_min`, `iptm_max`: Interface prediction confidence metrics
- `pdockq_e_max`, `pdockq_e_v2_max`: Docking quality scores
- `models_checked`: List of models that were analyzed
- `best_*`: Various metrics from the best-performing model
- `total_aa_length`: Total amino acid length
- `sequence`: Protein sequences

## Modes

### Standard Mode (Default)
- Uses all available models in each prediction folder
- Provides comprehensive analysis across multiple predictions
- Generally gives the most reliable SPOC scores

### Single Complex Mode (`--single-complex`)
- Uses only the best single model and repeats it 3 times
- Useful when you have folders with only one model
- Faster processing but potentially less comprehensive

## Performance

- Processing time varies based on complex size and number of models
- Typical range: 5-15 seconds per complex
- Docker overhead: ~2-3 seconds for image setup (one time)
- The script processes complexes sequentially

## Error Handling

- Skips folders that don't contain PDB files
- Skips analysis if output file already exists (delete to reprocess)
- Continues processing other folders if one fails
- Provides summary of successes and failures

## Docker Requirements

- Docker must be installed and running
- Script will attempt to pull `ghcr.io/jbardlab/spoc:v0.1` or build locally
- Uses `--platform linux/amd64` for ARM Mac compatibility (Rosetta)

## Troubleshooting

### No prediction folders found
- Check that your input directory contains subdirectories with `.pdb` files
- Verify the directory structure matches expected format

### Docker errors
- Ensure Docker is running
- Check that you have sufficient disk space for the Docker image
- On ARM Macs, ensure Rosetta is enabled for x86_64 emulation

### Analysis failures
- Use `--verbose` flag to see detailed error messages
- Check that prediction folders contain required files (PDB, JSON)
- Verify file permissions allow Docker to read the files

### Memory issues
- Large complexes may require more memory
- Adjust Docker memory limits if needed
- Consider processing smaller batches

## Integration

The script can be easily integrated into workflows:

```bash
# Process multiple batches
for batch in batch1 batch2 batch3; do
    ./scripts/run_spoc_batch.sh /data/$batch /results/$batch
done

# Combine results
cat /results/*/\*.csv > combined_results.csv
```
