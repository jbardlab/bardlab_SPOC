#!/bin/bash
# Source the conda initialization script
source activate spoc_venv
# Execute the command passed to the container
exec "$@"