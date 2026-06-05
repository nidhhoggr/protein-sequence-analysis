#!/bin/bash

# Initialize conda
source /root/miniconda3/etc/profile.d/conda.sh
conda activate bio

# Execute passed command or default shell
exec "$@"
