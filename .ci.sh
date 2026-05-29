#!/bin/bash
export PATH=/home/ksamuk/.local/bin:/home/ksamuk/miniconda3/envs/pixy-py311-test/bin:$PATH
cd /mnt/f/Dropbox/02_Projects/pixy
echo "=== ruff format --check ==="; ruff format --check 2>&1 | tail -5
echo "=== ruff check ==="; ruff check 2>&1 | tail -5
echo "=== mypy ==="; mypy 2>&1 | tail -8
echo "=== poetry check --lock ==="; poetry check --lock 2>&1 | tail -3
echo "=== pytest ==="; python -m pytest tests/ --no-cov -q 2>&1 | tail -8
