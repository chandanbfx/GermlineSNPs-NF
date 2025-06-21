#!/bin/bash
set -e

CONTAINER_DEF="containers/singularity.conda.def"
CONTAINER_SIF="containers/pipeline.sif"
NF_WORKFLOW="workflows/main.nf"

if [ ! -f "$CONTAINER_SIF" ]; then
  echo "Building Singularity container..."
  sudo singularity build "$CONTAINER_SIF" "$CONTAINER_DEF"
fi

echo "Starting pipeline..."
singularity exec --writable-tmpfs \
  --bind ${PWD}:/workspace:rw "$CONTAINER_SIF" \
  nextflow run "$NF_WORKFLOW" \
    -c "config/nextflow.config" \
    --samples "config/samples.csv" \
    --ref_genome "data/ref/chr1.fa" \
    -resume

echo "Pipeline completed!"
