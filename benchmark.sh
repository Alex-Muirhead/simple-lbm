#!/bin/bash

for i in examples/taylor_green/*; do
    export CURRENT_DIR=$i;
    echo "Calling sbatch at $CURRENT_DIR";
    sbatch ./goslurm.sh;
done