#!/bin/bash

#SBATCH -A csc359
#SBATCH -p batch

#SBATCH -t 10
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -n 1
#SBATCH --gpu-bind=closest

#SBATCH -J exasgd_opflow_profiling
#SBATCH -o profile_opflow_frontier_batch.%J
#SBATCH -e profile_opflow_frontier_batch.%J

#SBATCH --mail-user=sayefsakin@gmail.com
#SBATCH --mail-type=ALL

source buildsystem/clang-hip/frontier/base.sh
source buildsystem/spack/frontier/modules/dependencies.sh
source buildsystem/spack/frontier/modules/exago.sh

python3 performance_analysis/perf_pipeline.py performance_analysis/opflow_testsuite.toml
