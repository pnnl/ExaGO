#!/bin/bash

#SBATCH -A EXASGD
#SBATCH -p incline

#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -n 1
#SBATCH --gpu-bind=closest

#SBATCH -J exasgd_opflow_profiling
#SBATCH -o profile_opflow_incline_batch.%J
#SBATCH -e profile_opflow_incline_batch.%J

#SBATCH --mail-user=sayefsakin@gmail.com
#SBATCH --mail-type=ALL

source buildsystem/clang-hip/incline/base.sh
source buildsystem/spack/incline/modules/dependencies.sh
source buildsystem/spack/incline/modules/exago.sh

python3 performance_analysis/perf_pipeline.py performance_analysis/opflow_testsuite_incline.toml
