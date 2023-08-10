#!/bin/bash

#SBATCH -A csc359
#SBATCH -p batch
#SBATCH -C nvme

#SBATCH -t 2:00:00
#SBATCH -N 8000
#SBATCH --gpu-bind=closest
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=1

#SBATCH -J exasgd_sopflow_profiling_8000
#SBATCH -o profile_sopflow_frontier_batch_8000.%J
#SBATCH -e profile_sopflow_frontier_batch_8000.%J

#SBATCH --mail-user=sayefsakin@gmail.com
#SBATCH --mail-type=ALL

source buildsystem/clang-hip/frontier/base.sh
source buildsystem/spack/frontier/modules/dependencies.sh
#source buildsystem/spack/frontier/modules/exago.sh

export MPICH_ENV_DISPLAY=1
export MPICH_VERSION_DISPLAY=1
export MPICH_OFI_VERBOSE=1

export FI_CXI_RX_MATCH_MODE="software"
export FI_CXI_DEFAULT_CQ_SIZE=13107200
export FI_MR_CACHE_MONITOR=memhooks
export FI_OFI_RXM_RX_SIZE=13107200

export PMI_MMAP_SYNC_WAIT_TIME=1800

set -x

# sbcast
exe=sopflow
exe_fullpath=install/bin/$exe
echo "*****ldd ./${exe_fullpath}*****"
ldd ./${exe_fullpath}
echo "*************************"
sbcast --send-libs --exclude=NONE -pf ${exe_fullpath} /mnt/bb/$USER/${exe}
if [ ! "$?" == "0" ]; then
	echo "SBCAST failed!"
	exit 1
fi
echo "*****ls -lh /mnt/bb/$USER*****"
ls -lh /mnt/bb/$USER/
echo "*****ls -lh /mnt/bb/$USER/${exe}_libs*****"
ls -lh /mnt/bb/$USER/${exe}_libs
export LD_LIBRARY_PATH="/mnt/bb/$USER/${exe}_libs"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:$(pkg-config --variable=libdir libfabric)"
srun -N ${SLURM_NNODES} -n ${SLURM_NNODES} --ntasks-per-node=1 --label -D /mnt/bb/$USER/${exe}_libs \
	bash -c "if [ -f libhsa-runtime64.so.1 ]; then ln -s libhsa-runtime64.so.1 libhsa-runtime64.so; fi;
	if [ -f libamdhip64.so.5 ]; then ln -s libamdhip64.so.5 libamdhip64.so; fi"
export ROCBLAS_TENSILE_LIBPATH=${ROCM_PATH}/lib/rocblas/library
echo "*****ldd /mnt/bb/$USER/${exe}*****"
ldd /mnt/bb/$USER/${exe}
echo "*************************************"


python3 performance_analysis/perf_pipeline.py performance_analysis/sopflow_testsuite_frontier_8000.toml

set +x
