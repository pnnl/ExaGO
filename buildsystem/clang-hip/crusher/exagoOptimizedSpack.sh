export MY_CLUSTER=crusher
export PROJ_DIR=/autofs/nccs-svm1_proj/csc359

module purge

module use -a /autofs/nccs-svm1_proj/csc359/cameron/spack/share/spack/test-modules/cray-sles15-zen3/

# Spack modules
export SRCDIR=${SRCDIR:-$PWD}
source $SRCDIR/buildsystem/clang-hip/crusher/exago-optimized.sh

# System modules

module load rocm/5.2.0
module load libfabric/1.15.0.0

export CC=/opt/rocm-5.2.0/llvm/bin/clang
export CXX=/opt/rocm-5.2.0/llvm/bin/clang++
export FC=/opt/rocm-5.2.0/llvm/bin/flang

