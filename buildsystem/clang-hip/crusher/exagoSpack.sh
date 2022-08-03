export MY_CLUSTER=crusher
export PROJ_DIR=/autofs/nccs-svm1_proj/csc359

module purge

module use -a /autofs/nccs-svm1_proj/csc359/cameron/spack/share/spack/modules/cray-sles15-zen3/

# Spack modules
# exago@develop%clang@14.0.0~cuda+hiop~ipo+ipopt+mpi~python+raja+rocm amdgpu_target=gfx90a build_type=RelWithDebInfo dev_path=/ccs/home/rcruther/exago-git arch=cray-sles15-zen3
module load exago-develop-clang-14.0.0-xcsv3re

# System modules

module load rocm/5.2.0
module load libfabric/1.15.0.0

export CC=/opt/rocm-5.2.0/llvm/bin/clang
export CXX=/opt/rocm-5.2.0/llvm/bin/clang++
export FC=/opt/rocm-5.2.0/llvm/bin/flang

