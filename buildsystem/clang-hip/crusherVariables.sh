export MY_CLUSTER=crusher
export PROJ_DIR=/autofs/nccs-svm1_proj/csc359

module purge

module use -a /autofs/nccs-svm1_proj/csc359/cameron/spack/share/spack/modules/cray-sles15-zen3/

# Spack modules

# blt@0.4.1%clang@14.0.0 arch=cray-sles15-zen3
module load blt-0.4.1-clang-14.0.0-cgzud7l
# camp@0.2.3%clang@14.0.0~cuda~ipo+rocm~tests amdgpu_target=gfx90a build_type=RelWithDebInfo arch=cray-sles15-zen3
module load camp-0.2.3-clang-14.0.0-etswkoo
# cmake@3.20.6%clang@14.0.0~doc+ncurses+ownlibs~qt build_type=Release arch=cray-sles15-zen3
module load cmake-3.20.6-clang-14.0.0-kbznxux
# coinhsl@2019.05.21%clang@14.0.0+blas arch=cray-sles15-zen3
module load coinhsl-2019.05.21-clang-14.0.0-jmikoux
# cray-mpich@8.1.16%gcc@11.2.0+wrappers arch=cray-sles15-zen3
module load cray-mpich-8.1.16-gcc-11.2.0-yfbzrpg
# exago@develop%clang@14.0.0~cuda+hiop~ipo+ipopt+mpi~python+raja+rocm amdgpu_target=gfx90a build_type=RelWithDebInfo dev_path=/ccs/home/rcruther/exago-git arch=cray-sles15-zen3
module load exago-develop-clang-14.0.0-xcsv3re
# hiop@develop%clang@14.0.0~cuda~cusolver~deepchecking~ginkgo~ipo~jsrun+kron+mpi+raja+rocm~shared+sparse amdgpu_target=gfx90a build_type=RelWithDebInfo dev_path=/ccs/home/rcruther/hiop-git arch=cray-sles15-zen3
module load hiop-develop-clang-14.0.0-3n75cpk
# hip@5.2.0%clang@14.0.0~ipo build_type=Release patches=3b25db8 arch=cray-sles15-zen3
module load hip-5.2.0-clang-14.0.0-nwiinlu
# hipblas@5.2.0%clang@14.0.0~ipo build_type=Release arch=cray-sles15-zen3
module load hipblas-5.2.0-clang-14.0.0-o3wqq2c
# hipsparse@5.2.0%clang@14.0.0~ipo build_type=Release arch=cray-sles15-zen3
module load hipsparse-5.2.0-clang-14.0.0-d5iyoed
# hsa-rocr-dev@5.2.0%clang@14.0.0+image~ipo+shared build_type=Release patches=71e6851 arch=cray-sles15-zen3
module load hsa-rocr-dev-5.2.0-clang-14.0.0-q7lryu5
# ipopt@3.12.10%clang@14.0.0+coinhsl~debug~metis~mumps arch=cray-sles15-zen3
module load ipopt-3.12.10-clang-14.0.0-gebu4zk
# libiconv@1.16%gcc@11.2.0 libs=shared,static arch=cray-sles15-zen3
module load libiconv-1.16-gcc-11.2.0-po4cjxv
# llvm-amdgpu@5.2.0%clang@14.0.0~ipo~link_llvm_dylib~llvm_dylib~openmp+rocm-device-libs build_type=Release arch=cray-sles15-zen3
module load llvm-amdgpu-5.2.0-clang-14.0.0-khrhpxd
# magma@2.6.2%clang@14.0.0~cuda+fortran~ipo+rocm+shared amdgpu_target=gfx90a build_type=RelWithDebInfo arch=cray-sles15-zen3
module load magma-2.6.2-clang-14.0.0-ysldosl
# metis@5.1.0%clang@14.0.0~gdb~int64~real64+shared build_type=Release patches=4991da9 arch=cray-sles15-zen3
module load metis-5.1.0-clang-14.0.0-lky4kp3
# openblas@0.3.20%clang@14.0.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared symbol_suffix=none threads=openmp arch=cray-sles15-zen3
module load openblas-0.3.20-clang-14.0.0-dpeotz6
# perl@5.34.0%clang@14.0.0+cpanm+shared+threads arch=cray-sles15-zen3
module load perl-5.34.0-clang-14.0.0-3t3ib3e
# petsc@3.16.4-cpu%clang@14.0.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib+hdf5~hpddm~hwloc+hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind+metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse+superlu-dist~tetgen~trilinos~valgrind clanguage=C arch=cray-sles15-zen3
module load petsc-3.16.4-cpu-clang-14.0.0-gwgkmep
# pkgconf@1.8.0%clang@14.0.0 arch=cray-sles15-zen3
module load pkgconf-1.8.0-clang-14.0.0-67lrujf
# raja@0.14.0%clang@14.0.0~cuda~examples~exercises~ipo~openmp+rocm+shared~tests amdgpu_target=gfx90a build_type=RelWithDebInfo arch=cray-sles15-zen3
module load raja-0.14.0-clang-14.0.0-zbetmwb
# suite-sparse@4.5.6%clang@14.0.0~cuda~graphblas~openmp+pic~tbb arch=cray-sles15-zen3
module load suite-sparse-4.5.6-clang-14.0.0-vozkztd
# umpire@6.0.0%clang@14.0.0+c~cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp+rocm+shared amdgpu_target=gfx90a build_type=RelWithDebInfo tests=none arch=cray-sles15-zen3
module load umpire-6.0.0-clang-14.0.0-6kexnzi

# System modules

module load rocm/5.2.0
module load libfabric/1.15.0.0

export CC=/opt/rocm-5.2.0/llvm/bin/clang
export CXX=/opt/rocm-5.2.0/llvm/bin/clang++
export FC=/opt/rocm-5.2.0/llvm/bin/flang

[ -f $PWD/nvblas.conf ] && rm $PWD/nvblas.conf
cat > $PWD/nvblas.conf <<-EOD
NVBLAS_LOGFILE  nvblas.log
NVBLAS_CPU_BLAS_LIB $OPENBLAS_LIBRARY_DIR/libopenblas.so
NVBLAS_GPU_LIST ALL
NVBLAS_TILE_DIM 2048
NVBLAS_AUTOPIN_MEM_ENABLED
EOD
export NVBLAS_CONFIG_FILE=$PWD/nvblas.conf
echo "Generated $PWD/nvblas.conf"
