source /etc/profile.d/modules.sh
export MY_CLUSTER=marianas

module purge
module use -a /qfs/projects/exasgd/src/spack/share/spack/modules/linux-centos7-broadwell/
module use -a /qfs/projects/exasgd/src/spack/share/spack/modules/linux-centos7-x86_64/

# Load spack modules
# autoconf-archive@2019.01.06%gcc@7.3.0 arch=linux-centos7-x86_64
module load exasgd-autoconf-archive/2019.01.06/gcc-7.3.0-5a4lwsi
# blt@0.3.6%gcc@7.3.0 arch=linux-centos7-x86_64
module load exasgd-blt/0.3.6/gcc-7.3.0-frvx7ca
# camp@0.1.0%gcc@7.3.0+cuda~ipo~rocm~tests amdgpu_target=none build_type=RelWithDebInfo cuda_arch=none arch=linux-centos7-x86_64
module load exasgd-camp/0.1.0/cuda-10.2.89/gcc-7.3.0-a53tr7i
# cmake@3.20.3%gcc@7.3.0~doc+ncurses~openssl+ownlibs~qt build_type=Release arch=linux-centos7-x86_64
module load exasgd-cmake/3.20.3/gcc-7.3.0-qsnaili
# coinhsl@2015.06.23%gcc@7.3.0+blas arch=linux-centos7-x86_64
module load exasgd-coinhsl/2015.06.23/gcc-7.3.0-zv2yxcl
# cub@1.12.0-rc0%gcc@7.3.0 arch=linux-centos7-x86_64
module load exasgd-cub/1.12.0-rc0/gcc-7.3.0-cscwjgd
# gmp@6.2.1%gcc@7.3.0 arch=linux-centos7-x86_64
module load exasgd-gmp/6.2.1/gcc-7.3.0-sqgj2tv
# hdf5@1.10.7%gcc@7.3.0~cxx~debug~fortran~hl~java+mpi+pic+shared~szip~threadsafe api=none arch=linux-centos7-x86_64
module load exasgd-hdf5/1.10.7/openmpi-3.1.3/gcc-7.3.0-rpcf5xq
# hiop@0.4.6%gcc@8.1.0+cuda~deepchecking~ipo~jsrun+kron+mpi+raja+shared+sparse build_type=Release cuda_arch=80 arch=linux-centos7-zen
module load exasgd-hiop/0.4.6/cuda-10.2.89/openmpi-3.1.3/gcc-7.3.0-jeuhgas
# hypre@2.20.0%gcc@7.3.0~complex~cuda~debug~int64~internal-superlu~mixedint+mpi~openmp+shared~superlu-dist~unified-memory cuda_arch=none patches=6e3336b1d62155f6350dfe42b0f9ea25d4fa0af60c7e540959139deb93a26059 arch=linux-centos7-x86_64
module load exasgd-hypre/2.20.0/openmpi-3.1.3/gcc-7.3.0-kcfd7f4
# ipopt@3.12.10%gcc@7.3.0+coinhsl~debug~metis~mumps arch=linux-centos7-x86_64
module load exasgd-ipopt/3.12.10/gcc-7.3.0-pyjg6pd
# magma@2.6.1%gcc@8.1.0+cuda+fortran~ipo+shared build_type=RelWithDebInfo cuda_arch=80 arch=linux-centos7-zen
module load exasgd-magma/2.6.1/cuda-10.2.89/gcc-7.3.0-mb7xpns
# metis@5.1.0%gcc@7.3.0~gdb~int64~real64+shared build_type=Release patches=4991da938c1d3a1d3dea78e49bbebecba00273f98df2a656e38b83d55b281da1,b1225da886605ea558db7ac08dd8054742ea5afe5ed61ad4d0fe7a495b1270d2 arch=linux-centos7-x86_64
module load exasgd-metis/5.1.0/gcc-7.3.0-t2i65bz
# mpfr@4.0.2%gcc@7.3.0 patches=3f80b836948aa96f8d1cb9cc7f3f55973f19285482a96f9a4e1623d460bcccf0 arch=linux-centos7-x86_64
module load exasgd-mpfr/4.0.2/gcc-7.3.0-a72cnaj
# ncurses@6.2%gcc@7.3.0~symlinks+termlib abi=none arch=linux-centos7-x86_64
module load exasgd-ncurses/6.2/gcc-7.3.0-xus35rj
# openblas@0.3.15%gcc@4.8.5~bignuma~consistent_fpcsr~ilp64+locking+pic+shared threads=none arch=linux-centos7-x86_64
module load exasgd-openblas/0.3.15/gcc-4.8.5-hnfyfu5
# parmetis@4.0.3%gcc@7.3.0~gdb~int64~ipo+shared build_type=RelWithDebInfo patches=4f892531eb0a807eb1b82e683a416d3e35154a455274cf9b162fb02054d11a5b,50ed2081bc939269689789942067c58b3e522c269269a430d5d34c00edbc5870,704b84f7c7444d4372cb59cca6e1209df4ef3b033bc4ee3cf50f369bce972a9d arch=linux-centos7-x86_64
module load exasgd-parmetis/4.0.3/openmpi-3.1.3/gcc-7.3.0-y3hslfg
# petsc@3.14.6%gcc@7.3.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw~giflib+hdf5~hwloc+hypre~int64~jpeg~knl~libpng~libyaml~memkind+metis~mkl-pardiso~moab~mpfr+mpi~mumps~openmp~p4est~ptscotch~random123~rocm~saws+shared~suite-sparse+superlu-dist~trilinos~valgrind amdgpu_target=none clanguage=C cuda_arch=none arch=linux-centos7-x86_64
module load exasgd-petsc/3.14.6/openmpi-3.1.3/gcc-7.3.0-54tno2h
# pkgconf@1.7.4%gcc@7.3.0 arch=linux-centos7-x86_64
module load exasgd-pkgconf/1.7.4/gcc-7.3.0-6oic5tx
# raja@0.13.0%gcc@7.3.0+cuda~examples~exercises~ipo+openmp~rocm+shared~tests amdgpu_target=none build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-x86_64
module load exasgd-raja/0.13.0/cuda-10.2.89/gcc-7.3.0-dczmugk
# suite-sparse@5.10.1%gcc@7.3.0~cuda~openmp+pic~tbb arch=linux-centos7-x86_64
module load exasgd-suite-sparse/5.10.1/gcc-7.3.0-hyjqs5f
# superlu-dist@6.4.0%gcc@7.3.0~cuda~int64~ipo~openmp+shared build_type=RelWithDebInfo cuda_arch=none arch=linux-centos7-x86_64
module load exasgd-superlu-dist/6.4.0/openmpi-3.1.3/gcc-7.3.0-mt4a3ef
# umpire@4.1.2%gcc@7.3.0+c+cuda~deviceconst~examples~fortran~ipo~numa+openmp~rocm~shared amdgpu_target=none build_type=RelWithDebInfo cuda_arch=none patches=135bbc7d2f371531f432672b115ac0a407968aabfffc5b8a941db9b493dbf81f,7d912d31cd293df005ba74cb96c6f3e32dc3d84afff49b14509714283693db08 tests=none arch=linux-centos7-x86_64
module load exasgd-umpire/4.1.2/cuda-10.2.89/gcc-7.3.0-yjxsb2o
# zlib@1.2.11%gcc@7.3.0+optimize+pic+shared arch=linux-centos7-x86_64
module load exasgd-zlib/1.2.11/gcc-7.3.0-zknfwan

# Load system modules
module load gcc/7.3.0
module load openmpi/3.1.3
module load cuda/10.2.89
module load python/anaconda3.2019.3
source /share/apps/python/anaconda3.2019.3/etc/profile.d/conda.sh

export CC=/share/apps/gcc/7.3.0/bin/gcc CXX=/share/apps/gcc/7.3.0/bin/g++ FC=/share/apps/gcc/7.3.0/bin/gfortran

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

export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DCMAKE_CUDA_ARCHITECTURES=60"
