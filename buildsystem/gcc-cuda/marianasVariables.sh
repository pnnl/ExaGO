source /etc/profile.d/modules.sh
export MY_CLUSTER=marianas

module purge
module use -a /qfs/projects/exasgd/src/jaelyn-spack/spack/share/spack/modules/linux-centos7-zen2

# Load spack modules
# blt@0.4.1%gcc@10.2.0 arch=linux-centos7-zen2
module load blt-0.4.1-gcc-10.2.0-dfv7ktr
# camp@0.2.2%gcc@10.2.0+cuda~ipo~rocm~tests build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-zen2
module load camp-0.2.2-gcc-10.2.0-ooftid5
# cmake@3.21.4%gcc@10.2.0~doc+ncurses+ownlibs~qt build_type=Release arch=linux-centos7-zen2
module load cmake-3.21.4-gcc-10.2.0-mcdy53b
# coinhsl@2015.06.23%gcc@10.2.0+blas arch=linux-centos7-zen2
module load coinhsl-2015.06.23-gcc-10.2.0-rydelff
# cub@1.12.0%gcc@10.2.0 arch=linux-centos7-zen2
module load cub-1.12.0-gcc-10.2.0-5nooyca
# diffutils@3.8%gcc@10.2.0 arch=linux-centos7-zen2
module load diffutils-3.8-gcc-10.2.0-mjfwces
# hdf5@1.12.2%gcc@10.2.0~cxx~fortran~hl~ipo~java+mpi+shared~szip~threadsafe+tools api=default build_type=RelWithDebInfo arch=linux-centos7-zen2
module load hdf5-1.12.2-gcc-10.2.0-6dk553u
# hiop@develop%gcc@10.2.0+cuda~cusolver+deepchecking~ginkgo~ipo~jsrun~kron+mpi+raja~rocm~shared+sparse build_type=RelWithDebInfo cuda_arch=60 dev_path=/qfs/projects/exasgd/src/jaelyn-spack/dev-marianas/hiop arch=linux-centos7-zen2
module load hiop-develop-gcc-10.2.0-bgamkkk
# hypre@2.25.0%gcc@10.2.0~complex~cuda~debug+fortran~gptune~int64~internal-superlu~mixedint+mpi~openmp~rocm+shared~superlu-dist~unified-memory arch=linux-centos7-zen2
module load hypre-2.25.0-gcc-10.2.0-aq3hrgv
# ipopt@3.12.10%gcc@10.2.0+coinhsl+debug~metis~mumps arch=linux-centos7-zen2
module load ipopt-3.12.10-gcc-10.2.0-vzgarqw
# libiconv@1.16%gcc@10.2.0 libs=shared,static arch=linux-centos7-zen2
module load libiconv-1.16-gcc-10.2.0-gbg7l5p
# magma@2.6.2%gcc@10.2.0+cuda+fortran~ipo~rocm+shared build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-zen2
module load magma-2.6.2-gcc-10.2.0-bdcnq7a
# metis@5.1.0%gcc@10.2.0~gdb~int64~real64+shared build_type=Release patches=4991da9,b1225da arch=linux-centos7-zen2
module load metis-5.1.0-gcc-10.2.0-pysyrmt
# openblas@0.3.20%gcc@10.2.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared symbol_suffix=none threads=none arch=linux-centos7-zen2
module load openblas-0.3.20-gcc-10.2.0-qhcutll
# parmetis@4.0.3%gcc@10.2.0~gdb~int64~ipo+shared build_type=RelWithDebInfo patches=4f89253,50ed208,704b84f arch=linux-centos7-zen2
module load parmetis-4.0.3-gcc-10.2.0-5sl3pqo
# perl@5.26.0%gcc@10.2.0+cpanm+shared+threads patches=0eac10e,8cf4302 arch=linux-centos7-zen2
module load perl-5.26.0-gcc-10.2.0-l2yiybo
# petsc@3.16.6%gcc@10.2.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib+hdf5~hpddm~hwloc+hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind+metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse+superlu-dist~tetgen~trilinos~valgrind clanguage=C arch=linux-centos7-zen2
module load petsc-3.16.6-gcc-10.2.0-neo2xzl
# pkgconf@1.8.0%gcc@10.2.0 arch=linux-centos7-zen2
module load pkgconf-1.8.0-gcc-10.2.0-fuflwbl
# py-attrs@21.4.0%gcc@10.2.0 arch=linux-centos7-zen2
module load py-attrs-21.4.0-gcc-10.2.0-3lktnre
# py-flit-core@3.6.0%gcc@10.2.0 arch=linux-centos7-zen2
module load py-flit-core-3.6.0-gcc-10.2.0-5zrgn4g
# py-iniconfig@1.1.1%gcc@10.2.0 arch=linux-centos7-zen2
module load py-iniconfig-1.1.1-gcc-10.2.0-5js2b7d
# py-mpi4py@3.1.2%gcc@10.2.0 arch=linux-centos7-zen2
module load py-mpi4py-3.1.2-gcc-10.2.0-y5h43w5
# py-packaging@21.3%gcc@10.2.0 arch=linux-centos7-zen2
module load py-packaging-21.3-gcc-10.2.0-gumnw5v
# py-pip@21.3.1%gcc@10.2.0 arch=linux-centos7-zen2
module load py-pip-21.3.1-gcc-10.2.0-q6t5zmr
# py-pluggy@1.0.0%gcc@10.2.0 arch=linux-centos7-zen2
module load py-pluggy-1.0.0-gcc-10.2.0-r6wsldd
# py-py@1.11.0%gcc@10.2.0 arch=linux-centos7-zen2
module load py-py-1.11.0-gcc-10.2.0-jokegyi
# py-pyparsing@3.0.6%gcc@10.2.0 arch=linux-centos7-zen2
module load py-pyparsing-3.0.6-gcc-10.2.0-mq4ewhg
# py-pytest@6.2.5%gcc@10.2.0 arch=linux-centos7-zen2
module load py-pytest-6.2.5-gcc-10.2.0-ndntc24
# py-setuptools@63.0.0%gcc@10.2.0 arch=linux-centos7-zen2
module load py-setuptools-63.0.0-gcc-10.2.0-p3zj76s
# py-setuptools-scm@7.0.3%gcc@10.2.0+toml arch=linux-centos7-zen2
module load py-setuptools-scm-7.0.3-gcc-10.2.0-old2qf5
# py-toml@0.10.2%gcc@10.2.0 arch=linux-centos7-zen2
module load py-toml-0.10.2-gcc-10.2.0-nqxzcb5
# py-tomli@1.2.2%gcc@10.2.0 arch=linux-centos7-zen2
module load py-tomli-1.2.2-gcc-10.2.0-4y2zplr
# py-typing-extensions@4.3.0%gcc@10.2.0 arch=linux-centos7-zen2
module load py-typing-extensions-4.3.0-gcc-10.2.0-vyn77zo
# py-wheel@0.37.0%gcc@10.2.0 arch=linux-centos7-zen2
module load py-wheel-0.37.0-gcc-10.2.0-tqwco43
# raja@0.14.0%gcc@10.2.0+cuda~examples~exercises~ipo+openmp~rocm~shared~tests build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-zen2
module load raja-0.14.0-gcc-10.2.0-dwoe7eb
# superlu-dist@7.2.0%gcc@10.2.0~cuda~int64~ipo~openmp~rocm+shared build_type=RelWithDebInfo patches=8da9e21 arch=linux-centos7-zen2
module load superlu-dist-7.2.0-gcc-10.2.0-6hqy2pm
# umpire@6.0.0%gcc@10.2.0+c+cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp~rocm~shared build_type=RelWithDebInfo cuda_arch=60 tests=none arch=linux-centos7-zen2
module load umpire-6.0.0-gcc-10.2.0-ey5srfx
# zlib@1.2.12%gcc@10.2.0+optimize+pic+shared patches=0d38234 arch=linux-centos7-zen2
module load zlib-1.2.12-gcc-10.2.0-gnkqokp

# Load system modules
module load gcc/10.2.0
module load openmpi/4.1.0mlx5.0
module load cuda/11.4
module load python/miniconda3.8

source /share/apps/python/miniconda3.8/etc/profile.d/conda.sh

export CC=/share/apps/gcc/10.2.0/bin/gcc CXX=/share/apps/gcc/10.2.0/bin/g++ FC=/share/apps/gcc/10.2.0/bin/gfortran

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
