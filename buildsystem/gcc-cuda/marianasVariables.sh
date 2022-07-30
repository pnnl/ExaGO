source /etc/profile.d/modules.sh
export MY_CLUSTER=marianas

module purge
module use -a /qfs/projects/exasgd/src/jaelyn-spack/spack/share/spack/modules/linux-centos7-broadwell/
#module use -a /qfs/projects/exasgd/src/jaelyn-spack/spack/opt/spack/linux-centos7-broadwell/
#module use -a /qfs/projects/exasgd/src/jaelyn-spack/spack/opt/spack/linux-centos7-haswell/
# Load spack modules


# blt@0.4.1%gcc@7.3.0 arch=linux-centos7-broadwell
module load blt-0.4.1-gcc-7.3.0-urcrepw
# camp@0.2.2%gcc@7.3.0+cuda~ipo~rocm~tests build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-broadwell
module load camp-0.2.2-gcc-7.3.0-bsdzfpz
# cmake@3.21.4%gcc@7.3.0~doc+ncurses+ownlibs~qt build_type=Release arch=linux-centos7-broadwell
module load cmake-3.21.4-gcc-7.3.0-msd75ft
# coinhsl@2015.06.23%gcc@7.3.0+blas arch=linux-centos7-broadwell
module load coinhsl-2015.06.23-gcc-7.3.0-iju6fb2
# cub@1.12.0%gcc@7.3.0 arch=linux-centos7-broadwell
module load cub-1.12.0-gcc-7.3.0-d6bukcw
# diffutils@3.8%gcc@7.3.0 arch=linux-centos7-broadwell
module load diffutils-3.8-gcc-7.3.0-rosflkz
# hdf5@1.12.2%gcc@7.3.0~cxx~fortran~hl~ipo~java+mpi+shared~szip~threadsafe+tools api=default build_type=RelWithDebInfo arch=linux-centos7-broadwell
module load hdf5-1.12.2-gcc-7.3.0-oliewjd
# hiop@0.5.4%gcc@7.3.0+cuda~cusolver+deepchecking~ginkgo~ipo~jsrun~kron+mpi+raja~rocm~shared+sparse build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-broadwell
module load hiop-0.5.4-gcc-7.3.0-rpusc4l
# hypre@2.25.0%gcc@7.3.0~complex~cuda~debug+fortran~gptune~int64~internal-superlu~mixedint+mpi~openmp~rocm+shared~superlu-dist~unified-memory arch=linux-centos7-broadwell
module load hypre-2.25.0-gcc-7.3.0-ntiuixk
# ipopt@3.12.10%gcc@7.3.0+coinhsl+debug~metis~mumps arch=linux-centos7-broadwell
module load ipopt-3.12.10-gcc-7.3.0-ewk6mdo
# libiconv@1.16%gcc@7.3.0 libs=shared,static arch=linux-centos7-broadwell
module load libiconv-1.16-gcc-7.3.0-vdlzfgy
# magma@2.6.2%gcc@7.3.0+cuda+fortran~ipo~rocm+shared build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-broadwell
module load magma-2.6.2-gcc-7.3.0-66zopod
# metis@5.1.0%gcc@7.3.0~gdb~int64~real64+shared build_type=Release patches=4991da9,b1225da arch=linux-centos7-broadwell
module load metis-5.1.0-gcc-7.3.0-azhv4ox
# openblas@0.3.10%gcc@7.3.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared patches=865703b symbol_suffix=none threads=none arch=linux-centos7-broadwell
module load openblas-0.3.10-gcc-7.3.0-lqo2sao
# parmetis@4.0.3%gcc@7.3.0~gdb~int64~ipo+shared build_type=RelWithDebInfo patches=4f89253,50ed208,704b84f arch=linux-centos7-broadwell
module load parmetis-4.0.3-gcc-7.3.0-x5ao6hv
# perl@5.26.0%gcc@7.3.0+cpanm+shared+threads patches=0eac10e,8cf4302 arch=linux-centos7-broadwell
module load perl-5.26.0-gcc-7.3.0-mpxxfdv
# petsc@3.16.6%gcc@7.3.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib+hdf5~hpddm~hwloc+hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind+metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse+superlu-dist~tetgen~trilinos~valgrind clanguage=C arch=linux-centos7-broadwell
module load petsc-3.16.6-gcc-7.3.0-l5v2r4b
# pkgconf@1.8.0%gcc@7.3.0 arch=linux-centos7-broadwell
module load pkgconf-1.8.0-gcc-7.3.0-q2p7apv
# py-attrs@21.4.0%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-attrs-21.4.0-gcc-7.3.0-sidoqqb
# py-flit-core@3.6.0%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-flit-core-3.6.0-gcc-7.3.0-cy3kawd
# py-iniconfig@1.1.1%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-iniconfig-1.1.1-gcc-7.3.0-gnw4kxw
# py-mpi4py@3.1.2%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-mpi4py-3.1.2-gcc-7.3.0-njrmiuk
# py-packaging@21.3%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-packaging-21.3-gcc-7.3.0-5ieg3ri
# py-pip@21.3.1%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-pip-21.3.1-gcc-7.3.0-lp3w2wp
# py-pluggy@1.0.0%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-pluggy-1.0.0-gcc-7.3.0-xqtimvk
# py-py@1.11.0%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-py-1.11.0-gcc-7.3.0-r7ovhef
# py-pyparsing@3.0.6%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-pyparsing-3.0.6-gcc-7.3.0-2focqy6
# py-pytest@6.2.5%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-pytest-6.2.5-gcc-7.3.0-3pbbxrc
# py-setuptools@63.0.0%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-setuptools-63.0.0-gcc-7.3.0-3vai223
# py-setuptools-scm@7.0.3%gcc@7.3.0+toml arch=linux-centos7-broadwell
module load py-setuptools-scm-7.0.3-gcc-7.3.0-57n5gyr
# py-toml@0.10.2%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-toml-0.10.2-gcc-7.3.0-buq67gm
# py-tomli@1.2.2%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-tomli-1.2.2-gcc-7.3.0-ge2rqkc
# py-typing-extensions@4.3.0%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-typing-extensions-4.3.0-gcc-7.3.0-whnietr
# py-wheel@0.37.0%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-wheel-0.37.0-gcc-7.3.0-l7r3n6d
# raja@0.14.0%gcc@7.3.0+cuda+disable_blt_export~examples~exercises~ipo+openmp~rocm~shared~tests build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-broadwell
module load raja-0.14.0-gcc-7.3.0-nzn2nni
# superlu-dist@7.2.0%gcc@7.3.0~cuda~int64~ipo~openmp~rocm+shared build_type=RelWithDebInfo patches=8da9e21 arch=linux-centos7-broadwell
module load superlu-dist-7.2.0-gcc-7.3.0-tvimfqa
# umpire@6.0.0%gcc@7.3.0+c+cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp~rocm~shared build_type=RelWithDebInfo cuda_arch=60 tests=none arch=linux-centos7-broadwell
module load umpire-6.0.0-gcc-7.3.0-vppxkyu
# zlib@1.2.12%gcc@7.3.0+optimize+pic+shared patches=0d38234 arch=linux-centos7-broadwell
module load zlib-1.2.12-gcc-7.3.0-ywxvkcu


# Load system modules
module load gcc/7.3.0
module load openmpi/3.1.3
module load cuda/10.2.89

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
