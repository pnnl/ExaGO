export MY_CLUSTER=ascent

module purge

module use -a /gpfs/wolf/proj-shared/csc359/src/spack/share/spack/modules/linux-rhel8-power9le

# Load spack modules
# blt@0.4.1%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-blt/0.4.1/gcc-10.2.0-ndtzigt
# camp@0.2.2%gcc@10.2.0+cuda~ipo~rocm~tests build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel8-power9le
module load exasgd-camp/0.2.2/cuda-11.4.2/gcc-10.2.0-2xxjut5
# cmake@3.22.2%gcc@10.2.0~doc+ncurses+ownlibs~qt build_type=Release arch=linux-rhel8-power9le
module load exasgd-cmake/3.22.2/gcc-10.2.0-sj55k3d
# coinhsl@2015.06.23%gcc@10.2.0+blas arch=linux-rhel8-power9le
module load exasgd-coinhsl/2015.06.23/gcc-10.2.0-zr37tzv
# cub@1.16.0%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-cub/1.16.0/gcc-10.2.0-zeukz4z
# cuda@11.4.2%gcc@10.2.0~allow-unsupported-compilers~dev arch=linux-rhel8-power9le
module load exasgd-cuda/11.4.2/gcc-10.2.0-yf564wn
# diffutils@3.8%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-diffutils/3.8/gcc-10.2.0-n2gttub
# ginkgo@1.5.0.glu_experimental%gcc@10.2.0+cuda~develtools~full_optimizations~hwloc~ipo~oneapi+openmp~rocm+shared build_system=cmake build_type=Debug cuda_arch=70 dev_path=/gpfs/wolf/proj-shared/csc359/src/ginkgo arch=linux-rhel8-power9le
module load exasgd-ginkgo/1.5.0.glu_experimental/cuda-11.4.2/gcc-10.2.0-h4v5snw
# gmp@6.2.1%gcc@10.2.0 libs=shared,static arch=linux-rhel8-power9le
module load exasgd-gmp/6.2.1/gcc-10.2.0-6oxkdwz
# gnuconfig@2021-08-14%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-gnuconfig/2021-08-14/gcc-10.2.0-j2pelq3
# hdf5@1.12.2%gcc@10.2.0~cxx~fortran~hl~ipo~java+mpi+shared~szip~threadsafe+tools api=default build_type=RelWithDebInfo arch=linux-rhel8-power9le
module load exasgd-hdf5/1.12.2/spectrum-mpi-10.4.0.3-20210112/gcc-10.2.0-yh642ea
# hiop@develop%gcc@10.2.0+cuda+cusolver+deepchecking+ginkgo~ipo~jsrun+kron+mpi+raja~rocm~shared+sparse build_system=cmake build_type=RelWithDebInfo cuda_arch=70 dev_path=/gpfs/wolf/proj-shared/csc359/src/nicholson-hiop arch=linux-rhel8-power9le
module load exasgd-hiop/develop/cuda-11.4.2/spectrum-mpi-10.4.0.3-20210112/gcc-10.2.0-ucuclxs
# hypre@2.24.0%gcc@10.2.0~complex~cuda~debug+fortran~gptune~int64~internal-superlu~mixedint+mpi~openmp~rocm+shared~superlu-dist~unified-memory arch=linux-rhel8-power9le
module load exasgd-hypre/2.24.0/spectrum-mpi-10.4.0.3-20210112/gcc-10.2.0-oybqs3s
# ipopt@3.12.10%gcc@10.2.0+coinhsl~debug~metis~mumps build_system=autotools arch=linux-rhel8-power9le
module load exasgd-ipopt/3.12.10/gcc-10.2.0-i4gqiqg
# libiconv@1.16%gcc@10.2.0 libs=shared,static arch=linux-rhel8-power9le
module load exasgd-libiconv/1.16/gcc-10.2.0-cpn6euq
# magma@2.6.2%gcc@10.2.0+cuda+fortran~ipo~rocm+shared build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel8-power9le
module load exasgd-magma/2.6.2/cuda-11.4.2/gcc-10.2.0-ffujy7o
# metis@5.1.0%gcc@10.2.0~gdb~int64~real64+shared build_type=Release patches=4991da9,b1225da arch=linux-rhel8-power9le
module load exasgd-metis/5.1.0/gcc-10.2.0-kdt6x4i
# mpfr@4.1.0%gcc@10.2.0 libs=shared,static arch=linux-rhel8-power9le
module load exasgd-mpfr/4.1.0/gcc-10.2.0-ua5zak3
# openblas@0.3.20%gcc@10.2.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared symbol_suffix=none threads=none arch=linux-rhel8-power9le
module load exasgd-openblas/0.3.20/gcc-10.2.0-lkjyfws
# parmetis@4.0.3%gcc@10.2.0~gdb~int64~ipo+shared build_type=RelWithDebInfo patches=4f89253,50ed208,704b84f arch=linux-rhel8-power9le
module load exasgd-parmetis/4.0.3/spectrum-mpi-10.4.0.3-20210112/gcc-10.2.0-gqz32yk
# petsc@3.18.1%gcc@10.2.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib+hdf5~hpddm~hwloc+hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind+metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse+superlu-dist~tetgen~trilinos~valgrind build_system=generic clanguage=C arch=linux-rhel8-power9le
module load exasgd-petsc/3.18.1/spectrum-mpi-10.4.0.3-20210112/gcc-10.2.0-cwok3k4
# pkgconf@1.8.0%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-pkgconf/1.8.0/gcc-10.2.0-nngfigd
# py-attrs@21.4.0%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-py-attrs/21.4.0/gcc-10.2.0-av2tzhn
# py-flit-core@3.6.0%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-py-flit-core/3.6.0/gcc-10.2.0-wdoqgro
# py-importlib-metadata@4.8.2%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-py-importlib-metadata/4.8.2/gcc-10.2.0-klg3xix
# py-iniconfig@1.1.1%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-py-iniconfig/1.1.1/gcc-10.2.0-o6bmfak
# py-mpi4py@3.1.2%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-py-mpi4py/3.1.2/spectrum-mpi-10.4.0.3-20210112/gcc-10.2.0-q53tcvd
# py-packaging@21.3%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-py-packaging/21.3/gcc-10.2.0-ej3e6zm
# py-pip@21.3.1%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-py-pip/21.3.1/gcc-10.2.0-awkoczv
# py-pluggy@1.0.0%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-py-pluggy/1.0.0/gcc-10.2.0-tlqdjwe
# py-py@1.11.0%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-py-py/1.11.0/gcc-10.2.0-fwbtxbd
# py-pyparsing@3.0.6%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-py-pyparsing/3.0.6/gcc-10.2.0-qxfb6l5
# py-pytest@6.2.5%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-py-pytest/6.2.5/gcc-10.2.0-2vv5hdj
# py-setuptools@59.4.0%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-py-setuptools/59.4.0/gcc-10.2.0-elx3acw
# py-setuptools-scm@6.3.2%gcc@10.2.0+toml arch=linux-rhel8-power9le
module load exasgd-py-setuptools-scm/6.3.2/gcc-10.2.0-xdvaqbm
# py-toml@0.10.2%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-py-toml/0.10.2/gcc-10.2.0-6lhs4ix
# py-tomli@1.2.2%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-py-tomli/1.2.2/gcc-10.2.0-4glf4eg
# py-typing-extensions@4.1.1%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-py-typing-extensions/4.1.1/gcc-10.2.0-qg6k6g2
# py-wheel@0.37.0%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-py-wheel/0.37.0/gcc-10.2.0-tg3twdl
# py-zipp@3.6.0%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-py-zipp/3.6.0/gcc-10.2.0-4vxciug
# python@3.6.8%gcc@10.2.0+bz2+ctypes+dbm~debug+ensurepip+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tix~tkinter~ucs4+uuid+zlib patches=c129e34 arch=linux-rhel8-power9le
module load exasgd-python/3.6.8/gcc-10.2.0-273robh
# raja@0.14.0%gcc@10.2.0+cuda~examples~exercises~ipo+openmp~rocm+shared~tests build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel8-power9le
module load exasgd-raja/0.14.0/cuda-11.4.2/gcc-10.2.0-hkkig4v
# spectrum-mpi@10.4.0.3-20210112%gcc@10.2.0 arch=linux-rhel8-power9le
module load exasgd-spectrum-mpi/10.4.0.3-20210112/gcc-10.2.0-c5uqeqr
# suite-sparse@5.10.1%gcc@10.2.0~cuda~graphblas~openmp+pic~tbb arch=linux-rhel8-power9le
module load exasgd-suite-sparse/5.10.1/gcc-10.2.0-sq5isrj
# superlu-dist@7.2.0%gcc@10.2.0~cuda~int64~ipo~openmp~rocm+shared build_type=RelWithDebInfo patches=8da9e21 arch=linux-rhel8-power9le
module load exasgd-superlu-dist/7.2.0/spectrum-mpi-10.4.0.3-20210112/gcc-10.2.0-7cwmhkg
# umpire@6.0.0%gcc@10.2.0~c+cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp~rocm~shared build_type=RelWithDebInfo cuda_arch=70 tests=none arch=linux-rhel8-power9le
module load exasgd-umpire/6.0.0/cuda-11.4.2/gcc-10.2.0-selkipw
# zlib@1.2.12%gcc@10.2.0+optimize+pic+shared patches=0d38234 arch=linux-rhel8-power9le
module load exasgd-zlib/1.2.12/gcc-10.2.0-lfmc24z

# Load system modules
module load cuda/11.4.2
module load gcc/10.2.0
module load spectrum-mpi/10.4.0.3-20210112
module load cmake/3.22.2

export MY_PETSC_DIR=$PETSC_DIR
export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DEXAGO_CTEST_LAUNCH_COMMAND='jsrun -g 1 -n 1'"

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

export CC=$(which gcc)
export CXX=$(which g++)
export FC=$(which gfortran)
