export MY_CLUSTER=summit
export PROJ_DIR=/autofs/nccs-svm1_proj/csc359

module purge

module use -a /autofs/nccs-svm1_proj/csc359/shared_installs/spack/share/spack/modules/linux-rhel8-ppc64le

# blt@0.4.1%gcc@10.2.0 arch=linux-rhel8-ppc64le
module load exasgd-blt/0.4.1/gcc-10.2.0-bz5htrf
# camp@0.2.2%gcc@10.2.0+cuda~ipo~rocm~tests build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel8-ppc64le
module load exasgd-camp/0.2.2/cuda-11.5.2/gcc-10.2.0-dc42bvq
# cmake@3.21.3%gcc@10.2.0~doc+ncurses+ownlibs~qt build_type=Release arch=linux-rhel8-ppc64le
module load exasgd-cmake/3.21.3/gcc-10.2.0-jtqszri
# coinhsl@2015.06.23%gcc@10.2.0+blas arch=linux-rhel8-ppc64le
module load exasgd-coinhsl/2015.06.23/gcc-10.2.0-btplbob
# cuda@11.5.2%gcc@10.2.0~allow-unsupported-compilers~dev arch=linux-rhel8-ppc64le
module load exasgd-cuda/11.5.2/gcc-10.2.0-7gofref
# ginkgo@glu_experimental%gcc@10.2.0+cuda~develtools~full_optimizations~hwloc~ipo~oneapi+openmp~rocm+shared build_type=Release cuda_arch=70 arch=linux-rhel8-ppc64le
module load exasgd-ginkgo/glu_experimental/cuda-11.5.2/gcc-10.2.0-4ybxcrk
# hiop@develop%gcc@10.2.0+cuda+cusolver+deepchecking+ginkgo~ipo~jsrun~kron+mpi+raja~rocm~shared+sparse build_type=Release cuda_arch=70 arch=linux-rhel8-ppc64le
module load exasgd-hiop/develop/cuda-11.5.2/spectrum-mpi-10.4.0.3-20210112/gcc-10.2.0-nsn5koh
# ipopt@3.12.10%gcc@10.2.0+coinhsl~debug+metis~mumps patches=712f729 arch=linux-rhel8-ppc64le
module load exasgd-ipopt/3.12.10/gcc-10.2.0-5x34omc
# magma@2.6.2%gcc@10.2.0+cuda+fortran~ipo~rocm+shared build_type=Release cuda_arch=70 arch=linux-rhel8-ppc64le
module load exasgd-magma/2.6.2/cuda-11.5.2/gcc-10.2.0-pkzoqx7
# metis@5.1.0%gcc@10.2.0~gdb~int64~real64+shared build_type=Release patches=4991da9,b1225da arch=linux-rhel8-ppc64le
module load exasgd-metis/5.1.0/gcc-10.2.0-c5niajq
# openblas@0.17.0%gcc@10.2.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared symbol_suffix=none threads=none arch=linux-rhel8-ppc64le
module load exasgd-openblas/0.17.0/gcc-10.2.0-jd6kkan
# parmetis@4.0.3%gcc@10.2.0~gdb~int64~ipo+shared build_type=RelWithDebInfo patches=4f89253,50ed208,704b84f arch=linux-rhel8-ppc64le
module load exasgd-parmetis/4.0.3/spectrum-mpi-10.4.0.3-20210112/gcc-10.2.0-skqvpvq
# petsc@3.16.6%gcc@10.2.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind+metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~tetgen~trilinos~valgrind clanguage=C arch=linux-rhel8-ppc64le
module load exasgd-petsc/3.16.6/spectrum-mpi-10.4.0.3-20210112/gcc-10.2.0-cmc5awz
# py-attrs@21.4.0%gcc@10.2.0 arch=linux-rhel8-ppc64le
module load exasgd-py-attrs/21.4.0/gcc-10.2.0-5nsmg5d
# py-iniconfig@1.1.1%gcc@10.2.0 arch=linux-rhel8-ppc64le
module load exasgd-py-iniconfig/1.1.1/gcc-10.2.0-mrcocjs
# py-mpi4py@3.1.2%gcc@10.2.0 arch=linux-rhel8-ppc64le
module load exasgd-py-mpi4py/3.1.2/spectrum-mpi-10.4.0.3-20210112/gcc-10.2.0-isswj5r
# py-packaging@21.3%gcc@10.2.0 arch=linux-rhel8-ppc64le
module load exasgd-py-packaging/21.3/gcc-10.2.0-4bymaub
# py-pip@21.3.1%gcc@10.2.0 arch=linux-rhel8-ppc64le
module load exasgd-py-pip/21.3.1/gcc-10.2.0-ah6hfmn
# py-pluggy@1.0.0%gcc@10.2.0 arch=linux-rhel8-ppc64le
module load exasgd-py-pluggy/1.0.0/gcc-10.2.0-ywm32is
# py-py@1.11.0%gcc@10.2.0 arch=linux-rhel8-ppc64le
module load exasgd-py-py/1.11.0/gcc-10.2.0-agcnthq
# py-pyparsing@3.0.6%gcc@10.2.0 arch=linux-rhel8-ppc64le
module load exasgd-py-pyparsing/3.0.6/gcc-10.2.0-bvmgcsd
# py-pytest@6.2.5%gcc@10.2.0 arch=linux-rhel8-ppc64le
module load exasgd-py-pytest/6.2.5/gcc-10.2.0-4ybauco
# py-setuptools@59.4.0%gcc@10.2.0 arch=linux-rhel8-ppc64le
module load exasgd-py-setuptools/59.4.0/gcc-10.2.0-7z2lprw
# py-setuptools-scm@6.3.2%gcc@10.2.0+toml arch=linux-rhel8-ppc64le
module load exasgd-py-setuptools-scm/6.3.2/gcc-10.2.0-4dt7tcp
# py-toml@0.10.2%gcc@10.2.0 arch=linux-rhel8-ppc64le
module load exasgd-py-toml/0.10.2/gcc-10.2.0-n4na4o5
# py-tomli@1.2.2%gcc@10.2.0 arch=linux-rhel8-ppc64le
module load exasgd-py-tomli/1.2.2/gcc-10.2.0-aumjyxb
# py-wheel@0.37.0%gcc@10.2.0 arch=linux-rhel8-ppc64le
module load exasgd-py-wheel/0.37.0/gcc-10.2.0-qen5h6h
# python@3.8.10%gcc@10.2.0+bz2+ctypes+dbm~debug+ensurepip+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tix~tkinter~ucs4+uuid+zlib patches=0d98e93,4c24573,f2fd060 arch=linux-rhel8-ppc64le
module load exasgd-python/3.8.10/gcc-10.2.0-f72vdqp
# raja@0.14.0%gcc@10.2.0+cuda+examples+exercises~ipo+openmp~rocm+shared~tests build_type=Release cuda_arch=70 arch=linux-rhel8-ppc64le
module load exasgd-raja/0.14.0/cuda-11.5.2/gcc-10.2.0-ibkae64
# umpire@6.0.0%gcc@10.2.0~c+cuda~device_alloc~deviceconst+examples~fortran~ipo~numa~openmp~rocm~shared build_type=Release cuda_arch=70 tests=none arch=linux-rhel8-ppc64le
module load exasgd-umpire/6.0.0/cuda-11.5.2/gcc-10.2.0-wlm3tr7

# Load system modules
module load gcc/10.2.0
module load spectrum-mpi/10.4.0.3-20210112

export CC=/sw/summit/gcc/10.2.0-2/bin/gcc
export CXX=/sw/summit/gcc/10.2.0-2/bin/g++
export FC=/sw/summit/gcc/10.2.0-2/bin/gfotran

export MY_PETSC_DIR=$PETSC_DIR
export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DCMAKE_CUDA_ARCHITECTURES=70" # -DEXAGO_ENABLE_IPOPT=ON"
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
