source /etc/profile.d/modules.sh
export MY_CLUSTER=marianas

module purge
module use -a /qfs/projects/exasgd/src/jaelyn-spack/spack/share/spack/modules/linux-centos7-haswell/
module use -a /qfs/projects/exasgd/src/jaelyn-spack/spack/share/spack/modules/linux-centos7-broadwell/
#module use -a /qfs/projects/exasgd/src/jaelyn-spack/spack/opt/spack/linux-centos7-broadwell/
#module use -a /qfs/projects/exasgd/src/jaelyn-spack/spack/opt/spack/linux-centos7-haswell/
# Load spack modules

# blt@0.4.1%gcc@7.3.0 arch=linux-centos7-broadwell
module load blt-0.4.1-gcc-7.3.0-twuxein
# camp@0.2.2%gcc@7.3.0+cuda~ipo~rocm~tests build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-broadwell
module load camp-0.2.2-gcc-7.3.0-u3kctft
# cmake@3.21.4%gcc@7.3.0~doc+ncurses+openssl+ownlibs~qt build_type=Release arch=linux-centos7-broadwell
module load cmake-3.21.4-gcc-7.3.0-amjh4ql
# coinhsl@2015.06.23%gcc@7.3.0+blas arch=linux-centos7-broadwell
module load coinhsl-2015.06.23-gcc-7.3.0-x5z66iv
# cub@1.12.0-rc0%gcc@7.3.0 arch=linux-centos7-broadwell
module load cub-1.12.0-rc0-gcc-7.3.0-4zav6ns
# hdf5@1.12.1%gcc@7.3.0~cxx~fortran~hl~ipo~java+mpi+shared~szip~threadsafe+tools api=default build_type=RelWithDebInfo patches=ee351eb arch=linux-centos7-broadwell
module load hdf5-1.12.1-gcc-7.3.0-zlqhesj
# hiop@0.5.4%gcc@7.3.0+cuda+deepchecking~ipo~jsrun~kron+mpi+raja~rocm~shared+sparse build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-broadwell
module load hiop-0.5.4-gcc-7.3.0-huj7w7e
# hypre@2.24.0%gcc@7.3.0~complex~cuda~debug+fortran~gptune~int64~internal-superlu~mixedint+mpi~openmp~rocm+shared~superlu-dist~unified-memory arch=linux-centos7-broadwell
module load hypre-2.24.0-gcc-7.3.0-t54hlvi
# ipopt@3.12.10%gcc@7.3.0+coinhsl+debug~metis~mumps arch=linux-centos7-broadwell
module load ipopt-3.12.10-gcc-7.3.0-wfu2hcv
# magma@2.6.2rc1%gcc@7.3.0+cuda+fortran~ipo~rocm+shared build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-broadwell
module load magma-2.6.2rc1-gcc-7.3.0-aoq7o3h
# metis@5.1.0%gcc@7.3.0~gdb~int64~real64+shared build_type=Release patches=4991da9,b1225da arch=linux-centos7-broadwell
module load metis-5.1.0-gcc-7.3.0-v22vb3o
# openblas@0.3.10%gcc@7.3.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared patches=865703b symbol_suffix=none threads=none arch=linux-centos7-broadwell
module load openblas-0.3.10-gcc-7.3.0-me5hbx4
# parmetis@4.0.3%gcc@7.3.0~gdb~int64~ipo+shared build_type=RelWithDebInfo patches=4f89253,50ed208,704b84f arch=linux-centos7-broadwell
module load parmetis-4.0.3-gcc-7.3.0-4axmtgj
# petsc@3.16.0%gcc@7.3.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw~giflib+hdf5~hpddm~hwloc+hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind+metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse+superlu-dist~tetgen~trilinos~valgrind clanguage=C arch=linux-centos7-broadwell
module load petsc-3.16.0-gcc-7.3.0-urp7kkz
# pkgconf@1.8.0%gcc@7.3.0 arch=linux-centos7-broadwell
module load pkgconf-1.8.0-gcc-7.3.0-eoghvsn
# py-attrs@21.4.0%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-attrs-21.4.0-gcc-7.3.0-fingolm
# py-iniconfig@1.1.1%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-iniconfig-1.1.1-gcc-7.3.0-6ici4ro
# py-mpi4py@3.1.2%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-mpi4py-3.1.2-gcc-7.3.0-evfao4t
# py-packaging@21.3%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-packaging-21.3-gcc-7.3.0-ipe3rwk
# py-pluggy@1.0.0%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-pluggy-1.0.0-gcc-7.3.0-n7x6spr
# py-py@1.11.0%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-py-1.11.0-gcc-7.3.0-yzbe2hh
# py-pyparsing@3.0.6%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-pyparsing-3.0.6-gcc-7.3.0-wxtnhh3
# py-pytest@6.2.5%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-pytest-6.2.5-gcc-7.3.0-ysg63c4
# py-setuptools@59.4.0%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-setuptools-59.4.0-gcc-7.3.0-ysqfh5r
# py-toml@0.10.2%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-toml-0.10.2-gcc-7.3.0-e2z6qnb
# python@3.8.5%gcc@7.3.0+bz2+ctypes+dbm~debug+ensurepip+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tix~tkinter~ucs4+uuid+zlib patches=0d98e93,4c24573,f2fd060 arch=linux-centos7-broadwell
module load python-3.8.5-gcc-7.3.0-ahml2fl
# raja@0.14.0%gcc@7.3.0+cuda~examples~exercises~ipo+openmp~rocm~shared~tests build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-broadwell
module load raja-0.14.0-gcc-7.3.0-ui2yfv6
# superlu-dist@7.2.0%gcc@7.3.0~cuda~int64~ipo~openmp~rocm+shared build_type=RelWithDebInfo patches=8da9e21 arch=linux-centos7-broadwell
module load superlu-dist-7.2.0-gcc-7.3.0-3xzqgmz
# umpire@6.0.0%gcc@7.3.0+c+cuda~deviceconst~examples~fortran~ipo~numa~openmp~rocm~shared build_type=RelWithDebInfo cuda_arch=60 tests=none arch=linux-centos7-broadwell
module load umpire-6.0.0-gcc-7.3.0-pq7d6om
# zlib@1.2.11%gcc@7.3.0+optimize+pic+shared arch=linux-centos7-broadwell
module load zlib-1.2.11-gcc-7.3.0-43cnepi

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
