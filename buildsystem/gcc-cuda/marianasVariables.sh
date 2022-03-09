source /etc/profile.d/modules.sh
export MY_CLUSTER=marianas

module purge
module use -a /qfs/projects/exasgd/src/cameron-spack/share/spack/modules/linux-centos7-haswell/
module use -a /qfs/projects/exasgd/src/cameron-spack/share/spack/modules/linux-centos7-broadwell/

# Load spack modules
# autoconf@2.69%gcc@7.3.0 patches=35c449281546376449766f92d49fc121ca50e330e60fefcfc9be2af3253082c2,7793209b33013dc0f81208718c68440c5aae80e7a1c4b8d336e382525af791a7,a49dd5bac3b62daa0ff688ab4d508d71dbd2f4f8d7e2a02321926346161bf3ee arch=linux-centos7-broadwell
module load autoconf-2.69-gcc-7.3.0-nqowopj
# autoconf-archive@2019.01.06%gcc@7.3.0 arch=linux-centos7-broadwell
module load autoconf-archive-2019.01.06-gcc-7.3.0-genbzzs
# automake@1.16.5%gcc@7.3.0 arch=linux-centos7-broadwell
module load automake-1.16.5-gcc-7.3.0-e5337oh
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
# diffutils@3.8%gcc@7.3.0 arch=linux-centos7-broadwell
module load diffutils-3.8-gcc-7.3.0-2fxdbir
# gmp@6.2.1%gcc@7.3.0 arch=linux-centos7-broadwell
module load gmp-6.2.1-gcc-7.3.0-myfvhge
# hdf5@1.12.1%gcc@7.3.0~cxx~fortran~hl~ipo~java+mpi+shared~szip~threadsafe+tools api=default build_type=RelWithDebInfo arch=linux-centos7-broadwell
module load hdf5-1.12.1-gcc-7.3.0-57bsdkz
# hiop@0.5.3%gcc@7.3.0+cuda+deepchecking~ipo~jsrun+kron+mpi+raja~rocm~shared+sparse build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-broadwell
module load hiop-0.5.3-gcc-7.3.0-k4p3qky
# hypre@2.23.0%gcc@7.3.0~complex~cuda~debug+fortran~gptune~int64~internal-superlu~mixedint+mpi~openmp+shared~superlu-dist~unified-memory arch=linux-centos7-broadwell
module load hypre-2.23.0-gcc-7.3.0-zb2yfer
# ipopt@3.12.10%gcc@7.3.0+coinhsl+debug~metis~mumps arch=linux-centos7-broadwell
module load ipopt-3.12.10-gcc-7.3.0-wfu2hcv
# libiconv@1.16%gcc@7.3.0 libs=shared,static arch=linux-centos7-broadwell
module load libiconv-1.16-gcc-7.3.0-465thn2
# libsigsegv@2.13%gcc@7.3.0 arch=linux-centos7-broadwell
module load libsigsegv-2.13-gcc-7.3.0-7wzsrah
# libtool@2.4.6%gcc@7.3.0 arch=linux-centos7-broadwell
module load libtool-2.4.6-gcc-7.3.0-ftdeyga
# m4@1.4.19%gcc@7.3.0+sigsegv patches=9dc5fbd0d5cb1037ab1e6d0ecc74a30df218d0a94bdd5a02759a97f62daca573,bfdffa7c2eb01021d5849b36972c069693654ad826c1a20b53534009a4ec7a89 arch=linux-centos7-broadwell
module load m4-1.4.19-gcc-7.3.0-e6s2uj7
# magma@2.6.1%gcc@7.3.0+cuda+fortran~ipo~rocm+shared build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-broadwell
module load magma-2.6.1-gcc-7.3.0-wkihdjh
# metis@5.1.0%gcc@7.3.0~gdb~int64~real64+shared build_type=Release patches=4991da938c1d3a1d3dea78e49bbebecba00273f98df2a656e38b83d55b281da1,b1225da886605ea558db7ac08dd8054742ea5afe5ed61ad4d0fe7a495b1270d2 arch=linux-centos7-broadwell
module load metis-5.1.0-gcc-7.3.0-v22vb3o
# mpfr@4.1.0%gcc@7.3.0 arch=linux-centos7-broadwell
module load mpfr-4.1.0-gcc-7.3.0-7f747eu
# openblas@0.3.10%gcc@7.3.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared patches=865703b4f405543bbd583413fdeff2226dfda908be33639276c06e5aa7ae2cae symbol_suffix=none threads=none arch=linux-centos7-broadwell
module load openblas-0.3.10-gcc-7.3.0-me5hbx4
# openmpi@3.1.3%gcc@7.3.0~atomics~cuda~cxx~cxx_exceptions+gpfs~internal-hwloc~java~legacylaunchers~lustre~memchecker~pmi~pmix+romio~singularity~sqlite3+static~thread_multiple+vt+wrapper-rpath fabrics=none schedulers=none arch=linux-centos7-broadwell
module load openmpi-3.1.3-gcc-7.3.0-kbx2cvy
# parmetis@4.0.3%gcc@7.3.0~gdb~int64~ipo+shared build_type=RelWithDebInfo patches=4f892531eb0a807eb1b82e683a416d3e35154a455274cf9b162fb02054d11a5b,50ed2081bc939269689789942067c58b3e522c269269a430d5d34c00edbc5870,704b84f7c7444d4372cb59cca6e1209df4ef3b033bc4ee3cf50f369bce972a9d arch=linux-centos7-broadwell
module load parmetis-4.0.3-gcc-7.3.0-5t7icxb
# perl@5.26.0%gcc@7.3.0+cpanm+shared+threads patches=0eac10ed90aeb0459ad8851f88081d439a4e41978e586ec743069e8b059370ac,8cf4302ca8b480c60ccdcaa29ec53d9d50a71d4baf469ac8c6fca00ca31e58a2 arch=linux-centos7-broadwell
module load perl-5.26.0-gcc-7.3.0-glzxh6e
# petsc@3.16.0%gcc@7.3.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw~giflib+hdf5~hpddm~hwloc+hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind+metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse+superlu-dist~tetgen~trilinos~valgrind clanguage=C arch=linux-centos7-broadwell
module load petsc-3.16.0-gcc-7.3.0-4i7pnvi
# pkgconf@1.8.0%gcc@7.3.0 arch=linux-centos7-broadwell
module load pkgconf-1.8.0-gcc-7.3.0-eoghvsn
# py-mpi4py@3.1.2%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-mpi4py-3.1.2-gcc-7.3.0-jhzm4fm
# py-pip@21.3.1%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-pip-21.3.1-gcc-7.3.0-kqv6ze7
# py-setuptools@59.4.0%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-setuptools-59.4.0-gcc-7.3.0-ysqfh5r
# py-wheel@0.37.0%gcc@7.3.0 arch=linux-centos7-broadwell
module load py-wheel-0.37.0-gcc-7.3.0-ys7zinn
# python@3.8.5%gcc@7.3.0+bz2+ctypes+dbm~debug+ensurepip+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tix~tkinter~ucs4+uuid+zlib patches=0d98e93189bc278fbc37a50ed7f183bd8aaf249a8e1670a465f0db6bb4f8cf87,4c2457325f2b608b1b6a2c63087df8c26e07db3e3d493caf36a56f0ecf6fb768,f2fd060afc4b4618fe8104c4c5d771f36dc55b1db5a4623785a4ea707ec72fb4 arch=linux-centos7-broadwell
module load python-3.8.5-gcc-7.3.0-ahml2fl
# raja@0.14.0%gcc@7.3.0+cuda~examples~exercises~ipo+openmp~rocm~shared~tests build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-broadwell
module load raja-0.14.0-gcc-7.3.0-ui2yfv6
# suite-sparse@5.10.1%gcc@7.3.0~cuda~graphblas~openmp+pic~tbb arch=linux-centos7-broadwell
module load suite-sparse-5.10.1-gcc-7.3.0-opk6y26
# superlu-dist@7.2.0%gcc@7.3.0~cuda~int64~ipo~openmp~rocm+shared build_type=RelWithDebInfo patches=8da9e21f724e8f11a5782960cc0322e12f724bf93cead7df517901a788ea3d61 arch=linux-centos7-broadwell
module load superlu-dist-7.2.0-gcc-7.3.0-5rdw7ha
# texinfo@6.5%gcc@7.3.0 patches=12f6edb0c6b270b8c8dba2ce17998c580db01182d871ee32b7b6e4129bd1d23a,1732115f651cff98989cb0215d8f64da5e0f7911ebf0c13b064920f088f2ffe1 arch=linux-centos7-broadwell
module load texinfo-6.5-gcc-7.3.0-ae334fr
# umpire@6.0.0%gcc@7.3.0+c+cuda~deviceconst~examples~fortran~ipo~numa~openmp~rocm~shared build_type=RelWithDebInfo cuda_arch=60 tests=none arch=linux-centos7-broadwell
module load umpire-6.0.0-gcc-7.3.0-pq7d6om
# zlib@1.2.11%gcc@7.3.0+optimize+pic+shared arch=linux-centos7-broadwell
module load zlib-1.2.11-gcc-7.3.0-43cnepi

# Load system modules
module load gcc/7.3.0
module load openmpi/3.1.3
module load cuda/10.2.89

module load python/miniconda3.8
source /share/apps/python/miniconda3.8/etc/profile.d/conda.sh

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
