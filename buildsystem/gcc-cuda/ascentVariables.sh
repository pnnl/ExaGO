export MY_CLUSTER=ascent

module purge

module use -a /gpfs/wolf/proj-shared/csc359/src/cameron-spack/share/spack/modules/linux-rhel8-power9le

# Load spack modules
# autoconf@2.69%gcc@9.1.0 patches=35c4492,7793209,a49dd5b arch=linux-rhel8-power9le
module load autoconf-2.69-gcc-9.1.0-vnxzsnr
# autoconf-archive@2019.01.06%gcc@9.1.0 arch=linux-rhel8-power9le
module load autoconf-archive-2019.01.06-gcc-9.1.0-2kjmyyv
# automake@1.16.5%gcc@9.1.0 arch=linux-rhel8-power9le
module load automake-1.16.5-gcc-9.1.0-x5ndgg2
# blt@0.4.1%gcc@9.1.0 arch=linux-rhel8-power9le
module load blt-0.4.1-gcc-9.1.0-pvxpp36
# bzip2@1.0.8%gcc@9.1.0~debug~pic+shared arch=linux-rhel8-power9le
module load bzip2-1.0.8-gcc-9.1.0-nyx3bis
# camp@0.2.2%gcc@9.1.0+cuda~ipo~rocm~tests build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel8-power9le
module load camp-0.2.2-gcc-9.1.0-kd2joie
# coinhsl@2015.06.23%gcc@9.1.0+blas arch=linux-rhel8-power9le
module load coinhsl-2015.06.23-gcc-9.1.0-go24cay
# cub@1.12.0-rc0%gcc@9.1.0 arch=linux-rhel8-power9le
module load cub-1.12.0-rc0-gcc-9.1.0-ipbwuis
# diffutils@3.8%gcc@9.1.0 arch=linux-rhel8-power9le
module load diffutils-3.8-gcc-9.1.0-ob435zm
# expat@2.4.6%gcc@9.1.0+libbsd arch=linux-rhel8-power9le
module load expat-2.4.6-gcc-9.1.0-aqiqfga
# gdbm@1.23%gcc@9.1.0 arch=linux-rhel8-power9le
module load gdbm-1.23-gcc-9.1.0-jfuuhq4
# gettext@0.21%gcc@9.1.0+bzip2+curses+git~libunistring+libxml2+tar+xz arch=linux-rhel8-power9le
module load gettext-0.21-gcc-9.1.0-n6twic7
# gmp@6.2.1%gcc@9.1.0 arch=linux-rhel8-power9le
module load gmp-6.2.1-gcc-9.1.0-76za3mf
# gnuconfig@2021-08-14%gcc@9.1.0 arch=linux-rhel8-power9le
module load gnuconfig-2021-08-14-gcc-9.1.0-wt2yuir
# hdf5@1.12.1%gcc@9.1.0~cxx~fortran~hl~ipo~java+mpi+shared~szip~threadsafe+tools api=default build_type=RelWithDebInfo patches=ee351eb arch=linux-rhel8-power9le
module load hdf5-1.12.1-gcc-9.1.0-glesldy
# hiop@0.5.4%gcc@9.1.0+cuda+deepchecking~ipo~jsrun+kron+mpi+raja~rocm~shared+sparse build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel8-power9le
module load hiop-0.5.4-gcc-9.1.0-ymrqy2o
# hypre@2.24.0%gcc@9.1.0~complex~cuda~debug+fortran~gptune~int64~internal-superlu~mixedint+mpi~openmp+shared~superlu-dist~unified-memory arch=linux-rhel8-power9le
module load hypre-2.24.0-gcc-9.1.0-ld2fu6w
# ipopt@3.12.10%gcc@9.1.0+coinhsl~debug~metis~mumps arch=linux-rhel8-power9le
module load ipopt-3.12.10-gcc-9.1.0-jpxdtwv
# libbsd@0.11.5%gcc@9.1.0 arch=linux-rhel8-power9le
module load libbsd-0.11.5-gcc-9.1.0-g5gbwir
# libffi@3.4.2%gcc@9.1.0 arch=linux-rhel8-power9le
module load libffi-3.4.2-gcc-9.1.0-b75lkjq
# libiconv@1.16%gcc@9.1.0 libs=shared,static arch=linux-rhel8-power9le
module load libiconv-1.16-gcc-9.1.0-5upkvmm
# libmd@1.0.4%gcc@9.1.0 arch=linux-rhel8-power9le
module load libmd-1.0.4-gcc-9.1.0-7tohrlx
# libsigsegv@2.13%gcc@9.1.0 arch=linux-rhel8-power9le
module load libsigsegv-2.13-gcc-9.1.0-lrapsw3
# libtool@2.4.6%gcc@9.1.0 arch=linux-rhel8-power9le
module load libtool-2.4.6-gcc-9.1.0-c7ioa2h
# libxml2@2.9.12%gcc@9.1.0~python arch=linux-rhel8-power9le
module load libxml2-2.9.12-gcc-9.1.0-a2kj54a
# m4@1.4.19%gcc@9.1.0+sigsegv patches=9dc5fbd,bfdffa7 arch=linux-rhel8-power9le
module load m4-1.4.19-gcc-9.1.0-aw4chza
# magma@2.6.2rc1%gcc@9.1.0+cuda+fortran~ipo~rocm+shared build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel8-power9le
module load magma-2.6.2rc1-gcc-9.1.0-w7s57k6
# metis@5.1.0%gcc@9.1.0~gdb~int64~real64+shared build_type=Release patches=4991da9,b1225da arch=linux-rhel8-power9le
module load metis-5.1.0-gcc-9.1.0-sg2ca2t
# mpfr@4.1.0%gcc@9.1.0 arch=linux-rhel8-power9le
module load mpfr-4.1.0-gcc-9.1.0-2qjplyw
# ncurses@6.2%gcc@9.1.0~symlinks+termlib abi=none arch=linux-rhel8-power9le
module load ncurses-6.2-gcc-9.1.0-aa4mmvm
# openblas@0.3.19%gcc@9.1.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared symbol_suffix=none threads=none arch=linux-rhel8-power9le
module load openblas-0.3.19-gcc-9.1.0-c6nslyv
# openssl@1.1.1m%gcc@9.1.0~docs certs=system arch=linux-rhel8-power9le
module load openssl-1.1.1m-gcc-9.1.0-u6h6mbd
# parmetis@4.0.3%gcc@9.1.0~gdb~int64~ipo+shared build_type=RelWithDebInfo patches=4f89253,50ed208,704b84f arch=linux-rhel8-power9le
module load parmetis-4.0.3-gcc-9.1.0-nhshpa6
# perl@5.30.1%gcc@9.1.0+cpanm+shared+threads arch=linux-rhel8-power9le
module load perl-5.30.1-gcc-9.1.0-qmsmncp
# petsc@3.16.0%gcc@9.1.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw~giflib+hdf5~hpddm~hwloc+hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind+metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse+superlu-dist~tetgen~trilinos~valgrind clanguage=C arch=linux-rhel8-power9le
module load petsc-3.16.0-gcc-9.1.0-lnzdun7
# pkgconf@1.8.0%gcc@9.1.0 arch=linux-rhel8-power9le
module load pkgconf-1.8.0-gcc-9.1.0-ne6atqj
# py-mpi4py@3.1.2%gcc@9.1.0 arch=linux-rhel8-power9le
module load py-mpi4py-3.1.2-gcc-9.1.0-6v77rgc
# py-pip@21.3.1%gcc@9.1.0 arch=linux-rhel8-power9le
module load py-pip-21.3.1-gcc-9.1.0-dezo3ft
# py-setuptools@59.4.0%gcc@9.1.0 arch=linux-rhel8-power9le
module load py-setuptools-59.4.0-gcc-9.1.0-zqlhxba
# py-wheel@0.37.0%gcc@9.1.0 arch=linux-rhel8-power9le
module load py-wheel-0.37.0-gcc-9.1.0-usjkar2
# python@3.9.10%gcc@9.1.0+bz2+ctypes+dbm~debug+ensurepip+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tix~tkinter~ucs4+uuid+zlib patches=0d98e93,4c24573,f2fd060 arch=linux-rhel8-power9le
module load python-3.9.10-gcc-9.1.0-u5qrmqv
# raja@0.14.0%gcc@9.1.0+cuda~examples~exercises~ipo+openmp~rocm+shared~tests build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel8-power9le
module load raja-0.14.0-gcc-9.1.0-ili5h35
# readline@8.1%gcc@9.1.0 arch=linux-rhel8-power9le
module load readline-8.1-gcc-9.1.0-bp2hipf
# sqlite@3.37.2%gcc@9.1.0+column_metadata+dynamic_extensions+fts~functions+rtree arch=linux-rhel8-power9le
module load sqlite-3.37.2-gcc-9.1.0-doittuo
# suite-sparse@5.8.1%gcc@9.1.0~cuda~graphblas~openmp+pic~tbb arch=linux-rhel8-power9le
module load suite-sparse-5.8.1-gcc-9.1.0-r6zyrff
# superlu-dist@7.2.0%gcc@9.1.0~cuda~int64~ipo~openmp~rocm+shared build_type=RelWithDebInfo patches=8da9e21 arch=linux-rhel8-power9le
module load superlu-dist-7.2.0-gcc-9.1.0-isebdzh
# tar@1.34%gcc@9.1.0 arch=linux-rhel8-power9le
module load tar-1.34-gcc-9.1.0-pi3bs2h
# texinfo@6.5%gcc@9.1.0 patches=12f6edb,1732115 arch=linux-rhel8-power9le
module load texinfo-6.5-gcc-9.1.0-jkkfv4q
# umpire@6.0.0%gcc@9.1.0+c+cuda~deviceconst~examples~fortran~ipo~numa+openmp~rocm~shared build_type=RelWithDebInfo cuda_arch=70 tests=none arch=linux-rhel8-power9le
module load umpire-6.0.0-gcc-9.1.0-bosktbw
# util-linux-uuid@2.37.4%gcc@9.1.0 arch=linux-rhel8-power9le
module load util-linux-uuid-2.37.4-gcc-9.1.0-fzig4jl
# xz@5.2.5%gcc@9.1.0~pic libs=shared,static arch=linux-rhel8-power9le
module load xz-5.2.5-gcc-9.1.0-gkiw6gn
# zlib@1.2.11%gcc@9.1.0+optimize+pic+shared arch=linux-rhel8-power9le
module load zlib-1.2.11-gcc-9.1.0-rs4ceke

# Load system modules
module load cuda/11.4.2
module load gcc/9.1.0
module load spectrum-mpi/10.4.0.3-20210112
module load cmake/3.22.2

export MY_PETSC_DIR=$PETSC_DIR
export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DCMAKE_CUDA_ARCHITECTURES=70"
export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DEXAGO_CTEST_LAUNCH_COMMAND='jsrun -g 1'"

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
