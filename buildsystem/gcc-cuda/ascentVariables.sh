export MY_CLUSTER=ascent

module purge

module use -a /gpfs/wolf/proj-shared/csc359/src/spack/share/spack/modules/linux-rhel7-power9le

# Load spack modules
# arpack-ng@3.8.0%gcc@7.4.0+mpi+shared arch=linux-rhel7-power9le
module load exasgd-arpack-ng/3.8.0/spectrum-mpi-10.3.1.2-20200121/gcc-7.4.0-kxbr4jo
# autoconf@2.69%gcc@7.4.0 patches=35c449281546376449766f92d49fc121ca50e330e60fefcfc9be2af3253082c2,7793209b33013dc0f81208718c68440c5aae80e7a1c4b8d336e382525af791a7,a49dd5bac3b62daa0ff688ab4d508d71dbd2f4f8d7e2a02321926346161bf3ee arch=linux-rhel7-power9le
module load exasgd-autoconf/2.69/gcc-7.4.0-6vut5vy
# autoconf-archive@2019.01.06%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-autoconf-archive/2019.01.06/gcc-7.4.0-nn453cx
# automake@1.16.3%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-automake/1.16.3/gcc-7.4.0-hvhod3j
# berkeley-db@18.1.40%gcc@7.4.0+cxx~docs+stl patches=b231fcc4d5cff05e5c3a4814f6a5af0e9a966428dc2176540d2c05aff41de522 arch=linux-rhel7-power9le
module load exasgd-berkeley-db/18.1.40/gcc-7.4.0-ic4cqif
# blaspp@2021.04.01%gcc@7.4.0+cuda~ipo+openmp~rocm+shared amdgpu_target=none build_type=RelWithDebInfo cuda_arch=none arch=linux-rhel7-power9le
module load exasgd-blaspp/2021.04.01/cuda-11.0.194/gcc-7.4.0-6pxqmbi
# blt@0.4.1%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-blt/0.4.1/gcc-7.4.0-rb4qbqh
# butterflypack@2.0.0%gcc@7.4.0~ipo+shared build_type=RelWithDebInfo arch=linux-rhel7-power9le
module load exasgd-butterflypack/2.0.0/spectrum-mpi-10.3.1.2-20200121/gcc-7.4.0-efdckib
# bzip2@1.0.8%gcc@7.4.0~debug~pic+shared arch=linux-rhel7-power9le
module load exasgd-bzip2/1.0.8/gcc-7.4.0-jty62q7
# camp@0.2.2%gcc@7.4.0+cuda~ipo~rocm~tests amdgpu_target=none build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel7-power9le
module load exasgd-camp/0.2.2/cuda-11.0.194/gcc-7.4.0-xqz5l7f
# cmake@3.18.2%gcc@7.4.0~doc+ncurses+openssl+ownlibs~qt build_type=Release patches=bf695e3febb222da2ed94b3beea600650e4318975da90e4a71d6f31a6d5d8c3d arch=linux-rhel7-power9le
module load exasgd-cmake/3.18.2/gcc-7.4.0-m42tpk4
# coinhsl@2015.06.23%gcc@7.4.0+blas arch=linux-rhel7-power9le
module load exasgd-coinhsl/2015.06.23/gcc-7.4.0-ka7y5hb
# cub@1.12.0-rc0%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-cub/1.12.0-rc0/gcc-7.4.0-iwyj63t
# cuda@11.0.194%gcc@7.4.0~dev arch=linux-rhel7-power9le
module load exasgd-cuda/11.0.194/gcc-7.4.0-e5nl5fx
# diffutils@3.8%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-diffutils/3.8/gcc-7.4.0-cy55hsj
# gdbm@1.19%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-gdbm/1.19/gcc-7.4.0-ahdwucz
# gmp@6.2.1%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-gmp/6.2.1/gcc-7.4.0-ur2a3rb
# gnuconfig@2021-08-14%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-gnuconfig/2021-08-14/gcc-7.4.0-qr6nxuq
# hiop@0.5.1%gcc@7.4.0+cuda+deepchecking~ipo~jsrun+kron+mpi+raja~shared+sparse build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel7-power9le
module load exasgd-hiop/0.5.1/cuda-11.0.194/spectrum-mpi-10.3.1.2-20200121/gcc-7.4.0-nnnhkv2
# ipopt@3.12.10%gcc@7.4.0+coinhsl~debug~metis~mumps arch=linux-rhel7-power9le
module load exasgd-ipopt/3.12.10/gcc-7.4.0-yma4gh2
# lapackpp@2021.04.00%gcc@7.4.0~ipo+shared build_type=RelWithDebInfo arch=linux-rhel7-power9le
module load exasgd-lapackpp/2021.04.00/cuda-11.0.194/gcc-7.4.0-qhn44r6
# libiconv@1.16%gcc@7.4.0 libs=shared,static arch=linux-rhel7-power9le
module load exasgd-libiconv/1.16/gcc-7.4.0-idqno7d
# libsigsegv@2.13%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-libsigsegv/2.13/gcc-7.4.0-cbn4dja
# libtool@2.4.6%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-libtool/2.4.6/gcc-7.4.0-x5h54ly
# m4@1.4.19%gcc@7.4.0+sigsegv patches=9dc5fbd0d5cb1037ab1e6d0ecc74a30df218d0a94bdd5a02759a97f62daca573,bfdffa7c2eb01021d5849b36972c069693654ad826c1a20b53534009a4ec7a89 arch=linux-rhel7-power9le
module load exasgd-m4/1.4.19/gcc-7.4.0-nrrlksm
# magma@2.6.1%gcc@7.4.0+cuda+fortran~ipo~rocm+shared amdgpu_target=none build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel7-power9le
module load exasgd-magma/2.6.1/cuda-11.0.194/gcc-7.4.0-5ypmigk
# metis@5.1.0%gcc@7.4.0~gdb~int64~real64+shared build_type=Release patches=4991da938c1d3a1d3dea78e49bbebecba00273f98df2a656e38b83d55b281da1,b1225da886605ea558db7ac08dd8054742ea5afe5ed61ad4d0fe7a495b1270d2 arch=linux-rhel7-power9le
module load exasgd-metis/5.1.0/gcc-7.4.0-shhhyku
# mpfr@4.1.0%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-mpfr/4.1.0/gcc-7.4.0-33dvnf2
# ncurses@6.2%gcc@7.4.0~symlinks+termlib abi=none arch=linux-rhel7-power9le
module load exasgd-ncurses/6.2/gcc-7.4.0-kqhmmpv
# netlib-scalapack@2.1.0%gcc@7.4.0~ipo~pic+shared build_type=Release patches=1c9ce5fee1451a08c2de3cc87f446aeda0b818ebbce4ad0d980ddf2f2a0b2dc4,f2baedde688ffe4c20943c334f580eb298e04d6f35c86b90a1f4e8cb7ae344a2 arch=linux-rhel7-power9le
module load exasgd-netlib-scalapack/2.1.0/spectrum-mpi-10.3.1.2-20200121/gcc-7.4.0-j3a6l7a
# openblas@0.3.18%gcc@7.4.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared threads=openmp arch=linux-rhel7-power9le
module load exasgd-openblas/0.3.18/gcc-7.4.0-hd2vahd
# parmetis@4.0.3%gcc@7.4.0~gdb~int64~ipo+shared build_type=RelWithDebInfo patches=4f892531eb0a807eb1b82e683a416d3e35154a455274cf9b162fb02054d11a5b,50ed2081bc939269689789942067c58b3e522c269269a430d5d34c00edbc5870,704b84f7c7444d4372cb59cca6e1209df4ef3b033bc4ee3cf50f369bce972a9d arch=linux-rhel7-power9le
module load exasgd-parmetis/4.0.3/spectrum-mpi-10.3.1.2-20200121/gcc-7.4.0-g5uhy5m
# perl@5.34.0%gcc@7.4.0+cpanm+shared+threads arch=linux-rhel7-power9le
module load exasgd-perl/5.34.0/gcc-7.4.0-h45ivzd
# pkgconf@1.8.0%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-pkgconf/1.8.0/gcc-7.4.0-jfmmybn
# raja@0.14.0%gcc@7.4.0+cuda~examples~exercises~ipo+openmp~rocm+shared~tests amdgpu_target=none build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel7-power9le
module load exasgd-raja/0.14.0/cuda-11.0.194/gcc-7.4.0-pnjrr6w
# readline@8.1%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-readline/8.1/gcc-7.4.0-cszha3o
# slate@2021.05.02%gcc@7.4.0+cuda~ipo+mpi+openmp~rocm+shared amdgpu_target=none build_type=RelWithDebInfo cuda_arch=none arch=linux-rhel7-power9le
module load exasgd-slate/2021.05.02/cuda-11.0.194/spectrum-mpi-10.3.1.2-20200121/gcc-7.4.0-rcwazzh
# spectrum-mpi@10.3.1.2-20200121%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-spectrum-mpi/10.3.1.2-20200121/gcc-7.4.0-nqvpc62
# strumpack@6.1.0%gcc@7.4.0+butterflypack+c_interface~count_flops+cuda~ipo+mpi+openmp+parmetis~rocm~scotch+shared+slate~task_timers+zfp amdgpu_target=none build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel7-power9le
module load exasgd-strumpack/6.1.0/cuda-11.0.194/spectrum-mpi-10.3.1.2-20200121/gcc-7.4.0-cqtleae
# suite-sparse@5.8.1%gcc@7.4.0~cuda~openmp+pic~tbb arch=linux-rhel7-power9le
module load exasgd-suite-sparse/5.8.1/gcc-7.4.0-tu2ka4p
# texinfo@6.5%gcc@7.4.0 patches=12f6edb0c6b270b8c8dba2ce17998c580db01182d871ee32b7b6e4129bd1d23a,1732115f651cff98989cb0215d8f64da5e0f7911ebf0c13b064920f088f2ffe1 arch=linux-rhel7-power9le
module load exasgd-texinfo/6.5/gcc-7.4.0-do3jrb5
# umpire@6.0.0%gcc@7.4.0+c+cuda~deviceconst~examples~fortran~ipo~numa+openmp~rocm~shared amdgpu_target=none build_type=RelWithDebInfo cuda_arch=70 tests=none arch=linux-rhel7-power9le
module load exasgd-umpire/6.0.0/cuda-11.0.194/gcc-7.4.0-7ixxc5e
# zfp@0.5.5%gcc@7.4.0~aligned~c~cuda~fasthash~fortran~ipo~openmp~profile~python+shared~strided~twoway bsws=64 build_type=RelWithDebInfo cuda_arch=none arch=linux-rhel7-power9le
module load exasgd-zfp/0.5.5/gcc-7.4.0-hiozbgg
# zlib@1.2.11%gcc@7.4.0+optimize+pic+shared arch=linux-rhel7-power9le
module load exasgd-zlib/1.2.11/gcc-7.4.0-vnk3szs
# petsc@3.14.6%gcc@7.4.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw~giflib+hdf5~hip+hypre~int64~jpeg~knl~libpng~libyaml~memkind+metis~mkl-pardiso~moa
module load exasgd-petsc/3.14.6/spectrum-mpi-10.3.1.2-20200121/gcc-7.4.0-lgeqldi

# Load system modules
module load cuda/11.0.2
module load gcc/7.4.0
module load spectrum-mpi/10.3.1.2-20200121
module load cmake/3.18.2
module load python/3.7.0

export MY_PETSC_DIR=$PETSC_DIR
export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DCMAKE_CUDA_ARCHITECTURES=70"
export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DEXAGO_TEST_WITH_BSUB=ON"

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

export CC=/sw/ascent/gcc/7.4.0/bin/gcc
export CXX=/sw/ascent/gcc/7.4.0/bin/g++
export FC=/sw/ascent/gcc/7.4.0/bin/gfortran
