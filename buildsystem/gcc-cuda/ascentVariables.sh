export MY_CLUSTER=ascent

module purge

module use -a /gpfs/wolf/proj-shared/csc359/src/spack/share/spack/modules/linux-rhel7-power9le

# Load spack modules
# arpack-ng@3.8.0%gcc@7.4.0+mpi+shared arch=linux-rhel7-power9le
module load exasgd-arpack-ng/3.8.0/spectrum-mpi-10.3.1.2-20200121/gcc-7.4.0-e6d5zbw
# automake@1.16.3%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-automake/1.16.3/gcc-7.4.0-qprf3ee
# berkeley-db@18.1.40%gcc@7.4.0+cxx~docs+stl patches=b231fcc4d5cff05e5c3a4814f6a5af0e9a966428dc2176540d2c05aff41de522 arch=linux-rhel7-power9le
module load exasgd-berkeley-db/18.1.40/gcc-7.4.0-hn7z7sb
# blaspp@2021.04.01%gcc@7.4.0~cuda~ipo+openmp~rocm+shared amdgpu_target=none build_type=RelWithDebInfo cuda_arch=none arch=linux-rhel7-power9le
module load exasgd-blaspp/2021.04.01/gcc-7.4.0-acu7d4s
# butterflypack@1.2.1%gcc@7.4.0~ipo+shared build_type=RelWithDebInfo arch=linux-rhel7-power9le
module load exasgd-butterflypack/1.2.1/spectrum-mpi-10.3.1.2-20200121/gcc-7.4.0-bavjqs4
# coinhsl@2015.06.23%gcc@7.4.0+blas arch=linux-rhel7-power9le
module load exasgd-coinhsl/2015.06.23/gcc-7.4.0-ubfbspm
# gdbm@1.19%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-gdbm/1.19/gcc-7.4.0-v4st5ql
# gmp@6.2.1%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-gmp/6.2.1/gcc-7.4.0-6svtfvr
# hdf5@1.10.7%gcc@7.4.0~cxx~debug~fortran~hl~java+mpi+pic+shared~szip~threadsafe api=none arch=linux-rhel7-power9le
module load exasgd-hdf5/1.10.7/spectrum-mpi-10.3.1.2-20200121/gcc-7.4.0-take6gf
# hiop@0.4.1%gcc@7.4.0+cuda+deepchecking+gpu~hip~ipo~jsrun+kron+mpi+raja+shared+sparse~srun amdgpu_target=none build_type=Release cuda_arch=70 arch=linux-rhel7-power9le
module load exasgd-hiop/0.4.1/cuda-11.0.194/spectrum-mpi-10.3.1.2-20200121/gcc-7.4.0-4xarlsx
# hypre@2.20.0%gcc@7.4.0~complex~cuda~debug~int64~internal-superlu~mixedint+mpi~openmp+shared~superlu-dist~unified-memory cuda_arch=none patches=6e3336b1d62155f6350dfe42b0f9ea25d4fa0af60c7e540959139deb93a26059 arch=linux-rhel7-power9le
module load exasgd-hypre/2.20.0/spectrum-mpi-10.3.1.2-20200121/gcc-7.4.0-qcoeii7
# ipopt@3.12.10%gcc@7.4.0+coinhsl~debug~mumps arch=linux-rhel7-power9le
module load exasgd-ipopt/3.12.10/gcc-7.4.0-xs3yw6t
# lapackpp@2021.04.00%gcc@7.4.0~ipo+shared build_type=RelWithDebInfo arch=linux-rhel7-power9le
module load exasgd-lapackpp/2021.04.00/gcc-7.4.0-ikfu5fv
# libiconv@1.16%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-libiconv/1.16/gcc-7.4.0-z5ga2ok
# libsigsegv@2.13%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-libsigsegv/2.13/gcc-7.4.0-garv4jn
# libtool@2.4.6%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-libtool/2.4.6/gcc-7.4.0-lbasl6y
# magma@2.5.4%gcc@7.4.0+cuda+fortran~ipo+shared build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel7-power9le
module load exasgd-magma/2.5.4/cuda-11.0.194/gcc-7.4.0-pdoytfx
# metis@5.1.0%gcc@7.4.0~gdb~int64~real64+shared build_type=Release patches=4991da938c1d3a1d3dea78e49bbebecba00273f98df2a656e38b83d55b281da1,b1225da886605ea558db7ac08dd8054742ea5afe5ed61ad4d0fe7a495b1270d2 arch=linux-rhel7-power9le
module load exasgd-metis/5.1.0/gcc-7.4.0-7cjo5kb
# mpfr@4.1.0%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-mpfr/4.1.0/gcc-7.4.0-bkusb67
# ncurses@6.2%gcc@7.4.0~symlinks+termlib abi=none arch=linux-rhel7-power9le
module load exasgd-ncurses/6.2/gcc-7.4.0-rzufnf5
# netlib-scalapack@2.1.0%gcc@7.4.0~ipo~pic+shared build_type=Release patches=1c9ce5fee1451a08c2de3cc87f446aeda0b818ebbce4ad0d980ddf2f2a0b2dc4,f2baedde688ffe4c20943c334f580eb298e04d6f35c86b90a1f4e8cb7ae344a2 arch=linux-rhel7-power9le
module load exasgd-netlib-scalapack/2.1.0/spectrum-mpi-10.3.1.2-20200121/gcc-7.4.0-tyguivd
# openblas@0.3.10%gcc@7.4.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared patches=865703b4f405543bbd583413fdeff2226dfda908be33639276c06e5aa7ae2cae threads=none arch=linux-rhel7-power9le
module load exasgd-openblas/0.3.10/gcc-7.4.0-uynupjj
# parmetis@4.0.3%gcc@7.4.0~gdb~int64~ipo+shared build_type=RelWithDebInfo patches=4f892531eb0a807eb1b82e683a416d3e35154a455274cf9b162fb02054d11a5b,50ed2081bc939269689789942067c58b3e522c269269a430d5d34c00edbc5870,704b84f7c7444d4372cb59cca6e1209df4ef3b033bc4ee3cf50f369bce972a9d arch=linux-rhel7-power9le
module load exasgd-parmetis/4.0.3/spectrum-mpi-10.3.1.2-20200121/gcc-7.4.0-sdfju5f
# perl@5.34.0%gcc@7.4.0+cpanm+shared+threads arch=linux-rhel7-power9le
module load exasgd-perl/5.34.0/gcc-7.4.0-fh5ekqg
# petsc@3.14.6%gcc@7.4.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw~giflib+hdf5~hip+hypre~int64~jpeg~knl~libpng~libyaml~memkind+metis~mkl-pardiso~moab~mpfr+mpi~mumps~p4est~ptscotch~random123~saws+shared~suite-sparse+superlu-dist~trilinos~valgrind amdgpu_target=none clanguage=C arch=linux-rhel7-power9le
module load exasgd-petsc/3.14.6/spectrum-mpi-10.3.1.2-20200121/gcc-7.4.0-lgeqldi
# raja@0.13.0%gcc@7.4.0+cuda~examples~exercises~hip~ipo+openmp+shared~tests amdgpu_target=none build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel7-power9le
module load exasgd-raja/0.13.0/cuda-11.0.194/gcc-7.4.0-hhpqk3k
# readline@8.1%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-readline/8.1/gcc-7.4.0-vqamj6w
# slate@2021.05.02%gcc@7.4.0~cuda~ipo+mpi+openmp~rocm+shared amdgpu_target=none build_type=RelWithDebInfo cuda_arch=none arch=linux-rhel7-power9le
module load exasgd-slate/2021.05.02/spectrum-mpi-10.3.1.2-20200121/gcc-7.4.0-7xtdafw
# strumpack@master%gcc@7.4.0~build_dev_tests~build_tests+butterflypack+c_interface~count_flops~cuda~ipo+mpi~openmp+parmetis~rocm~scotch~shared+slate~task_timers+zfp amdgpu_target=none build_type=RelWithDebInfo cuda_arch=none arch=linux-rhel7-power9le
module load exasgd-strumpack/master/spectrum-mpi-10.3.1.2-20200121/gcc-7.4.0-d7nw5iz
# suite-sparse@5.10.1%gcc@7.4.0~cuda~openmp+pic~tbb arch=linux-rhel7-power9le
module load exasgd-suite-sparse/5.10.1/gcc-7.4.0-t5u7fh4
# superlu-dist@6.4.0%gcc@7.4.0~cuda~int64~ipo~openmp+shared build_type=RelWithDebInfo cuda_arch=none arch=linux-rhel7-power9le
module load exasgd-superlu-dist/6.4.0/spectrum-mpi-10.3.1.2-20200121/gcc-7.4.0-gd5qioe
# texinfo@6.5%gcc@7.4.0 patches=12f6edb0c6b270b8c8dba2ce17998c580db01182d871ee32b7b6e4129bd1d23a,1732115f651cff98989cb0215d8f64da5e0f7911ebf0c13b064920f088f2ffe1 arch=linux-rhel7-power9le
module load exasgd-texinfo/6.5/gcc-7.4.0-w4cx3wo
# umpire@4.1.2%gcc@7.4.0+c+cuda~deviceconst~examples~fortran~hip~ipo~numa+openmp~shared amdgpu_target=none build_type=RelWithDebInfo cuda_arch=none patches=7d912d31cd293df005ba74cb96c6f3e32dc3d84afff49b14509714283693db08 tests=none arch=linux-rhel7-power9le
module load exasgd-umpire/4.1.2/cuda-11.0.194/gcc-7.4.0-gii7nig
# zfp@0.5.5%gcc@7.4.0~aligned~c~cuda~fasthash~fortran~ipo~openmp~profile~python+shared~strided~twoway bsws=64 build_type=RelWithDebInfo cuda_arch=none arch=linux-rhel7-power9le
module load exasgd-zfp/0.5.5/gcc-7.4.0-av4beua
# zlib@1.2.11%gcc@7.4.0+optimize+pic+shared arch=linux-rhel7-power9le
module load exasgd-zlib/1.2.11/gcc-7.4.0-fkwiwuh

# Load system modules
module load cuda/11.0.2
module load gcc/7.4.0
module load spectrum-mpi/10.3.1.2-20200121
module load cmake/3.18.2

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
