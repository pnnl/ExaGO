source /etc/profile.d/modules.sh
export MY_CLUSTER=marianas

module purge
module use -a /qfs/projects/exasgd/src/spack/share/spack/modules/linux-centos7-broadwell/
module use -a /qfs/projects/exasgd/src/spack/share/spack/modules/linux-centos7-x86_64/

# Load spack modules
# autoconf-archive@2019.01.06%gcc@7.3.0 arch=linux-centos7-x86_64
module load exasgd-autoconf-archive/2019.01.06/gcc-7.3.0-5a4lwsi
# blt@0.3.6%gcc@7.3.0 arch=linux-centos7-x86_64
module load exasgd-blt/0.3.6/gcc-7.3.0-glzuyxk
# camp@0.1.0%gcc@7.3.0+cuda~ipo~rocm~tests amdgpu_target=none build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-x86_64
module load exasgd-camp/0.1.0/cuda-10.2.89/gcc-7.3.0-cvqyc5h
# cmake@3.19.6%gcc@7.3.0~doc+ncurses+openssl+ownlibs~qt build_type=Release patches=b48396c0e4f61756248156b6cebe9bc0d7a22228639b47b5aa77c9330588ce88 arch=linux-centos7-x86_64
module load exasgd-cmake/3.19.6/gcc-7.3.0-ihjiqv4
# coinhsl@2015.06.23%gcc@7.3.0+blas arch=linux-centos7-x86_64
module load exasgd-coinhsl/2015.06.23/gcc-7.3.0-zv2yxcl
# cub@1.12.0-rc0%gcc@7.3.0 arch=linux-centos7-x86_64
module load exasgd-cub/1.12.0-rc0/gcc-7.3.0-cscwjgd
# cuda@10.2.89%gcc@7.3.0~dev arch=linux-centos7-x86_64
module load exasgd-cuda/10.2.89/gcc-7.3.0-4ytm53w
# gmp@6.2.1%gcc@7.3.0 arch=linux-centos7-x86_64
module load exasgd-gmp/6.2.1/gcc-7.3.0-sqgj2tv
# hiop@0.4.6%gcc@7.3.0+cuda~deepchecking~ipo~jsrun+kron+mpi+raja+shared+sparse build_type=Release cuda_arch=60 arch=linux-centos7-x86_64
module load exasgd-hiop/0.4.6/cuda-10.2.89/openmpi-3.1.3/gcc-7.3.0-jeuhgas
# ipopt@3.12.10%gcc@7.3.0+coinhsl~debug~metis~mumps arch=linux-centos7-x86_64
module load exasgd-ipopt/3.12.10/gcc-7.3.0-pyjg6pd
# magma@2.6.1%gcc@7.3.0+cuda+fortran~ipo~rocm+shared amdgpu_target=none build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-x86_64
module load exasgd-magma/2.6.1/cuda-10.2.89/gcc-7.3.0-mb7xpns
# metis@5.1.0%gcc@7.3.0~gdb~int64~real64+shared build_type=Release patches=4991da938c1d3a1d3dea78e49bbebecba00273f98df2a656e38b83d55b281da1,b1225da886605ea558db7ac08dd8054742ea5afe5ed61ad4d0fe7a495b1270d2 arch=linux-centos7-x86_64
module load exasgd-metis/5.1.0/gcc-7.3.0-t2i65bz
# mpfr@4.1.0%gcc@7.3.0 arch=linux-centos7-x86_64
module load exasgd-mpfr/4.1.0/gcc-7.3.0-3sbafnc
# openblas@0.3.15%gcc@4.8.5~bignuma~consistent_fpcsr~ilp64+locking+pic+shared threads=none arch=linux-centos7-x86_64
module load exasgd-openblas/0.3.15/gcc-4.8.5-hnfyfu5
# openmpi@3.1.3%gcc@7.3.0~atomics~cuda~cxx~cxx_exceptions+gpfs~internal-hwloc~java~legacylaunchers~lustre~memchecker~pmi~singularity~sqlite3+static~thread_multiple+vt+wrapper-rpath fabrics=none schedulers=none arch=linux-centos7-x86_64
module load exasgd-openmpi/3.1.3/gcc-7.3.0-ccr6unl
# perl@5.26.0%gcc@4.8.5+cpanm~shared~threads patches=0eac10ed90aeb0459ad8851f88081d439a4e41978e586ec743069e8b059370ac arch=linux-centos7-x86_64
module load exasgd-perl/5.26.0/gcc-4.8.5-5cziri4
# raja@0.12.1%gcc@7.3.0+cuda~examples~exercises~ipo+openmp~rocm+shared~tests amdgpu_target=none build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-x86_64
module load exasgd-raja/0.12.1/cuda-10.2.89/gcc-7.3.0-wskztqm
# suite-sparse@5.10.1%gcc@7.3.0~cuda~openmp+pic~tbb arch=linux-centos7-x86_64
module load exasgd-suite-sparse/5.10.1/gcc-7.3.0-vgrw3e7
# texinfo@6.5%gcc@7.3.0 patches=12f6edb0c6b270b8c8dba2ce17998c580db01182d871ee32b7b6e4129bd1d23a,1732115f651cff98989cb0215d8f64da5e0f7911ebf0c13b064920f088f2ffe1 arch=linux-centos7-x86_64
module load exasgd-texinfo/6.5/gcc-7.3.0-yoz32zy
# umpire@4.1.2%gcc@7.3.0+c+cuda~deviceconst~examples~fortran~ipo~numa+openmp~rocm~shared amdgpu_target=none build_type=RelWithDebInfo cuda_arch=none patches=135bbc7d2f371531f432672b115ac0a407968aabfffc5b8a941db9b493dbf81f,7d912d31cd293df005ba74cb96c6f3e32dc3d84afff49b14509714283693db08 tests=none arch=linux-centos7-x86_64
module load exasgd-umpire/4.1.2/cuda-10.2.89/gcc-7.3.0-4bfftfv
# petsc@3.14.6%gcc@7.3.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw~giflib+hdf5~hwloc+hypre~int64~jpeg~knl~libpng~libyaml~memkind+metis~mkl-pardiso~moab~mpfr+mpi~mumps~openmp~p4est~ptscotch~random123~rocm~saws+shared~suite-sparse+superlu-dist~trilinos~valgrind amdgpu_target=none clanguage=C cuda_arch=none arch=linux-centos7-x86_64
module load exasgd-petsc/3.14.6/openmpi-3.1.3/gcc-7.3.0-54tno2h

# Load system modules
module load gcc/7.3.0
module load openmpi/3.1.3
module load cuda/10.2.89
module load python/anaconda3.2019.3

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
