source /etc/profile.d/modules.sh
export MY_CLUSTER=marianas

module purge
module use -a /qfs/projects/exasgd/src/cameron-spack/share/spack/modules/linux-centos7-haswell/
module use -a /qfs/projects/exasgd/src/cameron-spack/share/spack/modules/linux-centos7-broadwell/

# Load spack modules
# blt@0.4.1%gcc@7.3.0 arch=linux-centos7-haswell
module load blt-0.4.1-gcc-7.3.0-2u5wdnb
# camp@0.2.2%gcc@7.3.0+cuda~ipo~rocm~tests amdgpu_target=none build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-haswell
module load camp-0.2.2-gcc-7.3.0-em6qqth
# cmake@3.21.4%gcc@7.3.0~doc+ncurses+openssl+ownlibs~qt build_type=Release arch=linux-centos7-haswell
module load cmake-3.21.4-gcc-7.3.0-vh3thpl
# coinhsl@2015.06.23%gcc@7.3.0+blas arch=linux-centos7-haswell
module load coinhsl-2015.06.23-gcc-7.3.0-elcnqge
# cub@1.12.0-rc0%gcc@7.3.0 arch=linux-centos7-haswell
module load cub-1.12.0-rc0-gcc-7.3.0-sy6tp4j
# cuda@10.2.89%gcc@7.3.0~dev arch=linux-centos7-broadwell
module load cuda-10.2.89-gcc-7.3.0-2cpra4f
# diffutils@3.8%gcc@7.3.0 arch=linux-centos7-haswell
module load diffutils-3.8-gcc-7.3.0-sxh6j3c
# hdf5@1.10.7%gcc@7.3.0~cxx~fortran~hl~ipo~java+mpi+shared~szip~threadsafe+tools api=default build_type=RelWithDebInfo arch=linux-centos7-haswell
module load hdf5-1.10.7-gcc-7.3.0-qio2ewe
# hiop@0.5.1%gcc@7.3.0+cuda~deepchecking~ipo~jsrun~kron+mpi+raja~shared+sparse build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-haswell
module load hiop-0.5.1-gcc-7.3.0-pxdqnuo
# hypre@2.22.0%gcc@7.3.0~complex~cuda~debug+fortran~int64~internal-superlu~mixedint+mpi~openmp+shared~superlu-dist~unified-memory cuda_arch=none arch=linux-centos7-haswell
module load hypre-2.22.0-gcc-7.3.0-5kgy7wh
# ipopt@3.12.10%gcc@7.3.0+coinhsl+debug~metis~mumps arch=linux-centos7-haswell
module load ipopt-3.12.10-gcc-7.3.0-i6mp2xz
# libiconv@1.16%gcc@7.3.0 libs=shared,static arch=linux-centos7-haswell
module load libiconv-1.16-gcc-7.3.0-crirgwg
# magma@2.6.1%gcc@7.3.0+cuda+fortran~ipo~rocm+shared amdgpu_target=none build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-haswell
module load magma-2.6.1-gcc-7.3.0-fkzs6jx
# metis@5.1.0%gcc@7.3.0~gdb~int64~real64+shared build_type=Release patches=4991da938c1d3a1d3dea78e49bbebecba00273f98df2a656e38b83d55b281da1,b1225da886605ea558db7ac08dd8054742ea5afe5ed61ad4d0fe7a495b1270d2 arch=linux-centos7-haswell
module load metis-5.1.0-gcc-7.3.0-mp3iiiv
# ncurses@6.2%gcc@7.3.0~symlinks+termlib abi=none arch=linux-centos7-haswell
module load ncurses-6.2-gcc-7.3.0-6u5ymdz
# openblas@0.3.18%gcc@4.8.5~bignuma~consistent_fpcsr~ilp64+locking+pic+shared threads=none arch=linux-centos7-haswell
module load openblas-0.3.18-gcc-4.8.5-iqjleul
# openmpi@3.1.3%gcc@7.3.0~atomics~cuda~cxx~cxx_exceptions+gpfs~internal-hwloc~java~legacylaunchers~lustre~memchecker~pmi~pmix~singularity~sqlite3+static~thread_multiple+vt+wrapper-rpath fabrics=none schedulers=none arch=linux-centos7-haswell
module load openmpi-3.1.3-gcc-7.3.0-jkx4j7w
# openssl@1.1.1l%gcc@7.3.0~docs certs=system arch=linux-centos7-haswell
module load openssl-1.1.1l-gcc-7.3.0-kh4qqnb
# parmetis@4.0.3%gcc@7.3.0~gdb~int64~ipo+shared build_type=RelWithDebInfo patches=4f892531eb0a807eb1b82e683a416d3e35154a455274cf9b162fb02054d11a5b,50ed2081bc939269689789942067c58b3e522c269269a430d5d34c00edbc5870,704b84f7c7444d4372cb59cca6e1209df4ef3b033bc4ee3cf50f369bce972a9d arch=linux-centos7-haswell
module load parmetis-4.0.3-gcc-7.3.0-lj2ynno
# petsc@3.14.6%gcc@7.3.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw~giflib+hdf5~hpddm~hwloc+hypre~int64~jpeg~knl~libpng~libyaml~memkind+metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse+superlu-dist~tetgen~trilinos~valgrind amdgpu_target=none clanguage=C cuda_arch=none arch=linux-centos7-haswell
module load petsc-3.14.6-gcc-7.3.0-mk23bgh
# pkgconf@1.8.0%gcc@7.3.0 arch=linux-centos7-haswell
module load pkgconf-1.8.0-gcc-7.3.0-uyrzmnn
# raja@0.14.0%gcc@7.3.0+cuda~examples~exercises~ipo+openmp~rocm~shared~tests amdgpu_target=none build_type=RelWithDebInfo cuda_arch=60 arch=linux-centos7-haswell
module load raja-0.14.0-gcc-7.3.0-ug4wkba
# superlu-dist@7.1.1%gcc@7.3.0~cuda~int64~ipo~openmp+shared build_type=RelWithDebInfo cuda_arch=none arch=linux-centos7-haswell
module load superlu-dist-7.1.1-gcc-7.3.0-yjwh4el
# umpire@6.0.0%gcc@7.3.0+c+cuda~deviceconst~examples~fortran~ipo~numa~openmp~rocm~shared amdgpu_target=none build_type=RelWithDebInfo cuda_arch=60 tests=none arch=linux-centos7-haswell
module load umpire-6.0.0-gcc-7.3.0-l65sw4w
# zlib@1.2.11%gcc@7.3.0+optimize+pic+shared arch=linux-centos7-haswell
module load zlib-1.2.11-gcc-7.3.0-3m33c3a

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
