source /etc/profile.d/modules.sh

export OMP_CANCELLATION=true
export OMP_PROC_BIND=true
export OMPI_MCA_pml="ucx"
export UCX_NET_DEVICES=mlx5_1:1,mlx5_3:1
export MY_CLUSTER=newell

module purge
module use -a /qfs/projects/exasgd/src/jaelyn-spack/spack/share/spack/modules/linux-rhel7-power9le/

# blt@0.4.1%gcc@7.4.0 arch=linux-rhel7-power9le
module load blt-0.4.1-gcc-7.4.0-vbugf3i
# camp@0.2.2%gcc@7.4.0+cuda~ipo~rocm~tests build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel7-power9le
module load camp-0.2.2-gcc-7.4.0-tdeemi3
# cmake@3.22.2%gcc@7.4.0~doc+ncurses+openssl+ownlibs~qt build_type=Release arch=linux-rhel7-power9le
module load cmake-3.22.2-gcc-7.4.0-mt7gms5
# coinhsl@2015.06.23%gcc@7.4.0+blas arch=linux-rhel7-power9le
module load coinhsl-2015.06.23-gcc-7.4.0-ts5vjfq
# cub@1.12.0-rc0%gcc@7.4.0 arch=linux-rhel7-power9le
module load cub-1.12.0-rc0-gcc-7.4.0-iwyj63t
# hdf5@1.12.1%gcc@7.4.0~cxx~fortran~hl~ipo~java+mpi+shared~szip~threadsafe+tools api=default build_type=RelWithDebInfo patches=ee351eb arch=linux-rhel7-power9le
module load hdf5-1.12.1-gcc-7.4.0-iofagmz
# hiop@0.5.4%gcc@7.4.0+cuda~cusolver+deepchecking~ginkgo~ipo~jsrun~kron+mpi+raja~rocm~shared+sparse build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel7-power9le
module load hiop-0.5.4-gcc-7.4.0-i3py5on
# hypre@2.24.0%gcc@7.4.0~complex~cuda~debug+fortran~gptune~int64~internal-superlu~mixedint+mpi~openmp~rocm+shared~superlu-dist~unified-memory arch=linux-rhel7-power9le
module load hypre-2.24.0-gcc-7.4.0-accevxz
# ipopt@3.12.10%gcc@7.4.0+coinhsl+debug~metis~mumps arch=linux-rhel7-power9le
module load ipopt-3.12.10-gcc-7.4.0-srvg5ag
# magma@2.6.2%gcc@7.4.0+cuda+fortran~ipo~rocm+shared build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel7-power9le
module load magma-2.6.2-gcc-7.4.0-cqkns2l
# metis@5.1.0%gcc@7.4.0~gdb~int64~real64+shared build_type=Release patches=4991da9,b1225da arch=linux-rhel7-power9le
module load metis-5.1.0-gcc-7.4.0-shhhyku
# ncurses@6.2%gcc@7.4.0~symlinks+termlib abi=none arch=linux-rhel7-power9le
module load ncurses-6.2-gcc-7.4.0-kqhmmpv
# openblas@0.3.20%gcc@7.4.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared symbol_suffix=none threads=none arch=linux-rhel7-power9le
module load openblas-0.3.20-gcc-7.4.0-3zdqw2i
# openssl@1.1.1m%gcc@7.4.0~docs certs=system arch=linux-rhel7-power9le
module load openssl-1.1.1m-gcc-7.4.0-hvsl2y7
# parmetis@4.0.3%gcc@7.4.0~gdb~int64~ipo+shared build_type=RelWithDebInfo patches=4f89253,50ed208,704b84f arch=linux-rhel7-power9le
module load parmetis-4.0.3-gcc-7.4.0-hwey2q7
# petsc@3.16.0%gcc@7.4.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw~giflib+hdf5~hpddm~hwloc+hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind+metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse+superlu-dist~tetgen~trilinos~valgrind clanguage=C arch=linux-rhel7-power9le
module load petsc-3.16.0-gcc-7.4.0-4d36twq
# pkgconf@1.8.0%gcc@7.4.0 arch=linux-rhel7-power9le
module load pkgconf-1.8.0-gcc-7.4.0-jfmmybn
# py-attrs@21.4.0%gcc@7.4.0 arch=linux-rhel7-power9le
module load py-attrs-21.4.0-gcc-7.4.0-7op4yw4
# py-iniconfig@1.1.1%gcc@7.4.0 arch=linux-rhel7-power9le
module load py-iniconfig-1.1.1-gcc-7.4.0-r7rb7ix
# py-mpi4py@3.1.2%gcc@7.4.0 arch=linux-rhel7-power9le
module load py-mpi4py-3.1.2-gcc-7.4.0-b3qgz5v
# py-packaging@21.3%gcc@7.4.0 arch=linux-rhel7-power9le
module load py-packaging-21.3-gcc-7.4.0-v25zzrs
# py-pluggy@1.0.0%gcc@7.4.0 arch=linux-rhel7-power9le
module load py-pluggy-1.0.0-gcc-7.4.0-2xlrwwx
# py-py@1.11.0%gcc@7.4.0 arch=linux-rhel7-power9le
module load py-py-1.11.0-gcc-7.4.0-bw7c7xh
# py-pyparsing@3.0.6%gcc@7.4.0 arch=linux-rhel7-power9le
module load py-pyparsing-3.0.6-gcc-7.4.0-eoshvdt
# py-pytest@6.2.5%gcc@7.4.0 arch=linux-rhel7-power9le
module load py-pytest-6.2.5-gcc-7.4.0-kmtwdv2
# py-setuptools@59.4.0%gcc@7.4.0 arch=linux-rhel7-power9le
module load py-setuptools-59.4.0-gcc-7.4.0-utnhqzc
# py-toml@0.10.2%gcc@7.4.0 arch=linux-rhel7-power9le
module load py-toml-0.10.2-gcc-7.4.0-ovce6fl
# raja@0.14.0%gcc@7.4.0+cuda~examples~exercises~ipo+openmp~rocm~shared~tests build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel7-power9le
module load raja-0.14.0-gcc-7.4.0-oxqsaxf
# superlu-dist@7.2.0%gcc@7.4.0~cuda~int64~ipo~openmp~rocm+shared build_type=RelWithDebInfo patches=8da9e21 arch=linux-rhel7-power9le
module load superlu-dist-7.2.0-gcc-7.4.0-mg2t5sg
# umpire@6.0.0%gcc@7.4.0+c+cuda~deviceconst~examples~fortran~ipo~numa~openmp~rocm~shared build_type=RelWithDebInfo cuda_arch=70 tests=none arch=linux-rhel7-power9le
module load umpire-6.0.0-gcc-7.4.0-ncjv7xt
# zlib@1.2.11%gcc@7.4.0+optimize+pic+shared arch=linux-rhel7-power9le
module load zlib-1.2.11-gcc-7.4.0-vnk3szs

# Load system modules
module load gcc/7.4.0
module load openmpi/3.1.5
module load cuda/10.2
module load python/miniconda3.8

source /share/apps/python/miniconda3.8/etc/profile.d/conda.sh

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
export EXTRA_CMAKE_ARGS="$EXTRA_CMAKE_ARGS -DCMAKE_CUDA_ARCHITECTURES=70 -DEXAGO_ENABLE_IPOPT=ON"
