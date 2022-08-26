source /etc/profile.d/modules.sh

export OMP_CANCELLATION=true
export OMP_PROC_BIND=true
export OMPI_MCA_pml="ucx"
export UCX_NET_DEVICES=mlx5_1:1,mlx5_3:1
export MY_CLUSTER=newell

module purge
module use -a /qfs/projects/exasgd/src/jaelyn-spack/spack/share/spack/modules/linux-centos8-power9le/

# Load spack modules
# berkeley-db@18.1.40%gcc@8.5.0+cxx~docs+stl patches=b231fcc arch=linux-centos8-power9le
module load berkeley-db-18.1.40-gcc-8.5.0-cuzn6qn
# blt@0.4.1%gcc@8.5.0 arch=linux-centos8-power9le
module load blt-0.4.1-gcc-8.5.0-likpa4a
# bzip2@1.0.8%gcc@8.5.0~debug~pic+shared arch=linux-centos8-power9le
module load bzip2-1.0.8-gcc-8.5.0-tsweuon
# ca-certificates-mozilla@2022-03-29%gcc@8.5.0 arch=linux-centos8-power9le
module load ca-certificates-mozilla-2022-03-29-gcc-8.5.0-zyzfhdf
# camp@0.2.2%gcc@8.5.0+cuda~ipo~rocm~tests build_type=RelWithDebInfo cuda_arch=70 arch=linux-centos8-power9le
module load camp-0.2.2-gcc-8.5.0-5po5zoy
# cmake@3.23.2%gcc@8.5.0~doc+ncurses+ownlibs~qt build_type=Release arch=linux-centos8-power9le
module load cmake-3.23.2-gcc-8.5.0-pr3l2mn
# coinhsl@2015.06.23%gcc@8.5.0+blas arch=linux-centos8-power9le
module load coinhsl-2015.06.23-gcc-8.5.0-f6ka4rc
# cub@1.16.0%gcc@8.5.0 arch=linux-centos8-power9le
module load cub-1.16.0-gcc-8.5.0-p3cnthb
# diffutils@3.8%gcc@8.5.0 arch=linux-centos8-power9le
module load diffutils-3.8-gcc-8.5.0-ppyuisg
# gdbm@1.19%gcc@8.5.0 arch=linux-centos8-power9le
module load gdbm-1.19-gcc-8.5.0-unfo3x4
# ginkgo@glu_experimental%gcc@8.5.0+cuda~develtools~full_optimizations~hwloc~ipo~oneapi+openmp~rocm+shared build_type=Release cuda_arch=70 arch=linux-centos8-power9le
module load ginkgo-glu_experimental-gcc-8.5.0-tq3ravg
# gnuconfig@2021-08-14%gcc@8.5.0 arch=linux-centos8-power9le
module load gnuconfig-2021-08-14-gcc-8.5.0-qjyg7ls
# hdf5@1.12.2%gcc@8.5.0~cxx~fortran~hl~ipo~java+mpi+shared~szip~threadsafe+tools api=default build_type=RelWithDebInfo arch=linux-centos8-power9le
module load hdf5-1.12.2-gcc-8.5.0-7dqn43r
# hiop@develop%gcc@8.5.0+cuda~cusolver+deepchecking+ginkgo~ipo~jsrun~kron+mpi+raja~rocm~shared+sparse build_type=RelWithDebInfo cuda_arch=70 dev_path=/qfs/projects/exasgd/src/jaelyn-spack/dev-8newell/hiop arch=linux-centos8-power9le
module load hiop-develop-gcc-8.5.0-22buf3r
# hypre@2.25.0%gcc@8.5.0~complex~cuda~debug+fortran~gptune~int64~internal-superlu~mixedint+mpi~openmp~rocm+shared~superlu-dist~unified-memory arch=linux-centos8-power9le
module load hypre-2.25.0-gcc-8.5.0-d72airb
# ipopt@3.12.10%gcc@8.5.0+coinhsl+debug~metis~mumps arch=linux-centos8-power9le
module load ipopt-3.12.10-gcc-8.5.0-ajkmiw6
# libiconv@1.16%gcc@8.5.0 libs=shared,static arch=linux-centos8-power9le
module load libiconv-1.16-gcc-8.5.0-qqwmnok
# magma@2.6.2%gcc@8.5.0+cuda+fortran~ipo~rocm+shared build_type=RelWithDebInfo cuda_arch=70 arch=linux-centos8-power9le
module load magma-2.6.2-gcc-8.5.0-ee3572c
# metis@5.1.0%gcc@8.5.0~gdb~int64~real64+shared build_type=Release patches=4991da9,b1225da arch=linux-centos8-power9le
module load metis-5.1.0-gcc-8.5.0-ldsei63
# ncurses@6.2%gcc@8.5.0~symlinks+termlib abi=none arch=linux-centos8-power9le
module load ncurses-6.2-gcc-8.5.0-v24hmxo
# openblas@0.3.20%gcc@8.5.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared symbol_suffix=none threads=none arch=linux-centos8-power9le
module load openblas-0.3.20-gcc-8.5.0-rwstn2s
# openssl@1.1.1q%gcc@8.5.0~docs~shared certs=mozilla patches=3fdcf2d arch=linux-centos8-power9le
module load openssl-1.1.1q-gcc-8.5.0-xlfn3bw
# parmetis@4.0.3%gcc@8.5.0~gdb~int64~ipo+shared build_type=RelWithDebInfo patches=4f89253,50ed208,704b84f arch=linux-centos8-power9le
module load parmetis-4.0.3-gcc-8.5.0-67kg5ac
# perl@5.34.1%gcc@8.5.0+cpanm+shared+threads arch=linux-centos8-power9le
module load perl-5.34.1-gcc-8.5.0-fn534xj
# petsc@3.16.6%gcc@8.5.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib+hdf5~hpddm~hwloc+hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind+metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse+superlu-dist~tetgen~trilinos~valgrind clanguage=C arch=linux-centos8-power9le
module load petsc-3.16.6-gcc-8.5.0-vqix3dp
# pkgconf@1.8.0%gcc@8.5.0 arch=linux-centos8-power9le
module load pkgconf-1.8.0-gcc-8.5.0-imrnro2
# py-attrs@21.4.0%gcc@8.5.0 arch=linux-centos8-power9le
module load py-attrs-21.4.0-gcc-8.5.0-oejgn27
# py-flit-core@3.6.0%gcc@8.5.0 arch=linux-centos8-power9le
module load py-flit-core-3.6.0-gcc-8.5.0-3xhtf2d
# py-iniconfig@1.1.1%gcc@8.5.0 arch=linux-centos8-power9le
module load py-iniconfig-1.1.1-gcc-8.5.0-lcqxzow
# py-mpi4py@3.1.2%gcc@8.5.0 arch=linux-centos8-power9le
module load py-mpi4py-3.1.2-gcc-8.5.0-fcmbnkc
# py-packaging@21.3%gcc@8.5.0 arch=linux-centos8-power9le
module load py-packaging-21.3-gcc-8.5.0-pjxzk4f
# py-pip@21.3.1%gcc@8.5.0 arch=linux-centos8-power9le
module load py-pip-21.3.1-gcc-8.5.0-2edb3hu
# py-pluggy@1.0.0%gcc@8.5.0 arch=linux-centos8-power9le
module load py-pluggy-1.0.0-gcc-8.5.0-z35kd4h
# py-py@1.11.0%gcc@8.5.0 arch=linux-centos8-power9le
module load py-py-1.11.0-gcc-8.5.0-vo6wsnv
# py-pyparsing@3.0.6%gcc@8.5.0 arch=linux-centos8-power9le
module load py-pyparsing-3.0.6-gcc-8.5.0-6lqucvi
# py-pytest@6.2.5%gcc@8.5.0 arch=linux-centos8-power9le
module load py-pytest-6.2.5-gcc-8.5.0-b3pcszv
# py-setuptools@63.0.0%gcc@8.5.0 arch=linux-centos8-power9le
module load py-setuptools-63.0.0-gcc-8.5.0-iyruu6l
# py-setuptools-scm@7.0.3%gcc@8.5.0+toml arch=linux-centos8-power9le
module load py-setuptools-scm-7.0.3-gcc-8.5.0-pecvfim
# py-toml@0.10.2%gcc@8.5.0 arch=linux-centos8-power9le
module load py-toml-0.10.2-gcc-8.5.0-ptjfjva
# py-tomli@1.2.2%gcc@8.5.0 arch=linux-centos8-power9le
module load py-tomli-1.2.2-gcc-8.5.0-p7c7ssb
# py-typing-extensions@4.3.0%gcc@8.5.0 arch=linux-centos8-power9le
module load py-typing-extensions-4.3.0-gcc-8.5.0-odxczqu
# py-wheel@0.37.0%gcc@8.5.0 arch=linux-centos8-power9le
module load py-wheel-0.37.0-gcc-8.5.0-cfcufxy
# raja@0.14.0%gcc@8.5.0+cuda~examples~exercises~ipo+openmp~rocm~shared~tests build_type=RelWithDebInfo cuda_arch=70 arch=linux-centos8-power9le
module load raja-0.14.0-gcc-8.5.0-uq3jvzh
# readline@8.1.2%gcc@8.5.0 arch=linux-centos8-power9le
module load readline-8.1.2-gcc-8.5.0-l4hzlyf
# superlu-dist@7.2.0%gcc@8.5.0~cuda~int64~ipo~openmp~rocm+shared build_type=RelWithDebInfo patches=8da9e21 arch=linux-centos8-power9le
module load superlu-dist-7.2.0-gcc-8.5.0-xtm2lfm
# umpire@6.0.0%gcc@8.5.0+c+cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp~rocm~shared build_type=RelWithDebInfo cuda_arch=70 tests=none arch=linux-centos8-power9le
module load umpire-6.0.0-gcc-8.5.0-43tcar4
# zlib@1.2.12%gcc@8.5.0+optimize+pic+shared patches=0d38234 arch=linux-centos8-power9le
module load zlib-1.2.12-gcc-8.5.0-spb5k73

# Load system modules
module load gcc/8.5.0
module load openmpi/4.1.4
module load cuda/11.4
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
