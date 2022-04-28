source /etc/profile.d/modules.sh

export OMP_CANCELLATION=true
export OMP_PROC_BIND=true
export OMPI_MCA_pml="ucx"
export UCX_NET_DEVICES=mlx5_1:1,mlx5_3:1
export MY_CLUSTER=newell

module purge
module use -a /qfs/projects/exasgd/src/cameron-spack/share/spack/modules/linux-rhel7-power9le/

# Load spack modules
# autoconf@2.69%gcc@7.4.0 patches=35c4492,7793209,a49dd5b arch=linux-rhel7-power9le
module load autoconf-2.69-gcc-7.4.0-jcxtgdk
# autoconf-archive@2019.01.06%gcc@7.4.0 arch=linux-rhel7-power9le
module load autoconf-archive-2019.01.06-gcc-7.4.0-nn453cx
# automake@1.16.5%gcc@7.4.0 arch=linux-rhel7-power9le
module load automake-1.16.5-gcc-7.4.0-av2w6do
# berkeley-db@18.1.40%gcc@7.4.0+cxx~docs+stl patches=b231fcc arch=linux-rhel7-power9le
module load berkeley-db-18.1.40-gcc-7.4.0-ic4cqif
# blt@0.4.1%gcc@7.4.0 arch=linux-rhel7-power9le
module load blt-0.4.1-gcc-7.4.0-2th7jgq
# bzip2@1.0.8%gcc@7.4.0~debug~pic+shared arch=linux-rhel7-power9le
module load bzip2-1.0.8-gcc-7.4.0-jty62q7
# camp@0.2.2%gcc@7.4.0+cuda~ipo~rocm~tests build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel7-power9le
module load camp-0.2.2-gcc-7.4.0-vsu2jwh
# cmake@3.23.1%gcc@7.4.0~doc+ncurses+ownlibs~qt build_type=Release arch=linux-rhel7-power9le
module load cmake-3.23.1-gcc-7.4.0-ckfugtf
# coinhsl@2015.06.23%gcc@7.4.0+blas arch=linux-rhel7-power9le
module load coinhsl-2015.06.23-gcc-7.4.0-ts5vjfq
# cub@1.12.0%gcc@7.4.0 arch=linux-rhel7-power9le
module load cub-1.12.0-gcc-7.4.0-4qyvoqn
# diffutils@3.8%gcc@7.4.0 arch=linux-rhel7-power9le
module load diffutils-3.8-gcc-7.4.0-cy55hsj
# gdbm@1.19%gcc@7.4.0 arch=linux-rhel7-power9le
module load gdbm-1.19-gcc-7.4.0-ahdwucz
# gmp@6.2.1%gcc@7.4.0 libs=shared,static arch=linux-rhel7-power9le
module load gmp-6.2.1-gcc-7.4.0-oea2aet
# gnuconfig@2021-08-14%gcc@7.4.0 arch=linux-rhel7-power9le
module load gnuconfig-2021-08-14-gcc-7.4.0-qr6nxuq
# hdf5@1.12.2%gcc@7.4.0~cxx~fortran~hl~ipo~java+mpi+shared~szip~threadsafe+tools api=default build_type=RelWithDebInfo arch=linux-rhel7-power9le
module load hdf5-1.12.2-gcc-7.4.0-qtzvwle
# hiop@0.6.1%gcc@7.4.0+cuda+cusolver+deepchecking~ginkgo~ipo~jsrun+kron+mpi+raja~rocm~shared+sparse build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel7-power9le
module load hiop-0.6.1-gcc-7.4.0-75ogmba
# hypre@2.24.0%gcc@7.4.0~complex~cuda~debug+fortran~gptune~int64~internal-superlu~mixedint+mpi~openmp~rocm+shared~superlu-dist~unified-memory arch=linux-rhel7-power9le
module load hypre-2.24.0-gcc-7.4.0-2rkp5iw
# ipopt@3.12.10%gcc@7.4.0+coinhsl+debug~metis~mumps arch=linux-rhel7-power9le
module load ipopt-3.12.10-gcc-7.4.0-srvg5ag
# libiconv@1.16%gcc@7.4.0 libs=shared,static arch=linux-rhel7-power9le
module load libiconv-1.16-gcc-7.4.0-idqno7d
# libsigsegv@2.13%gcc@7.4.0 arch=linux-rhel7-power9le
module load libsigsegv-2.13-gcc-7.4.0-cbn4dja
# libtool@2.4.7%gcc@7.4.0 arch=linux-rhel7-power9le
module load libtool-2.4.7-gcc-7.4.0-p5juddc
# m4@1.4.19%gcc@7.4.0+sigsegv patches=9dc5fbd,bfdffa7 arch=linux-rhel7-power9le
module load m4-1.4.19-gcc-7.4.0-nrrlksm
# magma@2.6.2%gcc@7.4.0+cuda+fortran~ipo~rocm+shared build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel7-power9le
module load magma-2.6.2-gcc-7.4.0-6yuqfpm
# metis@5.1.0%gcc@7.4.0~gdb~int64~real64+shared build_type=Release patches=4991da9,b1225da arch=linux-rhel7-power9le
module load metis-5.1.0-gcc-7.4.0-shhhyku
# mpfr@4.1.0%gcc@7.4.0 libs=shared,static arch=linux-rhel7-power9le
module load mpfr-4.1.0-gcc-7.4.0-tz5esun
# ncurses@6.2%gcc@7.4.0~symlinks+termlib abi=none arch=linux-rhel7-power9le
module load ncurses-6.2-gcc-7.4.0-kqhmmpv
# openblas@0.3.20%gcc@7.4.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared symbol_suffix=none threads=none arch=linux-rhel7-power9le
module load openblas-0.3.20-gcc-7.4.0-3zdqw2i
# openssl@1.1.1n%gcc@7.4.0~docs~shared certs=system arch=linux-rhel7-power9le
module load openssl-1.1.1n-gcc-7.4.0-mxyibws
# parmetis@4.0.3%gcc@7.4.0~gdb~int64~ipo+shared build_type=RelWithDebInfo patches=4f89253,50ed208,704b84f arch=linux-rhel7-power9le
module load parmetis-4.0.3-gcc-7.4.0-c7vxddy
# perl@5.34.1%gcc@7.4.0+cpanm+shared+threads arch=linux-rhel7-power9le
module load perl-5.34.1-gcc-7.4.0-zy3h5bd
# petsc@3.16.6%gcc@7.4.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib+hdf5~hpddm~hwloc+hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind+metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse+superlu-dist~tetgen~trilinos~valgrind clanguage=C arch=linux-rhel7-power9le
module load petsc-3.16.6-gcc-7.4.0-rs3dccz
# pkgconf@1.8.0%gcc@7.4.0 arch=linux-rhel7-power9le
module load pkgconf-1.8.0-gcc-7.4.0-jfmmybn
# py-mpi4py@3.1.2%gcc@7.4.0 arch=linux-rhel7-power9le
module load py-mpi4py-3.1.2-gcc-7.4.0-nozxuwj
# py-pip@21.3.1%gcc@7.4.0 arch=linux-rhel7-power9le
module load py-pip-21.3.1-gcc-7.4.0-s7m26a2
# py-setuptools@59.4.0%gcc@7.4.0 arch=linux-rhel7-power9le
module load py-setuptools-59.4.0-gcc-7.4.0-utnhqzc
# py-wheel@0.37.0%gcc@7.4.0 arch=linux-rhel7-power9le
module load py-wheel-0.37.0-gcc-7.4.0-zsolxgr
# python@3.8.5%gcc@7.4.0+bz2+ctypes+dbm~debug+ensurepip+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tix~tkinter~ucs4+uuid+zlib patches=0d98e93,4c24573,f2fd060 arch=linux-rhel7-power9le
module load python-3.8.5-gcc-7.4.0-szartfv
# raja@0.14.0%gcc@7.4.0+cuda~examples~exercises~ipo+openmp~rocm+shared~tests build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel7-power9le
module load raja-0.14.0-gcc-7.4.0-sew5thv
# readline@8.1%gcc@7.4.0 arch=linux-rhel7-power9le
module load readline-8.1-gcc-7.4.0-cszha3o
# suite-sparse@5.10.1%gcc@7.4.0~cuda~graphblas~openmp+pic~tbb arch=linux-rhel7-power9le
module load suite-sparse-5.10.1-gcc-7.4.0-e5qockg
# superlu-dist@7.2.0%gcc@7.4.0~cuda~int64~ipo~openmp~rocm+shared build_type=RelWithDebInfo patches=8da9e21 arch=linux-rhel7-power9le
module load superlu-dist-7.2.0-gcc-7.4.0-34pmq2m
# texinfo@6.5%gcc@7.4.0 patches=12f6edb,1732115 arch=linux-rhel7-power9le
module load texinfo-6.5-gcc-7.4.0-t5dt7rq
# umpire@6.0.0%gcc@7.4.0+c+cuda~deviceconst~examples~fortran~ipo~numa~openmp~rocm~shared build_type=RelWithDebInfo cuda_arch=70 tests=none arch=linux-rhel7-power9le
module load umpire-6.0.0-gcc-7.4.0-rpwrj4p
# zlib@1.2.12%gcc@7.4.0+optimize+pic+shared patches=0d38234 arch=linux-rhel7-power9le
module load zlib-1.2.12-gcc-7.4.0-d6xlzc6

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
