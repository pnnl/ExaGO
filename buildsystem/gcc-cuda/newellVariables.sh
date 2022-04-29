source /etc/profile.d/modules.sh

export OMP_CANCELLATION=true
export OMP_PROC_BIND=true
export OMPI_MCA_pml="ucx"
export UCX_NET_DEVICES=mlx5_1:1,mlx5_3:1
export MY_CLUSTER=newell

module purge
module use -a /qfs/projects/exasgd/src/cameron-spack/share/spack/modules/linux-centos7-haswell/
module use -a /qfs/projects/exasgd/src/cameron-spack/share/spack/modules/linux-centos7-broadwell/
module use -a /qfs/projects/exasgd/src/cameron-spack/share/spack/modules/linux-rhel7-power9le/
module use -a /qfs/projects/exasgd/src/spack/share/spack/modules/linux-rhel7-power9le/

# Load spack modules
# berkeley-db@18.1.40%gcc@7.4.0+cxx~docs+stl patches=b231fcc4d5cff05e5c3a4814f6a5af0e9a966428dc2176540d2c05aff41de522 arch=linux-rhel7-power9le
module load berkeley-db-18.1.40-gcc-7.4.0-ic4cqif
# blt@0.4.1%gcc@7.4.0 arch=linux-rhel7-power9le
module load blt-0.4.1-gcc-7.4.0-vbugf3i
# bzip2@1.0.8%gcc@7.4.0~debug~pic+shared arch=linux-rhel7-power9le
module load bzip2-1.0.8-gcc-7.4.0-jty62q7
# camp@0.2.2%gcc@7.4.0+cuda~ipo~rocm~tests build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel7-power9le
module load camp-0.2.2-gcc-7.4.0-tdeemi3
# cmake@3.22.2%gcc@7.4.0~doc+ncurses+openssl+ownlibs~qt build_type=Release arch=linux-rhel7-power9le
module load cmake-3.22.2-gcc-7.4.0-mt7gms5
# coinhsl@2015.06.23%gcc@7.4.0+blas arch=linux-rhel7-power9le
module load coinhsl-2015.06.23-gcc-7.4.0-udagsad
# cub@1.12.0-rc0%gcc@7.4.0 arch=linux-rhel7-power9le
module load cub-1.12.0-rc0-gcc-7.4.0-iwyj63t
# diffutils@3.8%gcc@7.4.0 arch=linux-rhel7-power9le
module load diffutils-3.8-gcc-7.4.0-cy55hsj
# gdbm@1.19%gcc@7.4.0 arch=linux-rhel7-power9le
module load gdbm-1.19-gcc-7.4.0-ahdwucz
# gnuconfig@2021-08-14%gcc@7.4.0 arch=linux-rhel7-power9le
module load gnuconfig-2021-08-14-gcc-7.4.0-qr6nxuq
# hdf5@1.12.1%gcc@7.4.0~cxx~fortran~hl~ipo~java+mpi+shared~szip~threadsafe+tools api=default build_type=RelWithDebInfo arch=linux-rhel7-power9le
module load hdf5-1.12.1-gcc-7.4.0-hjyes2v
# hiop@0.5.3%gcc@7.4.0+cuda+deepchecking~ipo~jsrun~kron+mpi+raja~rocm~shared+sparse build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel7-power9le
module load hiop-0.5.3-gcc-7.4.0-7kgfuxh
# hypre@2.23.0%gcc@7.4.0~complex~cuda~debug+fortran~gptune~int64~internal-superlu~mixedint+mpi~openmp+shared~superlu-dist~unified-memory arch=linux-rhel7-power9le
module load hypre-2.23.0-gcc-7.4.0-7d3ugmb
# ipopt@3.12.10%gcc@7.4.0+coinhsl+debug~metis~mumps arch=linux-rhel7-power9le
module load ipopt-3.12.10-gcc-7.4.0-6qv5eb5
# libiconv@1.16%gcc@7.4.0 libs=shared,static arch=linux-rhel7-power9le
module load libiconv-1.16-gcc-7.4.0-idqno7d
# magma@2.6.1%gcc@7.4.0+cuda+fortran~ipo~rocm+shared build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel7-power9le
module load magma-2.6.1-gcc-7.4.0-huyhrb5
# metis@5.1.0%gcc@7.4.0~gdb~int64~real64+shared build_type=Release patches=4991da938c1d3a1d3dea78e49bbebecba00273f98df2a656e38b83d55b281da1,b1225da886605ea558db7ac08dd8054742ea5afe5ed61ad4d0fe7a495b1270d2 arch=linux-rhel7-power9le
module load metis-5.1.0-gcc-7.4.0-shhhyku
# ncurses@6.2%gcc@7.4.0~symlinks+termlib abi=none arch=linux-rhel7-power9le
module load ncurses-6.2-gcc-7.4.0-kqhmmpv
# openblas@0.3.19%gcc@7.4.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared symbol_suffix=none threads=none arch=linux-rhel7-power9le
module load openblas-0.3.19-gcc-7.4.0-w63lkax
# openmpi@3.1.5%gcc@7.4.0~atomics~cuda~cxx~cxx_exceptions+gpfs~internal-hwloc~java~legacylaunchers~lustre~memchecker~pmi~pmix+romio~singularity~sqlite3+static~thread_multiple+vt+wrapper-rpath fabrics=none schedulers=none arch=linux-rhel7-power9le
module load openmpi-3.1.5-gcc-7.4.0-e2d65hr
# openssl@1.1.1m%gcc@7.4.0~docs certs=system arch=linux-rhel7-power9le
module load openssl-1.1.1m-gcc-7.4.0-hvsl2y7
# parmetis@4.0.3%gcc@7.4.0~gdb~int64~ipo+shared build_type=RelWithDebInfo patches=4f892531eb0a807eb1b82e683a416d3e35154a455274cf9b162fb02054d11a5b,50ed2081bc939269689789942067c58b3e522c269269a430d5d34c00edbc5870,704b84f7c7444d4372cb59cca6e1209df4ef3b033bc4ee3cf50f369bce972a9d arch=linux-rhel7-power9le
module load parmetis-4.0.3-gcc-7.4.0-ia6pql6
# perl@5.34.0%gcc@7.4.0+cpanm+shared+threads arch=linux-rhel7-power9le
module load perl-5.34.0-gcc-7.4.0-h45ivzd
# petsc@3.16.0%gcc@7.4.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw~giflib+hdf5~hpddm~hwloc+hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind+metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse+superlu-dist~tetgen~trilinos~valgrind clanguage=C arch=linux-rhel7-power9le
module load petsc-3.16.0-gcc-7.4.0-3hto34j
# pkgconf@1.8.0%gcc@7.4.0 arch=linux-rhel7-power9le
module load pkgconf-1.8.0-gcc-7.4.0-jfmmybn
# py-mpi4py@3.1.2%gcc@7.4.0 arch=linux-rhel7-power9le
module load py-mpi4py-3.1.2-gcc-7.4.0-ogt62xc
# py-pip@21.3.1%gcc@7.4.0 arch=linux-rhel7-power9le
module load py-pip-21.3.1-gcc-7.4.0-s7m26a2
# py-setuptools@59.4.0%gcc@7.4.0 arch=linux-rhel7-power9le
module load py-setuptools-59.4.0-gcc-7.4.0-utnhqzc
# py-wheel@0.37.0%gcc@7.4.0 arch=linux-rhel7-power9le
module load py-wheel-0.37.0-gcc-7.4.0-zsolxgr
# python@3.8.5%gcc@7.4.0+bz2+ctypes+dbm~debug+ensurepip+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tix~tkinter~ucs4+uuid+zlib patches=0d98e93189bc278fbc37a50ed7f183bd8aaf249a8e1670a465f0db6bb4f8cf87,4c2457325f2b608b1b6a2c63087df8c26e07db3e3d493caf36a56f0ecf6fb768,f2fd060afc4b4618fe8104c4c5d771f36dc55b1db5a4623785a4ea707ec72fb4 arch=linux-rhel7-power9le
module load python-3.8.5-gcc-7.4.0-szartfv
# raja@0.14.0%gcc@7.4.0+cuda~examples~exercises~ipo+openmp~rocm~shared~tests build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel7-power9le
module load raja-0.14.0-gcc-7.4.0-oxqsaxf
# readline@8.1%gcc@7.4.0 arch=linux-rhel7-power9le
module load readline-8.1-gcc-7.4.0-cszha3o
# superlu-dist@7.2.0%gcc@7.4.0~cuda~int64~ipo~openmp~rocm+shared build_type=RelWithDebInfo patches=8da9e21f724e8f11a5782960cc0322e12f724bf93cead7df517901a788ea3d61 arch=linux-rhel7-power9le
module load superlu-dist-7.2.0-gcc-7.4.0-hy6qjol
# umpire@6.0.0%gcc@7.4.0+c+cuda~deviceconst~examples~fortran~ipo~numa~openmp~rocm~shared build_type=RelWithDebInfo cuda_arch=70 tests=none arch=linux-rhel7-power9le
module load umpire-6.0.0-gcc-7.4.0-ncjv7xt
# zlib@1.2.11%gcc@7.4.0+optimize+pic+shared arch=linux-rhel7-power9le
module load zlib-1.2.11-gcc-7.4.0-vnk3szs
# suite-sparse@5.10.1%gcc@7.4.0~cuda~graphblas~openmp+pic~tbb arch=linux-rhel7-power9le
module load suite-sparse-5.10.1-gcc-7.4.0-btuc2bk


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
