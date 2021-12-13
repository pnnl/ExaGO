export MY_CLUSTER=summit
export PROJ_DIR=/autofs/nccs-svm1_proj/csc359

module purge

module use -a /gpfs/alpine/csc359/proj-shared/asher/spack/share/spack/modules/linux-rhel8-power9le/

# ipopt@3.12.10%gcc@10.2.0+coinhsl~debug+metis~mumps arch=linux-rhel8-power9le
module load ipopt-3.12.10-gcc-10.2.0-shqkmmy
# coinhsl@2019.05.21%gcc@10.2.0+blas arch=linux-rhel8-power9le
module load coinhsl-2019.05.21-gcc-10.2.0-oxdsnrx
# berkeley-db@18.1.40%gcc@10.2.0+cxx~docs+stl patches=b231fcc4d5cff05e5c3a4814f6a5af0e9a966428dc2176540d2c05aff41de522 arch=linux-rhel8-power9le
module load berkeley-db-18.1.40-gcc-10.2.0-6bnatyc
# blt@0.4.1%gcc@10.2.0 arch=linux-rhel8-power9le
module load blt-0.4.1-gcc-10.2.0-y5x27ua
# bzip2@1.0.8%gcc@10.2.0~debug~pic+shared arch=linux-rhel8-power9le
module load bzip2-1.0.8-gcc-10.2.0-pdiku55
# camp@0.2.2%gcc@10.2.0+cuda~ipo~rocm~tests build_type=RelWithDebInfo cuda_arch=none arch=linux-rhel8-power9le
module load camp-0.2.2-gcc-10.2.0-hhtibao
# cmake@3.21.3%gcc@10.2.0~doc+ncurses+openssl+ownlibs~qt build_type=Release arch=linux-rhel8-power9le
module load cmake-3.21.3-gcc-10.2.0-affsbrr
# cub@1.12.0-rc0%gcc@10.2.0 arch=linux-rhel8-power9le
module load cub-1.12.0-rc0-gcc-10.2.0-nayq7qd
# cuda@11.1.1%gcc@10.2.0~dev arch=linux-rhel8-power9le
module load cuda-11.1.1-gcc-10.2.0-itd73gl
# diffutils@3.8%gcc@10.2.0 arch=linux-rhel8-power9le
module load diffutils-3.8-gcc-10.2.0-n2gttub
# gdbm@1.19%gcc@10.2.0 arch=linux-rhel8-power9le
module load gdbm-1.19-gcc-10.2.0-ab4rekw
# gnuconfig@2021-08-14%gcc@10.2.0 arch=linux-rhel8-power9le
module load gnuconfig-2021-08-14-gcc-10.2.0-j2pelq3
# hiop@0.5.3%gcc@10.2.0+cuda~deepchecking~ipo~jsrun~kron+mpi+raja~rocm~shared~sparse build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel8-power9le
module load hiop-0.5.3-gcc-10.2.0-pjhusig
# libiconv@1.16%gcc@10.2.0 libs=shared,static arch=linux-rhel8-power9le
module load libiconv-1.16-gcc-10.2.0-cpn6euq
# magma@2.6.1%gcc@10.2.0+cuda+fortran~ipo~rocm+shared build_type=RelWithDebInfo cuda_arch=none arch=linux-rhel8-power9le
module load magma-2.6.1-gcc-10.2.0-swf7uk7
# ncurses@6.2%gcc@10.2.0~symlinks+termlib abi=none arch=linux-rhel8-power9le
module load ncurses-6.2-gcc-10.2.0-y222xlf
# openblas@0.3.18%gcc@10.2.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared symbol_suffix=none threads=none arch=linux-rhel8-power9le
module load openblas-0.3.18-gcc-10.2.0-r4dabzt
# perl@5.34.0%gcc@10.2.0+cpanm+shared+threads arch=linux-rhel8-power9le
module load perl-5.34.0-gcc-10.2.0-q7baup2
# petsc@3.14.6%gcc@10.2.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw~giflib+hdf5~hpddm~hwloc+hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind+metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse+superlu-dist~tetgen~trilinos~valgrind clanguage=C arch=linux-rhel8-power9le
module load petsc-3.14.6-gcc-10.2.0-jwyafkj
# pkgconf@1.8.0%gcc@10.2.0 arch=linux-rhel8-power9le
module load pkgconf-1.8.0-gcc-10.2.0-nngfigd
# raja@0.14.0%gcc@10.2.0+cuda+examples+exercises~ipo+openmp~rocm+shared~tests build_type=RelWithDebInfo cuda_arch=none arch=linux-rhel8-power9le
module load raja-0.14.0-gcc-10.2.0-l2i7ckj
# readline@8.1%gcc@10.2.0 arch=linux-rhel8-power9le
module load readline-8.1-gcc-10.2.0-dr6yllv
# spectrum-mpi@10.4.0.3-20210112%gcc@10.2.0 arch=linux-rhel8-power9le
module load spectrum-mpi-10.4.0.3-20210112-gcc-10.2.0-nxf5xco
# umpire@6.0.0%gcc@10.2.0+c+cuda~deviceconst+examples~fortran~ipo~numa~openmp~rocm~shared build_type=RelWithDebInfo cuda_arch=none tests=none arch=linux-rhel8-power9le
module load umpire-6.0.0-gcc-10.2.0-t4xicvq
# zlib@1.2.11%gcc@10.2.0+optimize+pic+shared arch=linux-rhel8-power9le
module load zlib-1.2.11-gcc-10.2.0-va7dndx

module load gcc/10.2.0
module load cuda/11.1.1
module load spectrum-mpi/10.4.0.3-20210112
module load magma/2.6.1
module load petsc/3.14.6
module load raja/0.14.0
module load cmake/3.21.3
module load metis/5.1.0

export CC=/sw/summit/gcc/10.2.0-2/bin/gcc
export CXX=/sw/summit/gcc/10.2.0-2/bin/g++
export FC=/sw/summit/gcc/10.2.0-2/bin/gfotran

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
