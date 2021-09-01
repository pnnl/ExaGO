export MY_CLUSTER=summit
export PROJ_DIR=/autofs/nccs-svm1_proj/csc359

module use -a /autofs/nccs-svm1_proj/csc359/installs/spack/share/spack/modules/linux-rhel8-power9le/

module purge

# autoconf-archive@2019.01.06%gcc@9.1.0 arch=linux-rhel8-power9le
module load exasgd-autoconf-archive/2019.01.06/gcc-9.1.0-ofpeiyd
# berkeley-db@18.1.40%gcc@9.1.0+cxx~docs+stl patches=b231fcc4d5cff05e5c3a4814f6a5af0e9a966428dc2176540d2c05aff41de522 arch=linux-rhel8-power9le
module load exasgd-berkeley-db/18.1.40/gcc-9.1.0-qwxozkk
# bzip2@1.0.8%gcc@9.1.0~debug~pic+shared arch=linux-rhel8-power9le
module load exasgd-bzip2/1.0.8/gcc-9.1.0-sdo4em2
# coinhsl@2015.06.23%gcc@9.1.0+blas arch=linux-rhel8-power9le
module load exasgd-coinhsl/2015.06.23/gcc-9.1.0-sgfgayg
# cuda@11.0.3%gcc@9.1.0~allow-unsupported-compilers~dev arch=linux-rhel8-power9le
module load exasgd-cuda/11.0.3/gcc-9.1.0-bd5ppam
# exago@develop%gcc@9.1.0+cuda+gpu+hiop~hip~ipo+ipopt+mpi+petsc+raja amdgpu_target=none build_type=Release cuda_arch=70 arch=linux-rhel8-power9le
module load exasgd-exago/develop/cuda-11.0.3/spectrum-mpi-10.4.0.3-20210112/gcc-9.1.0-3j3p46s
# gdbm@1.19%gcc@9.1.0 arch=linux-rhel8-power9le
module load exasgd-gdbm/1.19/gcc-9.1.0-g4xlywa
# gmp@6.2.1%gcc@9.1.0 arch=linux-rhel8-power9le
module load exasgd-gmp/6.2.1/gcc-9.1.0-5ogv5uk
# hiop@0.4.1%gcc@9.1.0+cuda+deepchecking+gpu~hip~ipo~jsrun+kron+mpi+raja+shared~sparse~srun amdgpu_target=none build_type=Release cuda_arch=70 arch=linux-rhel8-power9le
module load exasgd-hiop/0.4.1/cuda-11.0.3/spectrum-mpi-10.4.0.3-20210112/gcc-9.1.0-5veh4g6
# ipopt@3.12.10%gcc@9.1.0+coinhsl~debug~mumps arch=linux-rhel8-power9le
module load exasgd-ipopt/3.12.10/gcc-9.1.0-l2ylu4z
# libsigsegv@2.13%gcc@9.1.0 arch=linux-rhel8-power9le
module load exasgd-libsigsegv/2.13/gcc-9.1.0-mxofdmp
# m4@1.4.19%gcc@9.1.0+sigsegv arch=linux-rhel8-power9le
module load exasgd-m4/1.4.19/gcc-9.1.0-onfjkyf
# magma@2.5.4%gcc@9.1.0+cuda+fortran~ipo~rocm+shared amdgpu_target=none build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel8-power9le
module load exasgd-magma/2.5.4/cuda-11.0.3/gcc-9.1.0-wxkchiu
# metis@5.1.0%gcc@9.1.0~gdb~int64~real64+shared build_type=Release patches=4991da938c1d3a1d3dea78e49bbebecba00273f98df2a656e38b83d55b281da1,b1225da886605ea558db7ac08dd8054742ea5afe5ed61ad4d0fe7a495b1270d2 arch=linux-rhel8-power9le
module load exasgd-metis/5.1.0/gcc-9.1.0-yueb4wp
# mpfr@4.1.0%gcc@9.1.0 arch=linux-rhel8-power9le
module load exasgd-mpfr/4.1.0/gcc-9.1.0-azpm7od
# openblas@0.3.17%gcc@9.1.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared threads=none arch=linux-rhel8-power9le
module load exasgd-openblas/0.3.17/gcc-9.1.0-q3wsicq
# perl@5.34.0%gcc@9.1.0+cpanm+shared+threads arch=linux-rhel8-power9le
module load exasgd-perl/5.34.0/gcc-9.1.0-sdm3z2p
# petsc@3.14.6%gcc@9.1.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw~giflib~hdf5~hip~hypre~int64~jpeg~knl~libpng~libyaml~memkind~metis~mkl-pardiso~moab~mpfr+mpi~mumps~p4est~ptscotch~random123~saws+shared~suite-sparse~superlu-dist~trilinos~valgrind amdgpu_target=none clanguage=C arch=linux-rhel8-power9le
module load exasgd-petsc/3.14.6/spectrum-mpi-10.4.0.3-20210112/gcc-9.1.0-ysnm7nb
# raja@0.12.1%gcc@9.1.0+cuda~examples~exercises~hip~ipo+openmp+shared~tests amdgpu_target=none build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel8-power9le
module load exasgd-raja/0.12.1/cuda-11.0.3/gcc-9.1.0-32tj3no
# readline@8.0%gcc@9.1.0 arch=linux-rhel8-power9le
module load exasgd-readline/8.0/gcc-9.1.0-75jtlrm
# spectrum-mpi@10.4.0.3-20210112%gcc@9.1.0 arch=linux-rhel8-power9le
module load exasgd-spectrum-mpi/10.4.0.3-20210112/gcc-9.1.0-pmqu4l5
# suite-sparse@5.10.1%gcc@9.1.0~cuda~openmp+pic~tbb arch=linux-rhel8-power9le
module load exasgd-suite-sparse/5.10.1/gcc-9.1.0-jqwvdui
# texinfo@6.5%gcc@9.1.0 patches=12f6edb0c6b270b8c8dba2ce17998c580db01182d871ee32b7b6e4129bd1d23a,1732115f651cff98989cb0215d8f64da5e0f7911ebf0c13b064920f088f2ffe1 arch=linux-rhel8-power9le
module load exasgd-texinfo/6.5/gcc-9.1.0-g7fbyu3
# umpire@4.1.2%gcc@9.1.0+c+cuda~deviceconst~examples~fortran~hip~ipo~numa+openmp~shared amdgpu_target=none build_type=RelWithDebInfo cuda_arch=none patches=7d912d31cd293df005ba74cb96c6f3e32dc3d84afff49b14509714283693db08 tests=none arch=linux-rhel8-power9le
module load exasgd-umpire/4.1.2/cuda-11.0.3/gcc-9.1.0-qhsp5nb
# zlib@1.2.11%gcc@9.1.0+optimize+pic+shared arch=linux-rhel8-power9le
module load exasgd-zlib/1.2.11/gcc-9.1.0-utn462j

# System modules
module load DefApps
module load gcc/9.1.0
module load cuda/11.0.3
module load spectrum-mpi/10.3.1.2-20200121
module load cmake

export CC=/sw/summit/gcc/9.3.0-2/bin/gcc CXX=/sw/summit/gcc/9.3.0-2/bin/g++

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
