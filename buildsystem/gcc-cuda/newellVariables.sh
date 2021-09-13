source /etc/profile.d/modules.sh
export OMP_CANCELLATION=true
export OMP_PROC_BIND=true
export OMPI_MCA_pml="ucx"
export UCX_NET_DEVICES=mlx5_1:1,mlx5_3:1
export MY_CLUSTER=newell

module purge
module use -a /qfs/projects/exasgd/src/spack/share/spack/modules/linux-rhel7-power9le/

# Load spack modules
# autoconf-archive@2019.01.06%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-autoconf-archive/2019.01.06/gcc-7.4.0-zr3h7p2
# bzip2@1.0.8%gcc@7.4.0~debug~pic+shared arch=linux-rhel7-power9le
module load exasgd-bzip2/1.0.8/gcc-7.4.0-zm2dl2f
# cmake@3.20.2%gcc@7.4.0~doc+ncurses+openssl+ownlibs~qt build_type=Release arch=linux-rhel7-power9le
module load exasgd-cmake/3.20.2/gcc-7.4.0-rfq4n6b
# coinhsl@2015.06.23%gcc@7.4.0+blas arch=linux-rhel7-power9le
module load exasgd-coinhsl/2015.06.23/gcc-7.4.0-k5ullmz
# expat@2.3.0%gcc@7.4.0+libbsd arch=linux-rhel7-power9le
module load exasgd-expat/2.3.0/gcc-7.4.0-kphxbv3
# gdbm@1.19%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-gdbm/1.19/gcc-7.4.0-v4st5ql
# gettext@0.21%gcc@7.4.0+bzip2+curses+git~libunistring+libxml2+tar+xz arch=linux-rhel7-power9le
module load exasgd-gettext/0.21/gcc-7.4.0-ig464ci
# gmp@6.2.1%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-gmp/6.2.1/gcc-7.4.0-6svtfvr
# hdf5@1.10.7%gcc@7.4.0~cxx~debug~fortran~hl~java+mpi+pic+shared~szip~threadsafe api=none arch=linux-rhel7-power9le
module load exasgd-hdf5/1.10.7/openmpi-3.1.5/gcc-7.4.0-g2s5bk6
# hiop@0.4.1%gcc@7.4.0+cuda~deepchecking+gpu~hip~ipo~jsrun+kron+mpi+raja+shared~sparse~srun amdgpu_target=none build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel7-power9le

# module load exasgd-hiop/0.4.1/cuda-10.2.89-system/openmpi-3.1.5/gcc-7.4.0-mwgv3uv

module load exasgd-hiop/0.4.6/cuda-10.2.89-system/openmpi-3.1.5/gcc-7.4.0-fbbiwc7

# hypre@2.20.0%gcc@7.4.0~complex~cuda~debug~int64~internal-superlu~mixedint+mpi~openmp+shared~superlu-dist~unified-memory cuda_arch=none patches=6e3336b1d62155f6350dfe42b0f9ea25d4fa0af60c7e540959139deb93a26059 arch=linux-rhel7-power9le
module load exasgd-hypre/2.20.0/openmpi-3.1.5/gcc-7.4.0-foavljr
# ipopt@3.12.10%gcc@7.4.0+coinhsl~debug~mumps arch=linux-rhel7-power9le
module load exasgd-ipopt/3.12.10/gcc-7.4.0-tj6jbm2
# libbsd@0.11.3%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-libbsd/0.11.3/gcc-7.4.0-4relejt
# libffi@3.3%gcc@7.4.0 patches=26f26c6f29a7ce9bf370ad3ab2610f99365b4bdd7b82e7c31df41a3370d685c0 arch=linux-rhel7-power9le
module load exasgd-libffi/3.3/gcc-7.4.0-6xqvxmj
# libiconv@1.16%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-libiconv/1.16/gcc-7.4.0-z5ga2ok
# libmd@1.0.3%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-libmd/1.0.3/gcc-7.4.0-dp4opcz
# libxml2@2.9.10%gcc@7.4.0~python arch=linux-rhel7-power9le
module load exasgd-libxml2/2.9.10/gcc-7.4.0-xl4bcrp
# magma@2.5.4%gcc@7.4.0+cuda+fortran~ipo+shared build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel7-power9le
# module load exasgd-magma/2.5.4/cuda-10.2.89-system/gcc-7.4.0-uypjvsx

# exasgd-magma/2.6.1/cuda-10.2.89-system/gcc-7.4.0-4bah4e6
module load exasgd-magma/2.6.1/cuda-10.2.89-system/gcc-7.4.0-4bah4e6
# metis@5.1.0%gcc@7.4.0~gdb~int64~real64+shared build_type=Release patches=4991da938c1d3a1d3dea78e49bbebecba00273f98df2a656e38b83d55b281da1,b1225da886605ea558db7ac08dd8054742ea5afe5ed61ad4d0fe7a495b1270d2 arch=linux-rhel7-power9le
module load exasgd-metis/5.1.0/gcc-7.4.0-7cjo5kb
# mpfr@4.1.0%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-mpfr/4.1.0/gcc-7.4.0-bkusb67
# ncurses@6.2%gcc@7.4.0~symlinks+termlib abi=none arch=linux-rhel7-power9le
module load exasgd-ncurses/6.2/gcc-7.4.0-rzufnf5
# openblas@0.3.15%gcc@7.4.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared threads=none arch=linux-rhel7-power9le
module load exasgd-openblas/0.3.15/gcc-7.4.0-i4ax3dw
# openssl@1.0.2k-fips%gcc@7.4.0~docs+systemcerts arch=linux-rhel7-power9le
module load exasgd-openssl/1.0.2k-fips/gcc-7.4.0-6izggwg
# parmetis@4.0.3%gcc@7.4.0~gdb~int64~ipo+shared build_type=RelWithDebInfo patches=4f892531eb0a807eb1b82e683a416d3e35154a455274cf9b162fb02054d11a5b,50ed2081bc939269689789942067c58b3e522c269269a430d5d34c00edbc5870,704b84f7c7444d4372cb59cca6e1209df4ef3b033bc4ee3cf50f369bce972a9d arch=linux-rhel7-power9le
module load exasgd-parmetis/4.0.3/openmpi-3.1.5/gcc-7.4.0-ym6nkdq
# perl@5.26.0%gcc@7.4.0+cpanm~shared~threads patches=0eac10ed90aeb0459ad8851f88081d439a4e41978e586ec743069e8b059370ac arch=linux-rhel7-power9le
module load exasgd-perl/5.26.0/gcc-7.4.0-j742qcl
# petsc@3.14.1%gcc@7.4.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw~giflib+hdf5~hip+hypre~int64~jpeg~knl~libpng~libyaml~memkind+metis~mkl-pardiso~moab~mpfr+mpi~mumps~p4est~ptscotch~random123~saws+shared~suite-sparse+superlu-dist~trilinos~valgrind amdgpu_target=none clanguage=C arch=linux-rhel7-power9le
module load exasgd-petsc/3.14.1/openmpi-3.1.5/gcc-7.4.0-243gohl
# pkgconf@1.7.4%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-pkgconf/1.7.4/gcc-7.4.0-5ios5aw
# python@3.8.6%gcc@7.4.0+bz2+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tix~tkinter~ucs4+uuid+zlib patches=0d98e93189bc278fbc37a50ed7f183bd8aaf249a8e1670a465f0db6bb4f8cf87 arch=linux-rhel7-power9le
module load exasgd-python/3.8.6/gcc-7.4.0-rs6cbgj
# raja@0.13.0%gcc@7.4.0+cuda~examples~exercises~hip~ipo+openmp+shared~tests amdgpu_target=none build_type=RelWithDebInfo cuda_arch=none arch=linux-rhel7-power9le
#module load exasgd-raja/0.13.0/cuda-10.2.89-system/gcc-7.4.0-haqrgbr
module load exasgd-raja/0.13.0/cuda-10.2.89-system/gcc-7.4.0-uj7zop4
# readline@8.1%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-readline/8.1/gcc-7.4.0-vqamj6w
# sqlite@3.35.5%gcc@7.4.0+column_metadata+fts~functions~rtree arch=linux-rhel7-power9le
module load exasgd-sqlite/3.35.5/gcc-7.4.0-nqk4jnf
# suite-sparse@5.8.1%gcc@7.4.0~cuda~openmp+pic~tbb arch=linux-rhel7-power9le
module load exasgd-suite-sparse/5.8.1/gcc-7.4.0-zggphid
# superlu-dist@6.4.0%gcc@7.4.0~cuda~int64~ipo~openmp+shared build_type=RelWithDebInfo cuda_arch=none arch=linux-rhel7-power9le
module load exasgd-superlu-dist/6.4.0/openmpi-3.1.5/gcc-7.4.0-bembyz7
# tar@1.34%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-tar/1.34/gcc-7.4.0-gkbivc5
# texinfo@6.5%gcc@7.4.0 patches=12f6edb0c6b270b8c8dba2ce17998c580db01182d871ee32b7b6e4129bd1d23a,1732115f651cff98989cb0215d8f64da5e0f7911ebf0c13b064920f088f2ffe1 arch=linux-rhel7-power9le
module load exasgd-texinfo/6.5/gcc-7.4.0-r24oezq
# umpire@4.1.2%gcc@7.4.0+c+cuda~deviceconst~examples~fortran~hip~ipo~numa+openmp~shared amdgpu_target=none build_type=RelWithDebInfo cuda_arch=none patches=7d912d31cd293df005ba74cb96c6f3e32dc3d84afff49b14509714283693db08 tests=none arch=linux-rhel7-power9le
#module load exasgd-umpire/4.1.2/cuda-10.2.89-system/gcc-7.4.0-hi2bbqe
# umpire@5.0.1%gcc@7.4.0+c+cuda~deviceconst~examples~fortran~hip~ipo~numa+openmp~shared amdgpu_target=none build_type=RelWithDebInfo cuda_arch=none patches=7d912d31cd293df005ba74cb96c6f3e32dc3d84afff49b14509714283693db08 tests=none arch=linux-rhel7-power9le
module load exasgd-umpire/5.0.1/cuda-10.2.89-system/gcc-7.4.0-mdo4fae
# Umpire 5.0.1 needs camp
module load exasgd-camp/0.1.0/cuda-10.2.89-system/gcc-7.4.0-fvkaniz

# util-linux-uuid@2.36.2%gcc@7.4.0 arch=linux-rhel7-power9le
module load exasgd-util-linux-uuid/2.36.2/gcc-7.4.0-oruzuyx
# xz@5.2.5%gcc@7.4.0~pic libs=shared,static arch=linux-rhel7-power9le
module load exasgd-xz/5.2.5/gcc-7.4.0-typzhoh
# zlib@1.2.11%gcc@7.4.0+optimize+pic+shared arch=linux-rhel7-power9le
module load exasgd-zlib/1.2.11/gcc-7.4.0-fkwiwuh

# Load system modules
module load gcc/7.4.0
module load openmpi/3.1.5
module load cuda/10.2

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

