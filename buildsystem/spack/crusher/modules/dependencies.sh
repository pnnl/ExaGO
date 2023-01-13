module use -a /gpfs/alpine/proj-shared/csc359/cameron/spack-install/test-modules/linux-sles15-zen3
# python@3.9.12%clang@14.0.2+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=0d98e93,4c24573,f2fd060 arch=linux-sles15-zen3
module load python-3.9.12-clang-14.0.2-ullptbr
# amdblis@3.2%clang@14.0.2+blas+cblas~ilp64 build_system=makefile libs=shared,static threads=openmp arch=linux-sles15-zen3
module load amdblis-3.2-clang-14.0.2-rxubipx
# pkgconf@1.8.0%clang@14.0.2 build_system=autotools arch=linux-sles15-zen3
module load pkgconf-1.8.0-clang-14.0.2-y2zkflv
# ncurses@6.3%clang@14.0.2~symlinks+termlib abi=none build_system=autotools arch=linux-sles15-zen3
module load ncurses-6.3-clang-14.0.2-hnyfwsu
# ca-certificates-mozilla@2022-10-11%clang@14.0.2 build_system=generic arch=linux-sles15-zen3
module load ca-certificates-mozilla-2022-10-11-clang-14.0.2-dxhn6js
# berkeley-db@18.1.40%clang@14.0.2+cxx~docs+stl build_system=autotools patches=26090f4,b231fcc arch=linux-sles15-zen3
module load berkeley-db-18.1.40-clang-14.0.2-z54xfrx
# libiconv@1.16%gcc@11.2.0 build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load libiconv-1.16-gcc-11.2.0-xfogkcu
# diffutils@3.8%clang@14.0.2 build_system=autotools arch=linux-sles15-zen3
module load diffutils-3.8-clang-14.0.2-aupdpnh
# bzip2@1.0.8%clang@14.0.2~debug~pic+shared build_system=generic arch=linux-sles15-zen3
module load bzip2-1.0.8-clang-14.0.2-wxqbbqw
# readline@8.2%clang@14.0.2 build_system=autotools patches=bbf97f1 arch=linux-sles15-zen3
module load readline-8.2-clang-14.0.2-y5ghqml
# gdbm@1.23%clang@14.0.2 build_system=autotools arch=linux-sles15-zen3
module load gdbm-1.23-clang-14.0.2-b7qjqnr
# zlib@1.2.13%clang@14.0.2+optimize+pic+shared build_system=makefile arch=linux-sles15-zen3
module load zlib-1.2.13-clang-14.0.2-qiyb7ui
# perl@5.36.0%clang@14.0.2+cpanm+shared+threads build_system=generic arch=linux-sles15-zen3
module load perl-5.36.0-clang-14.0.2-efnilou
# openssl@1.1.1s%clang@14.0.2~docs~shared build_system=generic certs=mozilla arch=linux-sles15-zen3
## module load openssl-1.1.1s-clang-14.0.2-vmges7t
# cmake@3.20.6%clang@14.0.2~doc+ncurses+ownlibs~qt build_system=generic build_type=Release arch=linux-sles15-zen3
module load cmake-3.20.6-clang-14.0.2-2s4lu26
# blt@0.4.1%clang@14.0.2 build_system=generic arch=linux-sles15-zen3
module load blt-0.4.1-clang-14.0.2-64jstoj
# hip@5.2.0%clang@14.0.2~ipo build_system=cmake build_type=Release patches=959d1fe arch=linux-sles15-zen3
module load hip-5.2.0-clang-14.0.2-bd2h4ej
# hsa-rocr-dev@5.2.0%clang@14.0.2+image~ipo+shared build_system=cmake build_type=Release patches=71e6851 arch=linux-sles15-zen3
module load hsa-rocr-dev-5.2.0-clang-14.0.2-2ujoy4f
# llvm-amdgpu@5.2.0%clang@14.0.2~ipo~link_llvm_dylib~llvm_dylib~openmp+rocm-device-libs build_system=cmake build_type=Release patches=a08bbe1 arch=linux-sles15-zen3
module load llvm-amdgpu-5.2.0-clang-14.0.2-j2ldgoy
# camp@0.2.3%clang@14.0.2~cuda~ipo~openmp+rocm~tests amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo arch=linux-sles15-zen3
module load camp-0.2.3-clang-14.0.2-zvzuiwe
# cray-mpich@8.1.23%gcc@11.2.0+wrappers build_system=generic arch=linux-sles15-zen3
module load cray-mpich-8.1.23-gcc-11.2.0-nnbdyog
# coinhsl@2019.05.21%clang@14.0.2+blas build_system=autotools arch=linux-sles15-zen3
module load coinhsl-2019.05.21-clang-14.0.2-577fc5o
# hipblas@5.2.0%clang@14.0.2~ipo build_system=cmake build_type=Release arch=linux-sles15-zen3
module load hipblas-5.2.0-clang-14.0.2-qqixtuf
# hipsparse@5.2.0%clang@14.0.2~ipo build_system=cmake build_type=Release arch=linux-sles15-zen3
module load hipsparse-5.2.0-clang-14.0.2-7mkrl5q
# rocprim@5.2.0%clang@14.0.2~ipo amdgpu_target=auto build_system=cmake build_type=Release arch=linux-sles15-zen3
module load rocprim-5.2.0-clang-14.0.2-q6h7yze
# rocrand@5.2.0%clang@14.0.2~ipo amdgpu_target=auto build_system=cmake build_type=Release patches=a35e689 arch=linux-sles15-zen3
module load rocrand-5.2.0-clang-14.0.2-7f2icji
# rocthrust@5.2.0%clang@14.0.2~ipo amdgpu_target=auto build_system=cmake build_type=Release arch=linux-sles15-zen3
module load rocthrust-5.2.0-clang-14.0.2-4fj4kzx
# ginkgo@1.5.0.glu_experimental%clang@14.0.2~cuda~develtools~full_optimizations~hwloc~ipo+mpi~oneapi+openmp+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=Debug patches=ba0956e arch=linux-sles15-zen3
module load ginkgo-1.5.0.glu_experimental-clang-14.0.2-wuvtnl5
# libflame@5.2.0%clang@14.0.2~debug+lapack2flame+shared+static build_system=autotools patches=bf75a5b,c5cae9e threads=none arch=linux-sles15-zen3
module load libflame-5.2.0-clang-14.0.2-alv4ovq
# magma@2.6.2%clang@14.0.2~cuda+fortran~ipo+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo arch=linux-sles15-zen3
module load magma-2.6.2-clang-14.0.2-c3frmxe
# metis@5.1.0%clang@14.0.2~gdb~int64~ipo~real64+shared build_system=cmake build_type=RelWithDebInfo patches=4991da9,93a7903 arch=linux-sles15-zen3
module load metis-5.1.0-clang-14.0.2-og6d6pl
# raja@0.14.0%clang@14.0.2~cuda~examples~exercises~ipo~openmp+rocm+shared~tests amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo arch=linux-sles15-zen3
module load raja-0.14.0-clang-14.0.2-hdnlv2f
# gmp@6.2.1%clang@14.0.2 build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load gmp-6.2.1-clang-14.0.2-judyqxq
# libsigsegv@2.13%clang@14.0.2 build_system=autotools arch=linux-sles15-zen3
module load libsigsegv-2.13-clang-14.0.2-h5n3cl3
# m4@1.4.19%clang@14.0.2+sigsegv build_system=autotools patches=9dc5fbd,bfdffa7 arch=linux-sles15-zen3
module load m4-1.4.19-clang-14.0.2-7ydmz4n
# autoconf@2.69%clang@14.0.2 build_system=autotools patches=35c4492,7793209,a49dd5b arch=linux-sles15-zen3
module load autoconf-2.69-clang-14.0.2-wgebbb5
# autoconf-archive@2022.02.11%clang@14.0.2 build_system=autotools patches=139214f arch=linux-sles15-zen3
module load autoconf-archive-2022.02.11-clang-14.0.2-xcr5olb
# automake@1.16.5%clang@14.0.2 build_system=autotools arch=linux-sles15-zen3
module load automake-1.16.5-clang-14.0.2-bvhew4s
# libtool@2.4.7%clang@14.0.2 build_system=autotools arch=linux-sles15-zen3
module load libtool-2.4.7-clang-14.0.2-lej6cj7
# texinfo@7.0%clang@14.0.2 build_system=autotools arch=linux-sles15-zen3
module load texinfo-7.0-clang-14.0.2-qinfgon
# mpfr@4.1.0%clang@14.0.2 build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load mpfr-4.1.0-clang-14.0.2-yh2d4ju
# suite-sparse@5.13.0%clang@14.0.2~cuda~graphblas~openmp+pic build_system=generic arch=linux-sles15-zen3
module load suite-sparse-5.13.0-clang-14.0.2-6yyi42p
# umpire@6.0.0%clang@14.0.2+c~cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo tests=none arch=linux-sles15-zen3
module load umpire-6.0.0-clang-14.0.2-4o5bvyt
# hiop@0.7.1%clang@14.0.2~cuda+deepchecking+ginkgo~ipo~jsrun+kron+mpi+raja+rocm~shared+sparse amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo arch=linux-sles15-zen3
module load hiop-0.7.1-clang-14.0.2-bgcu6dv
# ipopt@3.12.10%clang@14.0.2+coinhsl~debug~metis~mumps build_system=autotools arch=linux-sles15-zen3
module load ipopt-3.12.10-clang-14.0.2-3m7fado
# petsc@3.18.3%clang@14.0.2~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~tetgen~trilinos~valgrind build_system=generic clanguage=C arch=linux-sles15-zen3
module load petsc-3.18.3-clang-14.0.2-o6urmxv
# exago@develop%clang@14.0.2~cuda+hiop~ipo+ipopt+mpi~python+raja+rocm amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo dev_path=/ccs/home/rcruther/exago-git arch=linux-sles15-zen3
## module load exago-develop-clang-14.0.2-umknd2x
