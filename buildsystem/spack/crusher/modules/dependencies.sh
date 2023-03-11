module use -a /ccs/proj/csc359/spack-install/modules/linux-sles15-zen3
# pkgconf@1.8.0%cce@14.0.2-rocm5.2.0 build_system=autotools arch=linux-sles15-zen3
module load pkgconf-1.8.0-cce-14.0.2-rocm5.2.0-qsabajr
# ncurses@6.4%cce@14.0.2-rocm5.2.0~symlinks+termlib abi=none build_system=autotools arch=linux-sles15-zen3
module load ncurses-6.4-cce-14.0.2-rocm5.2.0-qb7tif2
# ca-certificates-mozilla@2023-01-10%cce@14.0.2-rocm5.2.0 build_system=generic arch=linux-sles15-zen3
module load ca-certificates-mozilla-2023-01-10-cce-14.0.2-rocm5.2.0-u66ogzl
# perl@5.34.0%cce@14.0.2-rocm5.2.0+cpanm+open+shared+threads build_system=generic arch=linux-sles15-zen3
module load perl-5.34.0-cce-14.0.2-rocm5.2.0-rh7u7qc
# zlib@1.2.13%cce@14.0.2-rocm5.2.0+optimize+pic+shared build_system=makefile arch=linux-sles15-zen3
module load zlib-1.2.13-cce-14.0.2-rocm5.2.0-mp62w2l
# openssl@1.1.1t%cce@14.0.2-rocm5.2.0~docs~shared build_system=generic certs=mozilla arch=linux-sles15-zen3
## module load openssl-1.1.1t-cce-14.0.2-rocm5.2.0-lgn2jdq
# cmake@3.20.6%cce@14.0.2-rocm5.2.0~doc+ncurses+ownlibs~qt build_system=generic build_type=Release arch=linux-sles15-zen3
module load cmake-3.20.6-cce-14.0.2-rocm5.2.0-2fcpsfj
# blt@0.4.1%cce@14.0.2-rocm5.2.0 build_system=generic arch=linux-sles15-zen3
module load blt-0.4.1-cce-14.0.2-rocm5.2.0-7b5hxdq
# hip@5.2.0%cce@14.0.2-rocm5.2.0~cuda~ipo+rocm build_system=cmake build_type=Release patches=959d1fe arch=linux-sles15-zen3
module load hip-5.2.0-cce-14.0.2-rocm5.2.0-xwqfcoi
# hsa-rocr-dev@5.2.0%cce@14.0.2-rocm5.2.0+image~ipo+shared build_system=cmake build_type=Release patches=71e6851 arch=linux-sles15-zen3
module load hsa-rocr-dev-5.2.0-cce-14.0.2-rocm5.2.0-caj2mid
# llvm-amdgpu@5.2.0%cce@14.0.2-rocm5.2.0~ipo~link_llvm_dylib~llvm_dylib~openmp+rocm-device-libs build_system=cmake build_type=Release patches=a08bbe1 arch=linux-sles15-zen3
module load llvm-amdgpu-5.2.0-cce-14.0.2-rocm5.2.0-vicpz7u
# camp@0.2.3%cce@14.0.2-rocm5.2.0~cuda~ipo~openmp+rocm~tests amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo arch=linux-sles15-zen3
module load camp-0.2.3-cce-14.0.2-rocm5.2.0-zm5udxk
# cray-mpich@8.1.23%cce@14.0.2-rocm5.2.0+wrappers build_system=generic arch=linux-sles15-zen3
module load cray-mpich-8.1.23-cce-14.0.2-rocm5.2.0-hu32rmr
# openblas@0.3.20%gcc@11.2.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared build_system=makefile patches=114f95f,9f12903 symbol_suffix=none threads=none arch=linux-sles15-zen3
module load openblas-0.3.20-gcc-11.2.0-rito7eu
# coinhsl@2019.05.21%gcc@11.2.0+blas build_system=autotools arch=linux-sles15-zen3
module load coinhsl-2019.05.21-gcc-11.2.0-wqskyh2
# hipblas@5.2.0%cce@14.0.2-rocm5.2.0~ipo build_system=cmake build_type=Release arch=linux-sles15-zen3
module load hipblas-5.2.0-cce-14.0.2-rocm5.2.0-s5cqs75
# hipsparse@5.2.0%cce@14.0.2-rocm5.2.0~ipo build_system=cmake build_type=Release arch=linux-sles15-zen3
module load hipsparse-5.2.0-cce-14.0.2-rocm5.2.0-huw5pz3
# rocprim@5.2.0%cce@14.0.2-rocm5.2.0~ipo amdgpu_target=auto build_system=cmake build_type=Release arch=linux-sles15-zen3
module load rocprim-5.2.0-cce-14.0.2-rocm5.2.0-3jsm2j7
# rocrand@5.2.0%cce@14.0.2-rocm5.2.0~ipo amdgpu_target=auto build_system=cmake build_type=Release patches=a35e689 arch=linux-sles15-zen3
module load rocrand-5.2.0-cce-14.0.2-rocm5.2.0-7rjxdlm
# rocthrust@5.2.0%cce@14.0.2-rocm5.2.0~ipo amdgpu_target=auto build_system=cmake build_type=Release arch=linux-sles15-zen3
module load rocthrust-5.2.0-cce-14.0.2-rocm5.2.0-fr3tez2
# ginkgo@1.5.0.glu_experimental%cce@14.0.2-rocm5.2.0~cuda~develtools~full_optimizations~hwloc~ipo~mpi~oneapi~openmp+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=Debug patches=ba0956e arch=linux-sles15-zen3
module load ginkgo-1.5.0.glu_experimental-cce-14.0.2-rocm5.2.0-j5tytoh
# magma@2.6.2%cce@14.0.2-rocm5.2.0~cuda+fortran~ipo+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo arch=linux-sles15-zen3
module load magma-2.6.2-cce-14.0.2-rocm5.2.0-hobxmfc
# metis@5.1.0%cce@14.0.2-rocm5.2.0~gdb~int64~ipo~real64+shared build_system=cmake build_type=RelWithDebInfo patches=4991da9,93a7903 arch=linux-sles15-zen3
module load metis-5.1.0-cce-14.0.2-rocm5.2.0-7ssoebe
# raja@0.14.0%cce@14.0.2-rocm5.2.0~cuda+examples+exercises~ipo~openmp+rocm+shared~tests amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo arch=linux-sles15-zen3
module load raja-0.14.0-cce-14.0.2-rocm5.2.0-bysdt3z
# libiconv@1.17%cce@14.0.2-rocm5.2.0 build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load libiconv-1.17-cce-14.0.2-rocm5.2.0-xcxmbl6
# diffutils@3.8%cce@14.0.2-rocm5.2.0 build_system=autotools arch=linux-sles15-zen3
module load diffutils-3.8-cce-14.0.2-rocm5.2.0-fbdahwf
# libsigsegv@2.13%cce@14.0.2-rocm5.2.0 build_system=autotools arch=linux-sles15-zen3
module load libsigsegv-2.13-cce-14.0.2-rocm5.2.0-3td32vd
# m4@1.4.19%cce@14.0.2-rocm5.2.0+sigsegv build_system=autotools patches=9dc5fbd,bfdffa7 arch=linux-sles15-zen3
module load m4-1.4.19-cce-14.0.2-rocm5.2.0-dlvx6nt
# autoconf@2.69%cce@14.0.2-rocm5.2.0 build_system=autotools patches=35c4492,7793209,a49dd5b arch=linux-sles15-zen3
module load autoconf-2.69-cce-14.0.2-rocm5.2.0-cuo5cnh
# automake@1.16.5%cce@14.0.2-rocm5.2.0 build_system=autotools arch=linux-sles15-zen3
module load automake-1.16.5-cce-14.0.2-rocm5.2.0-op7j6iy
# libtool@2.4.7%cce@14.0.2-rocm5.2.0 build_system=autotools arch=linux-sles15-zen3
module load libtool-2.4.7-cce-14.0.2-rocm5.2.0-q5k4f3s
# gmp@6.2.1%cce@14.0.2-rocm5.2.0+cxx build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load gmp-6.2.1-cce-14.0.2-rocm5.2.0-wmshlq6
# autoconf-archive@2022.02.11%cce@14.0.2-rocm5.2.0 build_system=autotools patches=139214f arch=linux-sles15-zen3
module load autoconf-archive-2022.02.11-cce-14.0.2-rocm5.2.0-ayr7xya
# bzip2@1.0.8%cce@14.0.2-rocm5.2.0~debug~pic+shared build_system=generic arch=linux-sles15-zen3
module load bzip2-1.0.8-cce-14.0.2-rocm5.2.0-ix24bpa
# xz@5.4.1%cce@14.0.2-rocm5.2.0~pic build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load xz-5.4.1-cce-14.0.2-rocm5.2.0-hta7qck
# libxml2@2.10.3%cce@14.0.2-rocm5.2.0~python build_system=autotools arch=linux-sles15-zen3
module load libxml2-2.10.3-cce-14.0.2-rocm5.2.0-5mabs3o
# pigz@2.7%cce@14.0.2-rocm5.2.0 build_system=makefile arch=linux-sles15-zen3
module load pigz-2.7-cce-14.0.2-rocm5.2.0-gwoyvsn
# zstd@1.5.2%cce@14.0.2-rocm5.2.0+programs build_system=makefile compression=none libs=shared,static arch=linux-sles15-zen3
module load zstd-1.5.2-cce-14.0.2-rocm5.2.0-ziptnei
# tar@1.34%cce@14.0.2-rocm5.2.0 build_system=autotools zip=pigz arch=linux-sles15-zen3
module load tar-1.34-cce-14.0.2-rocm5.2.0-siw66qi
# gettext@0.21.1%cce@14.0.2-rocm5.2.0+bzip2+curses+git~libunistring+libxml2+tar+xz build_system=autotools arch=linux-sles15-zen3
module load gettext-0.21.1-cce-14.0.2-rocm5.2.0-seigj5f
# texinfo@7.0%cce@14.0.2-rocm5.2.0 build_system=autotools arch=linux-sles15-zen3
module load texinfo-7.0-cce-14.0.2-rocm5.2.0-azboikb
# mpfr@4.2.0%cce@14.0.2-rocm5.2.0 build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load mpfr-4.2.0-cce-14.0.2-rocm5.2.0-tvir4fd
# suite-sparse@5.13.0%cce@14.0.2-rocm5.2.0~cuda~graphblas~openmp+pic build_system=generic arch=linux-sles15-zen3
module load suite-sparse-5.13.0-cce-14.0.2-rocm5.2.0-kwdmixn
# umpire@6.0.0%cce@14.0.2-rocm5.2.0+c~cuda~device_alloc~deviceconst+examples~fortran~ipo~numa~openmp+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo tests=none arch=linux-sles15-zen3
module load umpire-6.0.0-cce-14.0.2-rocm5.2.0-xx5s3oq
# hiop@develop%cce@14.0.2-rocm5.2.0~cuda+deepchecking+ginkgo~ipo~jsrun+kron+mpi+raja+rocm~shared+sparse amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo arch=linux-sles15-zen3
module load hiop-develop-cce-14.0.2-rocm5.2.0-iye5vwm
# ipopt@3.12.10%gcc@11.2.0+coinhsl~debug~metis~mumps build_system=autotools arch=linux-sles15-zen3
module load ipopt-3.12.10-gcc-11.2.0-4re2kvj
# python@3.9.12%cce@14.0.2-rocm5.2.0+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=0d98e93,4c24573,f2fd060 arch=linux-sles15-zen3
module load python-3.9.12-cce-14.0.2-rocm5.2.0-bhb3tpd
# petsc@3.18.3%cce@14.0.2-rocm5.2.0~X+batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~tetgen~trilinos~valgrind build_system=generic clanguage=C arch=linux-sles15-zen3
module load petsc-3.18.3-cce-14.0.2-rocm5.2.0-k7xb2ln
# exago@develop%cce@14.0.2-rocm5.2.0~cuda+hiop~ipo+ipopt+logging+mpi~python+raja+rocm amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo dev_path=/ccs/home/rcruther/exago-git arch=linux-sles15-zen3
## module load exago-develop-cce-14.0.2-rocm5.2.0-sy7mipq
