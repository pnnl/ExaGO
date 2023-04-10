module use -a /ccs/proj/csc359/spack-install/modules/linux-sles15-zen3
# libiconv@1.17%clang@14.0.0-rocm5.2.0 build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load libiconv-1.17-clang-14.0.0-rocm5.2.0-qx575m2
# diffutils@3.9%clang@14.0.0-rocm5.2.0 build_system=autotools arch=linux-sles15-zen3
module load diffutils-3.9-clang-14.0.0-rocm5.2.0-iymwfdp
# libsigsegv@2.14%clang@14.0.0-rocm5.2.0 build_system=autotools arch=linux-sles15-zen3
module load libsigsegv-2.14-clang-14.0.0-rocm5.2.0-nzvfsou
# m4@1.4.19%clang@14.0.0-rocm5.2.0+sigsegv build_system=autotools patches=9dc5fbd,bfdffa7 arch=linux-sles15-zen3
module load m4-1.4.19-clang-14.0.0-rocm5.2.0-5fr5trj
# perl@5.34.0%clang@14.0.0-rocm5.2.0+cpanm+open+shared+threads build_system=generic arch=linux-sles15-zen3
module load perl-5.34.0-clang-14.0.0-rocm5.2.0-irqjlkn
# autoconf@2.69%clang@14.0.0-rocm5.2.0 build_system=autotools patches=35c4492,7793209,a49dd5b arch=linux-sles15-zen3
module load autoconf-2.69-clang-14.0.0-rocm5.2.0-2vlbd33
# autoconf-archive@2023.02.20%clang@14.0.0-rocm5.2.0 build_system=autotools arch=linux-sles15-zen3
module load autoconf-archive-2023.02.20-clang-14.0.0-rocm5.2.0-orly6oh
# automake@1.16.5%clang@14.0.0-rocm5.2.0 build_system=autotools arch=linux-sles15-zen3
module load automake-1.16.5-clang-14.0.0-rocm5.2.0-rl2ekbh
# pkgconf@1.8.0%clang@14.0.0-rocm5.2.0 build_system=autotools arch=linux-sles15-zen3
module load pkgconf-1.8.0-clang-14.0.0-rocm5.2.0-7msa5ro
# ncurses@6.4%clang@14.0.0-rocm5.2.0~symlinks+termlib abi=none build_system=autotools arch=linux-sles15-zen3
module load ncurses-6.4-clang-14.0.0-rocm5.2.0-neuggfr
# ca-certificates-mozilla@2023-01-10%clang@14.0.0-rocm5.2.0 build_system=generic arch=linux-sles15-zen3
module load ca-certificates-mozilla-2023-01-10-clang-14.0.0-rocm5.2.0-qc2sxmw
# zlib@1.2.13%clang@14.0.0-rocm5.2.0+optimize+pic+shared build_system=makefile arch=linux-sles15-zen3
module load zlib-1.2.13-clang-14.0.0-rocm5.2.0-paehdy5
# openssl@1.1.1t%clang@14.0.0-rocm5.2.0~docs~shared build_system=generic certs=mozilla arch=linux-sles15-zen3
## module load openssl-1.1.1t-clang-14.0.0-rocm5.2.0-trx7c6p
# cmake@3.20.6%clang@14.0.0-rocm5.2.0~doc+ncurses+ownlibs~qt build_system=generic build_type=Release arch=linux-sles15-zen3
module load cmake-3.20.6-clang-14.0.0-rocm5.2.0-rfymesu
# blt@0.4.1%clang@14.0.0-rocm5.2.0 build_system=generic arch=linux-sles15-zen3
module load blt-0.4.1-clang-14.0.0-rocm5.2.0-lncshbc
# bzip2@1.0.8%clang@14.0.0-rocm5.2.0~debug~pic+shared build_system=generic arch=linux-sles15-zen3
module load bzip2-1.0.8-clang-14.0.0-rocm5.2.0-b723tqw
# gmake@4.4.1%clang@14.0.0-rocm5.2.0~guile build_system=autotools arch=linux-sles15-zen3
module load gmake-4.4.1-clang-14.0.0-rocm5.2.0-3lj6h5d
# hip@5.2.0%clang@14.0.0-rocm5.2.0~cuda~ipo+rocm build_system=cmake build_type=Release generator=make patches=959d1fe arch=linux-sles15-zen3
module load hip-5.2.0-clang-14.0.0-rocm5.2.0-t77exfq
# hsa-rocr-dev@5.2.0%clang@14.0.0-rocm5.2.0+image~ipo+shared build_system=cmake build_type=Release generator=make patches=71e6851 arch=linux-sles15-zen3
module load hsa-rocr-dev-5.2.0-clang-14.0.0-rocm5.2.0-inrdokj
# llvm-amdgpu@5.2.0%clang@14.0.0-rocm5.2.0~ipo~link_llvm_dylib~llvm_dylib~openmp+rocm-device-libs build_system=cmake build_type=Release generator=ninja patches=a08bbe1 arch=linux-sles15-zen3
module load llvm-amdgpu-5.2.0-clang-14.0.0-rocm5.2.0-dsqc66h
# camp@0.2.3%clang@14.0.0-rocm5.2.0~cuda~ipo~openmp+rocm~tests amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo generator=make arch=linux-sles15-zen3
module load camp-0.2.3-clang-14.0.0-rocm5.2.0-blsjn4h
# openblas@0.3.20%clang@14.0.0-rocm5.2.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared build_system=makefile patches=9f12903 symbol_suffix=none threads=none arch=linux-sles15-zen3
module load openblas-0.3.20-clang-14.0.0-rocm5.2.0-gpnoyst
# coinhsl@2019.05.21%clang@14.0.0-rocm5.2.0+blas build_system=autotools arch=linux-sles15-zen3
module load coinhsl-2019.05.21-clang-14.0.0-rocm5.2.0-hyeaxbz
# cray-mpich@8.1.23%clang@14.0.0-rocm5.2.0+wrappers build_system=generic arch=linux-sles15-zen3
module load cray-mpich-8.1.23-clang-14.0.0-rocm5.2.0-4bwvxwx
# xz@5.4.1%clang@14.0.0-rocm5.2.0~pic build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load xz-5.4.1-clang-14.0.0-rocm5.2.0-iwrlfwf
# libxml2@2.10.3%clang@14.0.0-rocm5.2.0~python build_system=autotools arch=linux-sles15-zen3
module load libxml2-2.10.3-clang-14.0.0-rocm5.2.0-ahoizys
# pigz@2.7%clang@14.0.0-rocm5.2.0 build_system=makefile arch=linux-sles15-zen3
module load pigz-2.7-clang-14.0.0-rocm5.2.0-iu4dbhw
# zstd@1.5.5%clang@14.0.0-rocm5.2.0+programs build_system=makefile compression=none libs=shared,static arch=linux-sles15-zen3
module load zstd-1.5.5-clang-14.0.0-rocm5.2.0-qknvav3
# tar@1.34%clang@14.0.0-rocm5.2.0 build_system=autotools zip=pigz arch=linux-sles15-zen3
module load tar-1.34-clang-14.0.0-rocm5.2.0-r7oesng
# gettext@0.21.1%clang@14.0.0-rocm5.2.0+bzip2+curses+git~libunistring+libxml2+tar+xz build_system=autotools arch=linux-sles15-zen3
module load gettext-0.21.1-clang-14.0.0-rocm5.2.0-imcqc2j
# hipblas@5.2.0%clang@14.0.0-rocm5.2.0~cuda~ipo+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load hipblas-5.2.0-clang-14.0.0-rocm5.2.0-4avahej
# hipsparse@5.2.0%clang@14.0.0-rocm5.2.0~cuda~ipo+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make patches=c447537 arch=linux-sles15-zen3
module load hipsparse-5.2.0-clang-14.0.0-rocm5.2.0-bokwexg
# rocprim@5.2.0%clang@14.0.0-rocm5.2.0~ipo amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load rocprim-5.2.0-clang-14.0.0-rocm5.2.0-doh7rym
# rocrand@5.2.0%clang@14.0.0-rocm5.2.0~ipo amdgpu_target=auto build_system=cmake build_type=Release generator=make patches=a35e689 arch=linux-sles15-zen3
module load rocrand-5.2.0-clang-14.0.0-rocm5.2.0-vozdoxx
# rocthrust@5.2.0%clang@14.0.0-rocm5.2.0~ipo amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load rocthrust-5.2.0-clang-14.0.0-rocm5.2.0-dyrvjhq
# ginkgo@1.5.0.glu_experimental%clang@14.0.0-rocm5.2.0~cuda~develtools~full_optimizations~hwloc~ipo~mpi~oneapi~openmp+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=Debug generator=make patches=ba0956e arch=linux-sles15-zen3
module load ginkgo-1.5.0.glu_experimental-clang-14.0.0-rocm5.2.0-ff6z3sw
# libtool@2.4.7%clang@14.0.0-rocm5.2.0 build_system=autotools arch=linux-sles15-zen3
module load libtool-2.4.7-clang-14.0.0-rocm5.2.0-67wv36s
# gmp@6.2.1%clang@14.0.0-rocm5.2.0+cxx build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load gmp-6.2.1-clang-14.0.0-rocm5.2.0-va37lbk
# magma@2.6.2%clang@14.0.0-rocm5.2.0~cuda+fortran~ipo+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo generator=make arch=linux-sles15-zen3
module load magma-2.6.2-clang-14.0.0-rocm5.2.0-epak2th
# metis@5.1.0%clang@14.0.0-rocm5.2.0~gdb~int64~ipo~real64+shared build_system=cmake build_type=RelWithDebInfo generator=make patches=4991da9,93a7903 arch=linux-sles15-zen3
module load metis-5.1.0-clang-14.0.0-rocm5.2.0-tyh6dse
# raja@0.14.0%clang@14.0.0-rocm5.2.0~cuda~examples~exercises~ipo~openmp+rocm+shared~tests amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo generator=make arch=linux-sles15-zen3
module load raja-0.14.0-clang-14.0.0-rocm5.2.0-e34ytuy
# texinfo@7.0%clang@14.0.0-rocm5.2.0 build_system=autotools arch=linux-sles15-zen3
module load texinfo-7.0-clang-14.0.0-rocm5.2.0-ihjde43
# mpfr@4.2.0%clang@14.0.0-rocm5.2.0 build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load mpfr-4.2.0-clang-14.0.0-rocm5.2.0-iio3dpj
# suite-sparse@5.13.0%clang@14.0.0-rocm5.2.0~cuda~graphblas~openmp+pic build_system=generic arch=linux-sles15-zen3
module load suite-sparse-5.13.0-clang-14.0.0-rocm5.2.0-7oo3se7
# umpire@6.0.0%clang@14.0.0-rocm5.2.0+c~cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo generator=make tests=none arch=linux-sles15-zen3
module load umpire-6.0.0-clang-14.0.0-rocm5.2.0-p5izspw
# hiop@0.7.2%clang@14.0.0-rocm5.2.0~cuda+deepchecking+ginkgo~ipo~jsrun+kron+mpi+raja+rocm~shared+sparse amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo generator=make arch=linux-sles15-zen3
module load hiop-0.7.2-clang-14.0.0-rocm5.2.0-qsoedde
# ipopt@3.12.10%clang@14.0.0-rocm5.2.0+coinhsl~debug~metis~mumps build_system=autotools arch=linux-sles15-zen3
module load ipopt-3.12.10-clang-14.0.0-rocm5.2.0-35vkate
# python@3.9.12%clang@14.0.0-rocm5.2.0+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=0d98e93,4c24573,f2fd060 arch=linux-sles15-zen3
module load python-3.9.12-clang-14.0.0-rocm5.2.0-w7cxnlm
# petsc@3.18.6%clang@14.0.0-rocm5.2.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~tetgen~trilinos~valgrind build_system=generic clanguage=C arch=linux-sles15-zen3
module load petsc-3.18.6-clang-14.0.0-rocm5.2.0-zgbxnan
