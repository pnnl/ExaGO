module use -a /ccs/proj/csc359/spack-install/modules/linux-sles15-zen3
# pkgconf@1.8.0%gcc@12.2.0 build_system=autotools arch=linux-sles15-zen3
module load pkgconf-1.8.0-gcc-12.2.0-5odxteb
# ncurses@6.4%gcc@12.2.0~symlinks+termlib abi=none build_system=autotools arch=linux-sles15-zen3
module load ncurses-6.4-gcc-12.2.0-qk2frcm
# ca-certificates-mozilla@2023-01-10%gcc@12.2.0 build_system=generic arch=linux-sles15-zen3
module load ca-certificates-mozilla-2023-01-10-gcc-12.2.0-2bmxlfd
# perl@5.34.0%gcc@12.2.0+cpanm+open+shared+threads build_system=generic arch=linux-sles15-zen3
module load perl-5.34.0-gcc-12.2.0-67efboo
# zlib@1.2.13%gcc@12.2.0+optimize+pic+shared build_system=makefile arch=linux-sles15-zen3
module load zlib-1.2.13-gcc-12.2.0-ed3qpa6
# openssl@1.1.1t%gcc@12.2.0~docs~shared build_system=generic certs=mozilla arch=linux-sles15-zen3
module load openssl-1.1.1t-gcc-12.2.0-3o22prg
# cmake@3.20.6%gcc@12.2.0~doc+ncurses+ownlibs~qt build_system=generic build_type=Release arch=linux-sles15-zen3
module load cmake-3.20.6-gcc-12.2.0-l5dw2jd
# blt@0.4.1%gcc@12.2.0 build_system=generic arch=linux-sles15-zen3
module load blt-0.4.1-gcc-12.2.0-h43c2lz
# gmake@4.4.1%gcc@12.2.0~guile build_system=autotools arch=linux-sles15-zen3
module load gmake-4.4.1-gcc-12.2.0-wh5zvat
# hip@5.2.0%gcc@12.2.0~cuda~ipo+rocm build_system=cmake build_type=Release generator=make patches=959d1fe arch=linux-sles15-zen3
module load hip-5.2.0-gcc-12.2.0-5trwotd
# hsa-rocr-dev@5.2.0%gcc@12.2.0+image~ipo+shared build_system=cmake build_type=Release generator=make patches=71e6851 arch=linux-sles15-zen3
module load hsa-rocr-dev-5.2.0-gcc-12.2.0-dzalb3i
# llvm-amdgpu@5.2.0%gcc@12.2.0~ipo~link_llvm_dylib~llvm_dylib~openmp+rocm-device-libs build_system=cmake build_type=Release generator=ninja patches=a08bbe1 arch=linux-sles15-zen3
module load llvm-amdgpu-5.2.0-gcc-12.2.0-tdohwvc
# camp@0.2.3%gcc@12.2.0~cuda~ipo~openmp+rocm~tests amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo generator=make arch=linux-sles15-zen3
module load camp-0.2.3-gcc-12.2.0-qjbw2vb
# cray-mpich@8.1.23%gcc@12.2.0+wrappers build_system=generic arch=linux-sles15-zen3
module load cray-mpich-8.1.23-gcc-12.2.0-xx5cby4
# openblas@0.3.20%gcc@12.2.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared build_system=makefile patches=9f12903 symbol_suffix=none threads=none arch=linux-sles15-zen3
module load openblas-0.3.20-gcc-12.2.0-43klv2l
# coinhsl@2019.05.21%gcc@12.2.0+blas build_system=autotools arch=linux-sles15-zen3
module load coinhsl-2019.05.21-gcc-12.2.0-z6eyt5a
# hipblas@5.2.0%gcc@12.2.0~cuda~ipo+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load hipblas-5.2.0-gcc-12.2.0-ozk6aif
# hipsparse@5.2.0%gcc@12.2.0~cuda~ipo+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make patches=c447537 arch=linux-sles15-zen3
module load hipsparse-5.2.0-gcc-12.2.0-vltixvc
# rocprim@5.2.0%gcc@12.2.0~ipo amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load rocprim-5.2.0-gcc-12.2.0-pdfmtz4
# rocrand@5.2.0%gcc@12.2.0~ipo amdgpu_target=auto build_system=cmake build_type=Release generator=make patches=a35e689 arch=linux-sles15-zen3
module load rocrand-5.2.0-gcc-12.2.0-i4cp7uf
# rocthrust@5.2.0%gcc@12.2.0~ipo amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load rocthrust-5.2.0-gcc-12.2.0-t6tfrnz
# ginkgo@1.5.0.glu_experimental%gcc@12.2.0~cuda~develtools~full_optimizations~hwloc~ipo~mpi~oneapi~openmp+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=Debug generator=make patches=ba0956e arch=linux-sles15-zen3
module load ginkgo-1.5.0.glu_experimental-gcc-12.2.0-eypvuvy
# magma@2.6.2%gcc@12.2.0~cuda+fortran~ipo+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo generator=make arch=linux-sles15-zen3
module load magma-2.6.2-gcc-12.2.0-gcbqpni
# metis@5.1.0%gcc@12.2.0~gdb~int64~ipo~real64+shared build_system=cmake build_type=RelWithDebInfo generator=make patches=4991da9,93a7903,b1225da arch=linux-sles15-zen3
module load metis-5.1.0-gcc-12.2.0-atnlhlk
# raja@0.14.0%gcc@12.2.0~cuda+examples+exercises~ipo~openmp+rocm+shared~tests amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo generator=make arch=linux-sles15-zen3
module load raja-0.14.0-gcc-12.2.0-dhab7us
# libiconv@1.17%gcc@12.2.0 build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load libiconv-1.17-gcc-12.2.0-4bmagws
# diffutils@3.9%gcc@12.2.0 build_system=autotools arch=linux-sles15-zen3
module load diffutils-3.9-gcc-12.2.0-adift3b
# libsigsegv@2.14%gcc@12.2.0 build_system=autotools arch=linux-sles15-zen3
module load libsigsegv-2.14-gcc-12.2.0-ooxw47i
# m4@1.4.19%gcc@12.2.0+sigsegv build_system=autotools patches=9dc5fbd,bfdffa7 arch=linux-sles15-zen3
module load m4-1.4.19-gcc-12.2.0-3zudw4i
# autoconf@2.69%gcc@12.2.0 build_system=autotools patches=35c4492,7793209,a49dd5b arch=linux-sles15-zen3
module load autoconf-2.69-gcc-12.2.0-xeruta7
# automake@1.16.5%gcc@12.2.0 build_system=autotools arch=linux-sles15-zen3
module load automake-1.16.5-gcc-12.2.0-jylxxbx
# libtool@2.4.7%gcc@12.2.0 build_system=autotools arch=linux-sles15-zen3
module load libtool-2.4.7-gcc-12.2.0-a3btfk2
# gmp@6.2.1%gcc@12.2.0+cxx build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load gmp-6.2.1-gcc-12.2.0-ugj2poz
# autoconf-archive@2023.02.20%gcc@12.2.0 build_system=autotools arch=linux-sles15-zen3
module load autoconf-archive-2023.02.20-gcc-12.2.0-vmrgrxz
# bzip2@1.0.8%gcc@12.2.0~debug~pic+shared build_system=generic arch=linux-sles15-zen3
module load bzip2-1.0.8-gcc-12.2.0-xl55jjs
# xz@5.4.1%gcc@12.2.0~pic build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load xz-5.4.1-gcc-12.2.0-s5x3a52
# libxml2@2.10.3%gcc@12.2.0~python build_system=autotools arch=linux-sles15-zen3
module load libxml2-2.10.3-gcc-12.2.0-2fyp3uf
# pigz@2.7%gcc@12.2.0 build_system=makefile arch=linux-sles15-zen3
module load pigz-2.7-gcc-12.2.0-s6eby7t
# zstd@1.5.5%gcc@12.2.0+programs build_system=makefile compression=none libs=shared,static arch=linux-sles15-zen3
module load zstd-1.5.5-gcc-12.2.0-lois2cm
# tar@1.34%gcc@12.2.0 build_system=autotools zip=pigz arch=linux-sles15-zen3
module load tar-1.34-gcc-12.2.0-s7k435z
# gettext@0.21.1%gcc@12.2.0+bzip2+curses+git~libunistring+libxml2+tar+xz build_system=autotools arch=linux-sles15-zen3
module load gettext-0.21.1-gcc-12.2.0-a3gdlpb
# texinfo@7.0%gcc@12.2.0 build_system=autotools arch=linux-sles15-zen3
module load texinfo-7.0-gcc-12.2.0-uq4qysv
# mpfr@4.2.0%gcc@12.2.0 build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load mpfr-4.2.0-gcc-12.2.0-t6szsnf
# suite-sparse@5.13.0%gcc@12.2.0~cuda~graphblas~openmp+pic build_system=generic arch=linux-sles15-zen3
module load suite-sparse-5.13.0-gcc-12.2.0-fb7dtqq
# umpire@6.0.0%gcc@12.2.0+c~cuda~device_alloc~deviceconst+examples~fortran~ipo~numa~openmp+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo generator=make tests=none arch=linux-sles15-zen3
module load umpire-6.0.0-gcc-12.2.0-lroaukw
# hiop@git.develop=0.7.2%gcc@12.2.0~cuda+deepchecking+ginkgo~ipo~jsrun+kron+mpi+raja+rocm~shared+sparse amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo generator=make arch=linux-sles15-zen3
module load hiop-git.develop=0.7.2-gcc-12.2.0-g4abspx
# ipopt@3.12.10%gcc@12.2.0+coinhsl+debug~metis~mumps build_system=autotools arch=linux-sles15-zen3
module load ipopt-3.12.10-gcc-12.2.0-natgn77
# python@3.9.12%gcc@12.2.0+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=0d98e93,4c24573,f2fd060 arch=linux-sles15-zen3
module load python-3.9.12-gcc-12.2.0-wcaect6
# petsc@3.18.6%gcc@12.2.0~X+batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~tetgen~trilinos~valgrind build_system=generic clanguage=C arch=linux-sles15-zen3
module load petsc-3.18.6-gcc-12.2.0-a6l4ghx
# exago@develop%gcc@12.2.0~cuda+hiop~ipo+ipopt+logging+mpi~python+raja+rocm amdgpu_target=gfx90a build_system=cmake build_type=RelWithDebInfo dev_path=/ccs/home/rcruther/exago-git generator=make arch=linux-sles15-zen3
## module load exago-develop-gcc-12.2.0-n5piuea
# camp@0.2.3%gcc@12.2.0~cuda~ipo~openmp+rocm~tests amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load camp-0.2.3-gcc-12.2.0-36e4p4t
# ginkgo@1.5.0.glu_experimental%gcc@12.2.0~cuda~develtools~full_optimizations~hwloc~ipo~mpi~oneapi~openmp+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make patches=ba0956e arch=linux-sles15-zen3
module load ginkgo-1.5.0.glu_experimental-gcc-12.2.0-jhukxse
# magma@2.6.2%gcc@12.2.0~cuda+fortran~ipo+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load magma-2.6.2-gcc-12.2.0-gmnuwn4
# raja@0.14.0%gcc@12.2.0~cuda+examples+exercises~ipo~openmp+rocm+shared~tests amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load raja-0.14.0-gcc-12.2.0-qtf27vk
# umpire@6.0.0%gcc@12.2.0+c~cuda~device_alloc~deviceconst+examples~fortran~ipo~numa~openmp+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make tests=none arch=linux-sles15-zen3
module load umpire-6.0.0-gcc-12.2.0-el2erch
# hiop@git.develop=0.7.2%gcc@12.2.0~cuda~deepchecking+ginkgo~ipo~jsrun+kron+mpi+raja+rocm~shared+sparse amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load hiop-git.develop=0.7.2-gcc-12.2.0-gikfw6w
# exago@develop%gcc@12.2.0~cuda+hiop~ipo+ipopt+logging+mpi~python+raja+rocm amdgpu_target=gfx90a build_system=cmake build_type=Release dev_path=/ccs/home/rcruther/exago-git generator=make arch=linux-sles15-zen3
## module load exago-develop-gcc-12.2.0-zihzuex
