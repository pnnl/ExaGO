module use -a /vast/projects/exasgd/spack/install/modules/linux-centos7-zen3
# pkgconf@1.8.0%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load pkgconf-1.8.0-clang-15.0.0-rocm5.3.0-lptbi7g
# ncurses@6.4%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~symlinks+termlib abi=none build_system=autotools arch=linux-centos7-zen
module load ncurses-6.4-clang-15.0.0-rocm5.3.0-wc7wnvv
# ca-certificates-mozilla@2023-01-10%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=generic arch=linux-centos7-zen
module load ca-certificates-mozilla-2023-01-10-clang-15.0.0-rocm5.3.0-6632r44
# perl@5.26.0%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +cpanm+open+shared+threads build_system=generic patches=0eac10e,8cf4302 arch=linux-centos7-zen
module load perl-5.26.0-clang-15.0.0-rocm5.3.0-igoczj4
# zlib@1.2.13%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +optimize+pic+shared build_system=makefile arch=linux-centos7-zen
module load zlib-1.2.13-clang-15.0.0-rocm5.3.0-pfa2yvw
# openssl@1.1.1t%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~docs~shared build_system=generic certs=mozilla arch=linux-centos7-zen
module load openssl-1.1.1t-clang-15.0.0-rocm5.3.0-kjvcx5a
# cmake@3.20.6%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~doc+ncurses+ownlibs~qt build_system=generic build_type=Release arch=linux-centos7-zen
module load cmake-3.20.6-clang-15.0.0-rocm5.3.0-h6jmxju
# blt@0.4.1%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=generic arch=linux-centos7-zen
module load blt-0.4.1-clang-15.0.0-rocm5.3.0-ohqfyq6
# gmake@4.4.1%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~guile build_system=autotools arch=linux-centos7-zen
module load gmake-4.4.1-clang-15.0.0-rocm5.3.0-x5krknj
# hip@5.3.0%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda~ipo+rocm build_system=cmake build_type=Release generator=make patches=ca523f1 arch=linux-centos7-zen
module load hip-5.3.0-clang-15.0.0-rocm5.3.0-flni4ra
# hsa-rocr-dev@5.3.0%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +image~ipo+shared build_system=cmake build_type=Release generator=make patches=71e6851 arch=linux-centos7-zen
module load hsa-rocr-dev-5.3.0-clang-15.0.0-rocm5.3.0-aoxomhi
# llvm-amdgpu@5.3.0%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~ipo~link_llvm_dylib~llvm_dylib~openmp+rocm-device-libs build_system=cmake build_type=Release generator=ninja patches=a08bbe1 arch=linux-centos7-zen
module load llvm-amdgpu-5.3.0-clang-15.0.0-rocm5.3.0-ofeglwz
# camp@0.2.3%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda~ipo~openmp+rocm~tests amdgpu_target=gfx908 build_system=cmake build_type=RelWithDebInfo generator=make arch=linux-centos7-zen
module load camp-0.2.3-clang-15.0.0-rocm5.3.0-skdmvd6
# openblas@0.3.20%gcc@8.4.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared build_system=makefile patches=9f12903 symbol_suffix=none threads=none arch=linux-centos7-zen
module load openblas-0.3.20-gcc-8.4.0-hhqjme4
# coinhsl@2019.05.21%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +blas build_system=autotools arch=linux-centos7-zen
module load coinhsl-2019.05.21-clang-15.0.0-rocm5.3.0-mimcunp
# hipblas@5.3.0%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda~ipo+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-centos7-zen
module load hipblas-5.3.0-clang-15.0.0-rocm5.3.0-3mepfxx
# hipsparse@5.3.0%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda~ipo+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make patches=c447537 arch=linux-centos7-zen
module load hipsparse-5.3.0-clang-15.0.0-rocm5.3.0-blmcuwn
# magma@2.6.2%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda+fortran~ipo+rocm+shared amdgpu_target=gfx908 build_system=cmake build_type=RelWithDebInfo generator=make arch=linux-centos7-zen
module load magma-2.6.2-clang-15.0.0-rocm5.3.0-h7ffwkf
# metis@5.1.0%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~gdb~int64~ipo~real64+shared build_system=cmake build_type=RelWithDebInfo generator=make patches=4991da9,93a7903 arch=linux-centos7-zen
module load metis-5.1.0-clang-15.0.0-rocm5.3.0-yvdpcgo
# openmpi@4.1.4%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~atomics~cuda~cxx~cxx_exceptions~gpfs~internal-hwloc~java~legacylaunchers~lustre~memchecker~orterunprefix+romio+rsh~singularity+static+vt+wrapper-rpath build_system=autotools fabrics=none schedulers=none arch=linux-centos7-zen
module load openmpi-4.1.4-clang-15.0.0-rocm5.3.0-dzsujir
# rocprim@5.3.0%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~ipo amdgpu_target=auto build_system=cmake build_type=Release generator=make patches=98506fb arch=linux-centos7-zen
module load rocprim-5.3.0-clang-15.0.0-rocm5.3.0-hhitx5k
# raja@0.14.0%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda~examples~exercises~ipo~openmp+rocm+shared~tests amdgpu_target=gfx908 build_system=cmake build_type=RelWithDebInfo generator=make arch=linux-centos7-zen
module load raja-0.14.0-clang-15.0.0-rocm5.3.0-4oicp3s
# libiconv@1.17%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools libs=shared,static arch=linux-centos7-zen
module load libiconv-1.17-clang-15.0.0-rocm5.3.0-bcqrqnt
# diffutils@3.9%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load diffutils-3.9-clang-15.0.0-rocm5.3.0-d7qy3yk
# libsigsegv@2.14%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load libsigsegv-2.14-clang-15.0.0-rocm5.3.0-2lrtpwx
# m4@1.4.19%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +sigsegv build_system=autotools patches=9dc5fbd,bfdffa7 arch=linux-centos7-zen
module load m4-1.4.19-clang-15.0.0-rocm5.3.0-m2hxqwv
# autoconf@2.69%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools patches=35c4492,7793209,a49dd5b arch=linux-centos7-zen
module load autoconf-2.69-clang-15.0.0-rocm5.3.0-svxd3jt
# automake@1.16.5%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load automake-1.16.5-clang-15.0.0-rocm5.3.0-e7iv77p
# libtool@2.4.7%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load libtool-2.4.7-clang-15.0.0-rocm5.3.0-dznqp7n
# gmp@6.2.1%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +cxx build_system=autotools libs=shared,static patches=69ad2e2 arch=linux-centos7-zen
module load gmp-6.2.1-clang-15.0.0-rocm5.3.0-vvrszuh
# autoconf-archive@2023.02.20%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load autoconf-archive-2023.02.20-clang-15.0.0-rocm5.3.0-vf36ebo
# bzip2@1.0.8%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~debug~pic+shared build_system=generic arch=linux-centos7-zen
module load bzip2-1.0.8-clang-15.0.0-rocm5.3.0-sn3czlz
# xz@5.4.1%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~pic build_system=autotools libs=shared,static arch=linux-centos7-zen
module load xz-5.4.1-clang-15.0.0-rocm5.3.0-rckubo2
# libxml2@2.10.3%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~python build_system=autotools arch=linux-centos7-zen
module load libxml2-2.10.3-clang-15.0.0-rocm5.3.0-z6hxybz
# pigz@2.7%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=makefile arch=linux-centos7-zen
module load pigz-2.7-clang-15.0.0-rocm5.3.0-w7ssmaq
# zstd@1.5.5%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +programs build_system=makefile compression=none libs=shared,static arch=linux-centos7-zen
module load zstd-1.5.5-clang-15.0.0-rocm5.3.0-baom5js
# tar@1.34%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools zip=pigz arch=linux-centos7-zen
module load tar-1.34-clang-15.0.0-rocm5.3.0-5dknjt2
# gettext@0.21.1%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +bzip2+curses+git~libunistring+libxml2+tar+xz build_system=autotools arch=linux-centos7-zen
module load gettext-0.21.1-clang-15.0.0-rocm5.3.0-yni5ply
# texinfo@7.0%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load texinfo-7.0-clang-15.0.0-rocm5.3.0-4rpqnes
# mpfr@4.2.0%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools libs=shared,static arch=linux-centos7-zen
module load mpfr-4.2.0-clang-15.0.0-rocm5.3.0-7gq36uq
# suite-sparse@5.13.0%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda~graphblas~openmp+pic build_system=generic arch=linux-centos7-zen
module load suite-sparse-5.13.0-clang-15.0.0-rocm5.3.0-qzvpnhs
# umpire@6.0.0%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +c~cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp+rocm+shared amdgpu_target=gfx908 build_system=cmake build_type=RelWithDebInfo generator=make tests=none arch=linux-centos7-zen
module load umpire-6.0.0-clang-15.0.0-rocm5.3.0-hsuiw34
# hiop@git.develop=0.7.2%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda+deepchecking~ginkgo~ipo~jsrun+kron+mpi+raja+rocm~shared+sparse amdgpu_target=gfx908 build_system=cmake build_type=RelWithDebInfo generator=make arch=linux-centos7-zen
module load hiop-git.develop=0.7.2-clang-15.0.0-rocm5.3.0-chceqoh
# ipopt@3.12.10%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +coinhsl~debug~metis~mumps build_system=autotools arch=linux-centos7-zen
module load ipopt-3.12.10-clang-15.0.0-rocm5.3.0-jnhzhnm
# python@3.7.0%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic arch=linux-centos7-zen
module load python-3.7.0-clang-15.0.0-rocm5.3.0-54xgmeh
# petsc@3.18.6%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~tetgen~trilinos~valgrind build_system=generic clanguage=C arch=linux-centos7-zen
module load petsc-3.18.6-clang-15.0.0-rocm5.3.0-jow3gvk
# exago@develop%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda+hiop~ipo+ipopt+logging+mpi~python+raja+rocm amdgpu_target=gfx908 build_system=cmake build_type=RelWithDebInfo dev_path=/vast/projects/exasgd/scratch/ruth521/exago generator=make arch=linux-centos7-zen
## module load exago-develop-clang-15.0.0-rocm5.3.0-d7flp76
# camp@0.2.3%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda~ipo~openmp+rocm~tests amdgpu_target=gfx908 build_system=cmake build_type=Release generator=make arch=linux-centos7-zen
module load camp-0.2.3-clang-15.0.0-rocm5.3.0-xayuovn
# magma@2.6.2%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda+fortran~ipo+rocm+shared amdgpu_target=gfx908 build_system=cmake build_type=Release generator=make arch=linux-centos7-zen
module load magma-2.6.2-clang-15.0.0-rocm5.3.0-failpgu
# raja@0.14.0%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda~examples~exercises~ipo~openmp+rocm+shared~tests amdgpu_target=gfx908 build_system=cmake build_type=Release generator=make arch=linux-centos7-zen
module load raja-0.14.0-clang-15.0.0-rocm5.3.0-x4u3jfh
# umpire@6.0.0%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +c~cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp+rocm+shared amdgpu_target=gfx908 build_system=cmake build_type=Release generator=make tests=none arch=linux-centos7-zen
module load umpire-6.0.0-clang-15.0.0-rocm5.3.0-tdgqmbo
# hiop@git.develop=0.7.2%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda~deepchecking~ginkgo~ipo~jsrun+kron+mpi+raja+rocm~shared+sparse amdgpu_target=gfx908 build_system=cmake build_type=Release generator=make arch=linux-centos7-zen
module load hiop-git.develop=0.7.2-clang-15.0.0-rocm5.3.0-bosbzo5
# exago@develop%clang@15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda+hiop~ipo+ipopt+logging+mpi~python+raja+rocm amdgpu_target=gfx908 build_system=cmake build_type=Release dev_path=/vast/projects/exasgd/scratch/ruth521/exago generator=make arch=linux-centos7-zen
## module load exago-develop-clang-15.0.0-rocm5.3.0-t7cuwxu
