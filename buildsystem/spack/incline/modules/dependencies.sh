module use -a /qfs/projects/exasgd/src/ci-incline/install/modules/linux-centos7-zen
# gmake@=4.4.1%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~guile build_system=generic arch=linux-centos7-zen
module load gmake/4.4.1-clang-15.0.0-rocm5.3.0-gzvkpew
# pkgconf@=1.9.5%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load pkgconf/1.9.5-clang-15.0.0-rocm5.3.0-7jofo2w
# nghttp2@=1.57.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load nghttp2/1.57.0-clang-15.0.0-rocm5.3.0-fik2ce4
# ca-certificates-mozilla@=2023-05-30%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=generic arch=linux-centos7-zen
module load ca-certificates-mozilla/2023-05-30-clang-15.0.0-rocm5.3.0-3h22nit
# perl@=5.26.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +cpanm+opcode+open+shared+threads build_system=generic patches=0eac10e,8cf4302 arch=linux-centos7-zen
module load perl/5.26.0-clang-15.0.0-rocm5.3.0-caldjac
# zlib-ng@=2.1.5%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +compat+opt build_system=autotools arch=linux-centos7-zen
module load zlib-ng/2.1.5-clang-15.0.0-rocm5.3.0-526o4rf
# openssl@=3.1.3%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~docs+shared build_system=generic certs=mozilla arch=linux-centos7-zen
## module load openssl/3.1.3-clang-15.0.0-rocm5.3.0-i7yqqwi
# curl@=8.4.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~gssapi~ldap~libidn2~librtmp~libssh~libssh2+nghttp2 build_system=autotools libs=shared,static tls=openssl arch=linux-centos7-zen
module load curl/8.4.0-clang-15.0.0-rocm5.3.0-jsvdhcc
# ncurses@=6.4%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~symlinks+termlib abi=none build_system=autotools arch=linux-centos7-zen
module load ncurses/6.4-clang-15.0.0-rocm5.3.0-w7n6gfi
# cmake@=3.20.6%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~doc+ncurses+ownlibs build_system=generic build_type=Release arch=linux-centos7-zen
module load cmake/3.20.6-clang-15.0.0-rocm5.3.0-o42dvyt
# blt@=0.4.1%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=generic arch=linux-centos7-zen
module load blt/0.4.1-clang-15.0.0-rocm5.3.0-ghyk7hz
# hip@=5.3.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda+rocm build_system=cmake build_type=Release generator=make patches=c2ee21c,ca523f1 arch=linux-centos7-zen
module load hip/5.3.0-clang-15.0.0-rocm5.3.0-h2m6ftv
# hsa-rocr-dev@=5.3.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +image+shared build_system=cmake build_type=Release generator=make patches=9267179 arch=linux-centos7-zen
module load hsa-rocr-dev/5.3.0-clang-15.0.0-rocm5.3.0-6u5zxpw
# llvm-amdgpu@=5.3.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~link_llvm_dylib~llvm_dylib~openmp+rocm-device-libs build_system=cmake build_type=Release generator=ninja patches=a08bbe1 arch=linux-centos7-zen
module load llvm-amdgpu/5.3.0-clang-15.0.0-rocm5.3.0-yggd56o
# camp@=0.2.3%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda~ipo~openmp+rocm~tests amdgpu_target=gfx908 build_system=cmake build_type=Release generator=make arch=linux-centos7-zen
module load camp/0.2.3-clang-15.0.0-rocm5.3.0-lukemzi
# gmake@=4.4.1%gcc@=8.4.0~guile build_system=generic arch=linux-centos7-zen
module load gmake/4.4.1-gcc-8.4.0-4fgniv2
# openblas@=0.3.20%gcc@=8.4.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared build_system=makefile patches=9f12903 symbol_suffix=none threads=none arch=linux-centos7-zen
module load openblas/0.3.20-gcc-8.4.0-ii6ckx5
# coinhsl@=2019.05.21%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +blas build_system=autotools arch=linux-centos7-zen
module load coinhsl/2019.05.21-clang-15.0.0-rocm5.3.0-o5rf7bo
# hipblas@=5.3.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-centos7-zen
module load hipblas/5.3.0-clang-15.0.0-rocm5.3.0-zwlmr5u
# hipsparse@=5.3.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make patches=c447537 arch=linux-centos7-zen
module load hipsparse/5.3.0-clang-15.0.0-rocm5.3.0-ngj46ki
# magma@=2.6.2%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda+fortran~ipo+rocm+shared amdgpu_target=gfx908 build_system=cmake build_type=Release generator=make arch=linux-centos7-zen
module load magma/2.6.2-clang-15.0.0-rocm5.3.0-zc53wqt
# metis@=5.1.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~gdb~int64~ipo~real64+shared build_system=cmake build_type=Release generator=make patches=4991da9,93a7903 arch=linux-centos7-zen
module load metis/5.1.0-clang-15.0.0-rocm5.3.0-4vgjgia
# openmpi@=4.1.4%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~atomics~cuda~cxx~cxx_exceptions~gpfs~internal-hwloc~internal-pmix~java~legacylaunchers~lustre~memchecker~openshmem~orterunprefix+romio+rsh~singularity+static+vt+wrapper-rpath build_system=autotools fabrics=none schedulers=none arch=linux-centos7-zen
module load openmpi/4.1.4-clang-15.0.0-rocm5.3.0-2ry77mh
# rocprim@=5.3.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  amdgpu_target=auto build_system=cmake build_type=Release generator=make patches=98506fb arch=linux-centos7-zen
module load rocprim/5.3.0-clang-15.0.0-rocm5.3.0-z5i2hbr
# raja@=0.14.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda~examples~exercises~ipo~openmp~plugins+rocm+shared~tests amdgpu_target=gfx908 build_system=cmake build_type=Release generator=make arch=linux-centos7-zen
module load raja/0.14.0-clang-15.0.0-rocm5.3.0-mubleoh
# libiconv@=1.17%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools libs=shared,static arch=linux-centos7-zen
module load libiconv/1.17-clang-15.0.0-rocm5.3.0-anlnmwt
# diffutils@=3.9%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load diffutils/3.9-clang-15.0.0-rocm5.3.0-q66q3nb
# libsigsegv@=2.14%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load libsigsegv/2.14-clang-15.0.0-rocm5.3.0-7ervaqq
# m4@=1.4.19%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +sigsegv build_system=autotools patches=9dc5fbd,bfdffa7 arch=linux-centos7-zen
module load m4/1.4.19-clang-15.0.0-rocm5.3.0-q35v7cq
# autoconf@=2.69%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools patches=35c4492,7793209,a49dd5b arch=linux-centos7-zen
module load autoconf/2.69-clang-15.0.0-rocm5.3.0-ykvmhrb
# automake@=1.16.5%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load automake/1.16.5-clang-15.0.0-rocm5.3.0-w5anlbv
# libtool@=2.4.7%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load libtool/2.4.7-clang-15.0.0-rocm5.3.0-bj5s5e5
# gmp@=6.2.1%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +cxx build_system=autotools libs=shared,static patches=69ad2e2 arch=linux-centos7-zen
module load gmp/6.2.1-clang-15.0.0-rocm5.3.0-kshhft2
# autoconf-archive@=2023.02.20%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load autoconf-archive/2023.02.20-clang-15.0.0-rocm5.3.0-of2i77c
# bzip2@=1.0.8%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~debug~pic+shared build_system=generic arch=linux-centos7-zen
module load bzip2/1.0.8-clang-15.0.0-rocm5.3.0-l3whoev
# xz@=5.4.1%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~pic build_system=autotools libs=shared,static arch=linux-centos7-zen
module load xz/5.4.1-clang-15.0.0-rocm5.3.0-mvnfwal
# libxml2@=2.10.3%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +pic~python+shared build_system=autotools arch=linux-centos7-zen
module load libxml2/2.10.3-clang-15.0.0-rocm5.3.0-3yavliv
# pigz@=2.7%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=makefile arch=linux-centos7-zen
module load pigz/2.7-clang-15.0.0-rocm5.3.0-w6hqfb7
# zstd@=1.5.5%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +programs build_system=makefile compression=none libs=shared,static arch=linux-centos7-zen
module load zstd/1.5.5-clang-15.0.0-rocm5.3.0-5mdcukt
# tar@=1.34%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools zip=pigz arch=linux-centos7-zen
module load tar/1.34-clang-15.0.0-rocm5.3.0-wrpvflu
# gettext@=0.22.4%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +bzip2+curses+git~libunistring+libxml2+pic+shared+tar+xz build_system=autotools arch=linux-centos7-zen
module load gettext/0.22.4-clang-15.0.0-rocm5.3.0-cb5rwnh
# texinfo@=7.0.3%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load texinfo/7.0.3-clang-15.0.0-rocm5.3.0-khbvgj6
# mpfr@=4.2.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools libs=shared,static arch=linux-centos7-zen
module load mpfr/4.2.0-clang-15.0.0-rocm5.3.0-sg4eb4y
# suite-sparse@=5.13.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda~graphblas~openmp+pic build_system=generic arch=linux-centos7-zen
module load suite-sparse/5.13.0-clang-15.0.0-rocm5.3.0-e2gf3ub
# umpire@=6.0.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +c~cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp+rocm+shared amdgpu_target=gfx908 build_system=cmake build_type=Release generator=make tests=none arch=linux-centos7-zen
module load umpire/6.0.0-clang-15.0.0-rocm5.3.0-jlfcdqu
# hiop@=develop%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda~deepchecking~ginkgo~ipo~jsrun+kron+mpi+raja+rocm~shared+sparse amdgpu_target=gfx908 build_system=cmake build_type=Release dev_path=/people/svcexasgd/gitlab/26536/spack_incline/hiop_dev generator=make arch=linux-centos7-zen
module load hiop/develop-clang-15.0.0-rocm5.3.0-rpx6cfv
# ipopt@=3.12.10%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +coinhsl~debug~metis~mumps build_system=autotools arch=linux-centos7-zen
module load ipopt/3.12.10-clang-15.0.0-rocm5.3.0-fyusbut
# python@=3.11.4%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=13fa8bf,b0615b2,ebdca64,f2fd060 arch=linux-centos7-zen
module load python/3.11.4-clang-15.0.0-rocm5.3.0-bjornb6
# petsc@=3.20.1%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~sycl~tetgen~trilinos~valgrind build_system=generic clanguage=C memalign=none arch=linux-centos7-zen
module load petsc/3.20.1-clang-15.0.0-rocm5.3.0-qwd245y
# exago@=develop%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda+hiop~ipo+ipopt+logging+mpi~python+raja+rocm amdgpu_target=gfx908 build_system=cmake build_type=Release dev_path=/people/svcexasgd/gitlab/26536/spack_incline generator=make arch=linux-centos7-zen
## module load exago/develop-clang-15.0.0-rocm5.3.0-pnku3g5
