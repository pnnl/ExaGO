module use -a /qfs/projects/exasgd/src/ci-incline/install/modules/linux-centos7-zen
# pkgconf@=1.9.5%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load pkgconf/1.9.5-clang-15.0.0-rocm5.3.0-l6m7w52
# nghttp2@=1.52.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load nghttp2/1.52.0-clang-15.0.0-rocm5.3.0-hcgjz26
# ca-certificates-mozilla@=2023-05-30%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=generic arch=linux-centos7-zen
module load ca-certificates-mozilla/2023-05-30-clang-15.0.0-rocm5.3.0-3h22nit
# perl@=5.26.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +cpanm+opcode+open+shared+threads build_system=generic patches=0eac10e,8cf4302 arch=linux-centos7-zen
module load perl/5.26.0-clang-15.0.0-rocm5.3.0-43pwa35
# zlib-ng@=2.1.3%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +compat+opt build_system=autotools arch=linux-centos7-zen
module load zlib-ng/2.1.3-clang-15.0.0-rocm5.3.0-nglc55e
# openssl@=3.1.2%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~docs+shared build_system=generic certs=mozilla arch=linux-centos7-zen
## module load openssl/3.1.2-clang-15.0.0-rocm5.3.0-myldvct
# curl@=8.1.2%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~gssapi~ldap~libidn2~librtmp~libssh~libssh2+nghttp2 build_system=autotools libs=shared,static tls=openssl arch=linux-centos7-zen
module load curl/8.1.2-clang-15.0.0-rocm5.3.0-oea3xhj
# ncurses@=6.4%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~symlinks+termlib abi=none build_system=autotools arch=linux-centos7-zen
module load ncurses/6.4-clang-15.0.0-rocm5.3.0-xfyzjrp
# cmake@=3.20.6%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~doc+ncurses+ownlibs build_system=generic build_type=Release arch=linux-centos7-zen
module load cmake/3.20.6-clang-15.0.0-rocm5.3.0-vfpr6me
# blt@=0.4.1%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=generic arch=linux-centos7-zen
module load blt/0.4.1-clang-15.0.0-rocm5.3.0-anjmmal
# gmake@=4.4.1%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~guile build_system=autotools arch=linux-centos7-zen
module load gmake/4.4.1-clang-15.0.0-rocm5.3.0-aasl32j
# hip@=5.3.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda+rocm build_system=cmake build_type=Release generator=make patches=c2ee21c,ca523f1 arch=linux-centos7-zen
module load hip/5.3.0-clang-15.0.0-rocm5.3.0-clc43s5
# hsa-rocr-dev@=5.3.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +image+shared build_system=cmake build_type=Release generator=make patches=71e6851 arch=linux-centos7-zen
module load hsa-rocr-dev/5.3.0-clang-15.0.0-rocm5.3.0-iq6jll3
# llvm-amdgpu@=5.3.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~link_llvm_dylib~llvm_dylib~openmp+rocm-device-libs build_system=cmake build_type=Release generator=ninja patches=a08bbe1 arch=linux-centos7-zen
module load llvm-amdgpu/5.3.0-clang-15.0.0-rocm5.3.0-ubc6533
# camp@=0.2.3%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda~ipo~openmp+rocm~tests amdgpu_target=gfx908 build_system=cmake build_type=Release generator=make arch=linux-centos7-zen
module load camp/0.2.3-clang-15.0.0-rocm5.3.0-zygugu6
# openblas@=0.3.20%gcc@=8.4.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared build_system=makefile patches=9f12903 symbol_suffix=none threads=none arch=linux-centos7-zen
module load openblas/0.3.20-gcc-8.4.0-jgwl2xj
# coinhsl@=2019.05.21%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +blas build_system=autotools arch=linux-centos7-zen
module load coinhsl/2019.05.21-clang-15.0.0-rocm5.3.0-ssnjpee
# hipblas@=5.3.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-centos7-zen
module load hipblas/5.3.0-clang-15.0.0-rocm5.3.0-d77frzb
# hipsparse@=5.3.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make patches=c447537 arch=linux-centos7-zen
module load hipsparse/5.3.0-clang-15.0.0-rocm5.3.0-ngj46ki
# magma@=2.6.2%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda+fortran~ipo+rocm+shared amdgpu_target=gfx908 build_system=cmake build_type=Release generator=make arch=linux-centos7-zen
module load magma/2.6.2-clang-15.0.0-rocm5.3.0-kmxsz43
# metis@=5.1.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~gdb~int64~ipo~real64+shared build_system=cmake build_type=Release generator=make patches=4991da9,93a7903 arch=linux-centos7-zen
module load metis/5.1.0-clang-15.0.0-rocm5.3.0-wauxd7r
# openmpi@=4.1.4%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~atomics~cuda~cxx~cxx_exceptions~gpfs~internal-hwloc~internal-pmix~java~legacylaunchers~lustre~memchecker~openshmem~orterunprefix+romio+rsh~singularity+static+vt+wrapper-rpath build_system=autotools fabrics=none schedulers=none arch=linux-centos7-zen
module load openmpi/4.1.4-clang-15.0.0-rocm5.3.0-2ry77mh
# rocprim@=5.3.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  amdgpu_target=auto build_system=cmake build_type=Release generator=make patches=98506fb arch=linux-centos7-zen
module load rocprim/5.3.0-clang-15.0.0-rocm5.3.0-z5i2hbr
# raja@=0.14.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda~examples~exercises~ipo~openmp+rocm+shared~tests amdgpu_target=gfx908 build_system=cmake build_type=Release generator=make arch=linux-centos7-zen
module load raja/0.14.0-clang-15.0.0-rocm5.3.0-cy57u2e
# libiconv@=1.17%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools libs=shared,static arch=linux-centos7-zen
module load libiconv/1.17-clang-15.0.0-rocm5.3.0-56oowkk
# diffutils@=3.9%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load diffutils/3.9-clang-15.0.0-rocm5.3.0-64ipopm
# libsigsegv@=2.14%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load libsigsegv/2.14-clang-15.0.0-rocm5.3.0-k2jw6uv
# m4@=1.4.19%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +sigsegv build_system=autotools patches=9dc5fbd,bfdffa7 arch=linux-centos7-zen
module load m4/1.4.19-clang-15.0.0-rocm5.3.0-a3pehkf
# autoconf@=2.69%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools patches=35c4492,7793209,a49dd5b arch=linux-centos7-zen
module load autoconf/2.69-clang-15.0.0-rocm5.3.0-yb25icj
# automake@=1.16.5%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load automake/1.16.5-clang-15.0.0-rocm5.3.0-mk7pljt
# libtool@=2.4.7%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load libtool/2.4.7-clang-15.0.0-rocm5.3.0-gnnghi4
# gmp@=6.2.1%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +cxx build_system=autotools libs=shared,static patches=69ad2e2 arch=linux-centos7-zen
module load gmp/6.2.1-clang-15.0.0-rocm5.3.0-spxx6ht
# autoconf-archive@=2023.02.20%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load autoconf-archive/2023.02.20-clang-15.0.0-rocm5.3.0-6qvs3k6
# bzip2@=1.0.8%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~debug~pic+shared build_system=generic arch=linux-centos7-zen
module load bzip2/1.0.8-clang-15.0.0-rocm5.3.0-agwvllw
# xz@=5.4.1%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~pic build_system=autotools libs=shared,static arch=linux-centos7-zen
module load xz/5.4.1-clang-15.0.0-rocm5.3.0-ke7l3hr
# libxml2@=2.10.3%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~python build_system=autotools arch=linux-centos7-zen
module load libxml2/2.10.3-clang-15.0.0-rocm5.3.0-f6bprlx
# pigz@=2.7%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=makefile arch=linux-centos7-zen
module load pigz/2.7-clang-15.0.0-rocm5.3.0-x3qqgg3
# zstd@=1.5.5%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +programs build_system=makefile compression=none libs=shared,static arch=linux-centos7-zen
module load zstd/1.5.5-clang-15.0.0-rocm5.3.0-lv6xa5q
# tar@=1.34%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools zip=pigz arch=linux-centos7-zen
module load tar/1.34-clang-15.0.0-rocm5.3.0-yyiph5k
# gettext@=0.21.1%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +bzip2+curses+git~libunistring+libxml2+tar+xz build_system=autotools arch=linux-centos7-zen
module load gettext/0.21.1-clang-15.0.0-rocm5.3.0-6miw6zm
# texinfo@=7.0.3%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools arch=linux-centos7-zen
module load texinfo/7.0.3-clang-15.0.0-rocm5.3.0-nbp2bo4
# mpfr@=4.2.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/"  build_system=autotools libs=shared,static arch=linux-centos7-zen
module load mpfr/4.2.0-clang-15.0.0-rocm5.3.0-f4qxift
# suite-sparse@=5.13.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda~graphblas~openmp+pic build_system=generic arch=linux-centos7-zen
module load suite-sparse/5.13.0-clang-15.0.0-rocm5.3.0-oyjwo5g
# umpire@=6.0.0%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +c~cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp+rocm+shared amdgpu_target=gfx908 build_system=cmake build_type=Release generator=make tests=none arch=linux-centos7-zen
module load umpire/6.0.0-clang-15.0.0-rocm5.3.0-xlgdssp
# hiop@=develop%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda~deepchecking~ginkgo~ipo~jsrun+kron+mpi+raja+rocm~shared+sparse amdgpu_target=gfx908 build_system=cmake build_type=Release generator=make arch=linux-centos7-zen
module load hiop/develop-clang-15.0.0-rocm5.3.0-msynpq2
# ipopt@=3.12.10%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +coinhsl~debug~metis~mumps build_system=autotools arch=linux-centos7-zen
module load ipopt/3.12.10-clang-15.0.0-rocm5.3.0-bc34jlg
# python@=3.11.4%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" +bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=13fa8bf,b0615b2,f2fd060 arch=linux-centos7-zen
module load python/3.11.4-clang-15.0.0-rocm5.3.0-n2tvv22
# petsc@=3.18.6%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~tetgen~trilinos~valgrind build_system=generic clanguage=C arch=linux-centos7-zen
module load petsc/3.18.6-clang-15.0.0-rocm5.3.0-jk3rfiw
# exago@=develop%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda+hiop~ipo+ipopt+logging+mpi~python+raja+rocm amdgpu_target=gfx908 build_system=cmake build_type=Release dev_path=/people/svcexasgd/gitlab/22333/spack_incline generator=make arch=linux-centos7-zen
## module load exago/develop-clang-15.0.0-rocm5.3.0-z2sfasr
