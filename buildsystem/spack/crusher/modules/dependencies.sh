module use -a /lustre/orion/eng145/world-shared/spack-install/modules/linux-sles15-x86_64
# gmake@=4.3%clang@=17.0.0-rocm5.7.1-mixed~guile build_system=generic patches=599f134 arch=linux-sles15-x86_64
module load gmake/4.3-clang-17.0.0-rocm5.7.1-mixed-x2fskbm
# pkgconf@=1.9.5%clang@=17.0.0-rocm5.7.1-mixed build_system=autotools arch=linux-sles15-x86_64
module load pkgconf/1.9.5-clang-17.0.0-rocm5.7.1-mixed-m4urn4t
# nghttp2@=1.57.0%clang@=17.0.0-rocm5.7.1-mixed build_system=autotools arch=linux-sles15-x86_64
module load nghttp2/1.57.0-clang-17.0.0-rocm5.7.1-mixed-rrynpn7
# ca-certificates-mozilla@=2023-05-30%clang@=17.0.0-rocm5.7.1-mixed build_system=generic arch=linux-sles15-x86_64
module load ca-certificates-mozilla/2023-05-30-clang-17.0.0-rocm5.7.1-mixed-lhn5rvo
# perl@=5.34.0%clang@=17.0.0-rocm5.7.1-mixed+cpanm+opcode+open+shared+threads build_system=generic arch=linux-sles15-x86_64
module load perl/5.34.0-clang-17.0.0-rocm5.7.1-mixed-7yifv3t
# zlib-ng@=2.1.6%clang@=17.0.0-rocm5.7.1-mixed+compat+new_strategies+opt build_system=autotools arch=linux-sles15-x86_64
module load zlib-ng/2.1.6-clang-17.0.0-rocm5.7.1-mixed-menugrv
# openssl@=3.2.1%clang@=17.0.0-rocm5.7.1-mixed~docs+shared build_system=generic certs=mozilla arch=linux-sles15-x86_64
## module load openssl/3.2.1-clang-17.0.0-rocm5.7.1-mixed-3ha6f2t
# curl@=8.6.0%clang@=17.0.0-rocm5.7.1-mixed~gssapi~ldap~libidn2~librtmp~libssh~libssh2+nghttp2 build_system=autotools libs=shared,static tls=openssl arch=linux-sles15-x86_64
module load curl/8.6.0-clang-17.0.0-rocm5.7.1-mixed-cq4gjom
# ncurses@=6.4%clang@=17.0.0-rocm5.7.1-mixed~symlinks+termlib abi=none build_system=autotools arch=linux-sles15-x86_64
module load ncurses/6.4-clang-17.0.0-rocm5.7.1-mixed-5fu3rph
# cmake@=3.20.6%clang@=17.0.0-rocm5.7.1-mixed~doc+ncurses+ownlibs build_system=generic build_type=Release arch=linux-sles15-x86_64
module load cmake/3.20.6-clang-17.0.0-rocm5.7.1-mixed-tsnlb3w
# blt@=0.4.1%clang@=17.0.0-rocm5.7.1-mixed build_system=generic arch=linux-sles15-x86_64
module load blt/0.4.1-clang-17.0.0-rocm5.7.1-mixed-yv7w5ka
# hip@=5.7.1%clang@=17.0.0-rocm5.7.1-mixed~cuda+rocm build_system=cmake build_type=Release generator=make patches=5bb9b0e,7668b2a,aee7249,b589a02,c2ee21c arch=linux-sles15-x86_64
module load hip/5.7.1-clang-17.0.0-rocm5.7.1-mixed-aokghsk
# hsa-rocr-dev@=5.7.1%clang@=17.0.0-rocm5.7.1-mixed~asan+image+shared build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load hsa-rocr-dev/5.7.1-clang-17.0.0-rocm5.7.1-mixed-4cjwtst
# llvm-amdgpu@=5.7.1%clang@=17.0.0-rocm5.7.1-mixed~link_llvm_dylib~llvm_dylib+rocm-device-libs build_system=cmake build_type=Release generator=ninja patches=53f9500,9a97712,b66529f arch=linux-sles15-x86_64
module load llvm-amdgpu/5.7.1-clang-17.0.0-rocm5.7.1-mixed-4xxmwvj
# camp@=0.2.3%clang@=17.0.0-rocm5.7.1-mixed~cuda~ipo~openmp+rocm~tests amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make patches=cb9e25b arch=linux-sles15-x86_64
module load camp/0.2.3-clang-17.0.0-rocm5.7.1-mixed-s4wguwz
# cray-mpich@=8.1.28%clang@=17.0.0-rocm5.7.1-mixed+wrappers build_system=generic arch=linux-sles15-x86_64
module load cray-mpich/8.1.28-clang-17.0.0-rocm5.7.1-mixed-i3f7um4
# gcc-runtime@=12.3-mixed%gcc@=12.3-mixed build_system=generic arch=linux-sles15-x86_64
module load gcc-runtime/12.3-mixed-gcc-12.3-mixed-mu6frky
## gmake@=4.3%gcc@=12.3-mixed~guile build_system=generic patches=599f134 arch=linux-sles15-x86_64
#module load gmake/4.3-gcc-12.3-mixed-xell4pg
## perl@=5.34.0%gcc@=12.3-mixed+cpanm+opcode+open+shared+threads build_system=generic arch=linux-sles15-x86_64
#module load perl/5.34.0-gcc-12.3-mixed-ukbfqpz
# openblas@=0.3.20%gcc@=12.3-mixed~bignuma~consistent_fpcsr~ilp64+locking+pic+shared build_system=makefile patches=9f12903 symbol_suffix=none threads=none arch=linux-sles15-x86_64
module load openblas/0.3.20-gcc-12.3-mixed-7iowws6
# coinhsl@=2019.05.21%gcc@=12.3-mixed+blas build_system=autotools arch=linux-sles15-x86_64
module load coinhsl/2019.05.21-gcc-12.3-mixed-4re5qz3
# hipblas@=5.7.1%clang@=17.0.0-rocm5.7.1-mixed~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make patches=f1b0687 arch=linux-sles15-x86_64
module load hipblas/5.7.1-clang-17.0.0-rocm5.7.1-mixed-majtcms
# hipsparse@=5.7.1%clang@=17.0.0-rocm5.7.1-mixed~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load hipsparse/5.7.1-clang-17.0.0-rocm5.7.1-mixed-mr2uw2o
# magma@=2.7.2%clang@=17.0.0-rocm5.7.1-mixed~cuda+fortran~ipo+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load magma/2.7.2-clang-17.0.0-rocm5.7.1-mixed-d6z2a7t
# metis@=5.1.0%clang@=17.0.0-rocm5.7.1-mixed~gdb~int64~ipo~real64+shared build_system=cmake build_type=Release generator=make patches=4991da9,93a7903 arch=linux-sles15-x86_64
module load metis/5.1.0-clang-17.0.0-rocm5.7.1-mixed-466rdfs
# rocprim@=5.7.1%clang@=17.0.0-rocm5.7.1-mixed amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load rocprim/5.7.1-clang-17.0.0-rocm5.7.1-mixed-eigcftr
# raja@=0.14.0%clang@=17.0.0-rocm5.7.1-mixed~cuda~examples~exercises~ipo~openmp~plugins+rocm+shared~tests amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load raja/0.14.0-clang-17.0.0-rocm5.7.1-mixed-6yzko55
# libiconv@=1.17%clang@=17.0.0-rocm5.7.1-mixed build_system=autotools libs=shared,static arch=linux-sles15-x86_64
module load libiconv/1.17-clang-17.0.0-rocm5.7.1-mixed-xa3wi7t
# diffutils@=3.10%clang@=17.0.0-rocm5.7.1-mixed build_system=autotools arch=linux-sles15-x86_64
module load diffutils/3.10-clang-17.0.0-rocm5.7.1-mixed-7v327yj
# libsigsegv@=2.14%clang@=17.0.0-rocm5.7.1-mixed build_system=autotools arch=linux-sles15-x86_64
module load libsigsegv/2.14-clang-17.0.0-rocm5.7.1-mixed-qa6g3ep
# m4@=1.4.19%clang@=17.0.0-rocm5.7.1-mixed+sigsegv build_system=autotools patches=9dc5fbd,bfdffa7 arch=linux-sles15-x86_64
module load m4/1.4.19-clang-17.0.0-rocm5.7.1-mixed-mmkqsk5
# autoconf@=2.72%clang@=17.0.0-rocm5.7.1-mixed build_system=autotools arch=linux-sles15-x86_64
module load autoconf/2.72-clang-17.0.0-rocm5.7.1-mixed-3xud5lx
# automake@=1.16.5%clang@=17.0.0-rocm5.7.1-mixed build_system=autotools arch=linux-sles15-x86_64
module load automake/1.16.5-clang-17.0.0-rocm5.7.1-mixed-fm7gg6r
# findutils@=4.9.0%clang@=17.0.0-rocm5.7.1-mixed build_system=autotools patches=440b954 arch=linux-sles15-x86_64
module load findutils/4.9.0-clang-17.0.0-rocm5.7.1-mixed-ks25vfq
# libtool@=2.4.7%clang@=17.0.0-rocm5.7.1-mixed build_system=autotools arch=linux-sles15-x86_64
module load libtool/2.4.7-clang-17.0.0-rocm5.7.1-mixed-hdagdfy
# gmp@=6.2.1%clang@=17.0.0-rocm5.7.1-mixed+cxx build_system=autotools libs=shared,static patches=69ad2e2 arch=linux-sles15-x86_64
module load gmp/6.2.1-clang-17.0.0-rocm5.7.1-mixed-ypvpkcj
# autoconf-archive@=2023.02.20%clang@=17.0.0-rocm5.7.1-mixed build_system=autotools arch=linux-sles15-x86_64
module load autoconf-archive/2023.02.20-clang-17.0.0-rocm5.7.1-mixed-badlwxr
# bzip2@=1.0.8%clang@=17.0.0-rocm5.7.1-mixed~debug~pic+shared build_system=generic arch=linux-sles15-x86_64
module load bzip2/1.0.8-clang-17.0.0-rocm5.7.1-mixed-fcspibr
# xz@=5.4.6%clang@=17.0.0-rocm5.7.1-mixed~pic build_system=autotools libs=shared,static arch=linux-sles15-x86_64
module load xz/5.4.6-clang-17.0.0-rocm5.7.1-mixed-lrxwsek
# libxml2@=2.10.3%clang@=17.0.0-rocm5.7.1-mixed+pic~python+shared build_system=autotools arch=linux-sles15-x86_64
module load libxml2/2.10.3-clang-17.0.0-rocm5.7.1-mixed-e6fl67c
# pigz@=2.8%clang@=17.0.0-rocm5.7.1-mixed build_system=makefile arch=linux-sles15-x86_64
module load pigz/2.8-clang-17.0.0-rocm5.7.1-mixed-gzgnvtm
# zstd@=1.5.5%clang@=17.0.0-rocm5.7.1-mixed+programs build_system=makefile compression=none libs=shared,static arch=linux-sles15-x86_64
module load zstd/1.5.5-clang-17.0.0-rocm5.7.1-mixed-cvopzgs
# tar@=1.34%clang@=17.0.0-rocm5.7.1-mixed build_system=autotools zip=pigz arch=linux-sles15-x86_64
module load tar/1.34-clang-17.0.0-rocm5.7.1-mixed-lzhjsn2
# gettext@=0.22.4%clang@=17.0.0-rocm5.7.1-mixed+bzip2+curses+git~libunistring+libxml2+pic+shared+tar+xz build_system=autotools arch=linux-sles15-x86_64
module load gettext/0.22.4-clang-17.0.0-rocm5.7.1-mixed-mmpz2s4
# texinfo@=7.0.3%clang@=17.0.0-rocm5.7.1-mixed build_system=autotools arch=linux-sles15-x86_64
module load texinfo/7.0.3-clang-17.0.0-rocm5.7.1-mixed-gcudjaq
# mpfr@=4.2.1%clang@=17.0.0-rocm5.7.1-mixed build_system=autotools libs=shared,static arch=linux-sles15-x86_64
module load mpfr/4.2.1-clang-17.0.0-rocm5.7.1-mixed-5kolmk3
# suite-sparse@=5.13.0%clang@=17.0.0-rocm5.7.1-mixed~cuda~graphblas~openmp+pic build_system=generic arch=linux-sles15-x86_64
module load suite-sparse/5.13.0-clang-17.0.0-rocm5.7.1-mixed-jk4j4xm
# umpire@=6.0.0%clang@=17.0.0-rocm5.7.1-mixed+c~cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make tests=none arch=linux-sles15-x86_64
module load umpire/6.0.0-clang-17.0.0-rocm5.7.1-mixed-dfwqrzx
# hiop@=develop%clang@=17.0.0-rocm5.7.1-mixed~cuda~deepchecking~ginkgo~ipo~jsrun+kron+mpi+raja+rocm~shared+sparse amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load hiop/develop-clang-17.0.0-rocm5.7.1-mixed-5rqwq45
# ipopt@=3.12.10%clang@=17.0.0-rocm5.7.1-mixed+coinhsl~debug~metis~mumps build_system=autotools arch=linux-sles15-x86_64
module load ipopt/3.12.10-clang-17.0.0-rocm5.7.1-mixed-2zdjszo
# python@=3.11.5%clang@=17.0.0-rocm5.7.1-mixed+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=13fa8bf,b0615b2,ebdca64,f2fd060 arch=linux-sles15-x86_64
module load python/3.11.5-clang-17.0.0-rocm5.7.1-mixed-ubwwl7x
# petsc@=3.20.4%clang@=17.0.0-rocm5.7.1-mixed~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~sycl~tetgen~trilinos~valgrind~zoltan build_system=generic clanguage=C memalign=none arch=linux-sles15-x86_64
module load petsc/3.20.4-clang-17.0.0-rocm5.7.1-mixed-2btkvwy
# exago@=develop%clang@=17.0.0-rocm5.7.1-mixed~cuda+hiop~ipo+ipopt+logging+mpi~python+raja+rocm amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
## module load exago/develop-clang-17.0.0-rocm5.7.1-mixed-6l3vqyp
