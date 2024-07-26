module use -a /lustre/orion/eng145/world-shared/spack-install/modules/linux-sles15-x86_64
# gmake@=4.3%clang@=17.0.0-rocm5.7.1~guile build_system=generic patches=599f134 arch=linux-sles15-x86_64
module load gmake/4.3-clang-17.0.0-rocm5.7.1-nddwmha
# pkgconf@=1.9.5%clang@=17.0.0-rocm5.7.1 build_system=autotools arch=linux-sles15-x86_64
module load pkgconf/1.9.5-clang-17.0.0-rocm5.7.1-nwo2jra
# nghttp2@=1.57.0%clang@=17.0.0-rocm5.7.1 build_system=autotools arch=linux-sles15-x86_64
module load nghttp2/1.57.0-clang-17.0.0-rocm5.7.1-4lexdyk
# ca-certificates-mozilla@=2023-05-30%clang@=17.0.0-rocm5.7.1 build_system=generic arch=linux-sles15-x86_64
module load ca-certificates-mozilla/2023-05-30-clang-17.0.0-rocm5.7.1-ndecgeb
# perl@=5.34.0%clang@=17.0.0-rocm5.7.1+cpanm+opcode+open+shared+threads build_system=generic arch=linux-sles15-x86_64
module load perl/5.34.0-clang-17.0.0-rocm5.7.1-5wz23qm
# zlib-ng@=2.1.6%clang@=17.0.0-rocm5.7.1+compat+new_strategies+opt build_system=autotools arch=linux-sles15-x86_64
module load zlib-ng/2.1.6-clang-17.0.0-rocm5.7.1-ram6oz7
# openssl@=3.2.1%clang@=17.0.0-rocm5.7.1~docs+shared build_system=generic certs=mozilla arch=linux-sles15-x86_64
## module load openssl/3.2.1-clang-17.0.0-rocm5.7.1-pdstdgu
# curl@=8.6.0%clang@=17.0.0-rocm5.7.1~gssapi~ldap~libidn2~librtmp~libssh~libssh2+nghttp2 build_system=autotools libs=shared,static tls=openssl arch=linux-sles15-x86_64
module load curl/8.6.0-clang-17.0.0-rocm5.7.1-avjuvwe
# ncurses@=6.4%clang@=17.0.0-rocm5.7.1~symlinks+termlib abi=none build_system=autotools arch=linux-sles15-x86_64
module load ncurses/6.4-clang-17.0.0-rocm5.7.1-bkyrvxv
# cmake@=3.20.6%clang@=17.0.0-rocm5.7.1~doc+ncurses+ownlibs build_system=generic build_type=Release arch=linux-sles15-x86_64
module load cmake/3.20.6-clang-17.0.0-rocm5.7.1-3r6do3u
# blt@=0.4.1%clang@=17.0.0-rocm5.7.1 build_system=generic arch=linux-sles15-x86_64
module load blt/0.4.1-clang-17.0.0-rocm5.7.1-ob4xxpn
# hip@=5.7.1%clang@=17.0.0-rocm5.7.1~cuda+rocm build_system=cmake build_type=Release generator=make patches=5bb9b0e,7668b2a,aee7249,b589a02,c2ee21c arch=linux-sles15-x86_64
module load hip/5.7.1-clang-17.0.0-rocm5.7.1-75dn5yy
# hsa-rocr-dev@=5.7.1%clang@=17.0.0-rocm5.7.1~asan+image+shared build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load hsa-rocr-dev/5.7.1-clang-17.0.0-rocm5.7.1-mizwkrm
# llvm-amdgpu@=5.7.1%clang@=17.0.0-rocm5.7.1~link_llvm_dylib~llvm_dylib+rocm-device-libs build_system=cmake build_type=Release generator=ninja patches=53f9500,9a97712,b66529f arch=linux-sles15-x86_64
module load llvm-amdgpu/5.7.1-clang-17.0.0-rocm5.7.1-mi32544
# camp@=0.2.3%clang@=17.0.0-rocm5.7.1~cuda~ipo~openmp+rocm~tests amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make patches=cb9e25b arch=linux-sles15-x86_64
module load camp/0.2.3-clang-17.0.0-rocm5.7.1-c55efcq
# cray-mpich@=8.1.28%clang@=17.0.0-rocm5.7.1+wrappers build_system=generic arch=linux-sles15-x86_64
module load cray-mpich/8.1.28-clang-17.0.0-rocm5.7.1-hlzs2oy
# openblas@=0.3.20%clang@=17.0.0-rocm5.7.1~bignuma~consistent_fpcsr~ilp64+locking+pic+shared build_system=makefile patches=9968625,9f12903 symbol_suffix=none threads=none arch=linux-sles15-x86_64
module load openblas/0.3.20-clang-17.0.0-rocm5.7.1-slikqnh
# coinhsl@=2019.05.21%clang@=17.0.0-rocm5.7.1+blas build_system=autotools arch=linux-sles15-x86_64
module load coinhsl/2019.05.21-clang-17.0.0-rocm5.7.1-wkelt3i
# hipblas@=5.7.1%clang@=17.0.0-rocm5.7.1~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make patches=f1b0687 arch=linux-sles15-x86_64
module load hipblas/5.7.1-clang-17.0.0-rocm5.7.1-7dopdxg
# hipsparse@=5.7.1%clang@=17.0.0-rocm5.7.1~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load hipsparse/5.7.1-clang-17.0.0-rocm5.7.1-dl3ghr6
# magma@=2.7.2%clang@=17.0.0-rocm5.7.1~cuda+fortran~ipo+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load magma/2.7.2-clang-17.0.0-rocm5.7.1-apopwgf
# metis@=5.1.0%clang@=17.0.0-rocm5.7.1~gdb~int64~ipo~real64+shared build_system=cmake build_type=Release generator=make patches=4991da9,93a7903 arch=linux-sles15-x86_64
module load metis/5.1.0-clang-17.0.0-rocm5.7.1-psr5tyz
# rocprim@=5.7.1%clang@=17.0.0-rocm5.7.1 amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load rocprim/5.7.1-clang-17.0.0-rocm5.7.1-gakrls3
# raja@=0.14.0%clang@=17.0.0-rocm5.7.1~cuda~examples~exercises~ipo~openmp~plugins+rocm+shared~tests amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load raja/0.14.0-clang-17.0.0-rocm5.7.1-onmgi2d
# libiconv@=1.17%clang@=17.0.0-rocm5.7.1 build_system=autotools libs=shared,static arch=linux-sles15-x86_64
module load libiconv/1.17-clang-17.0.0-rocm5.7.1-hnl6m24
# diffutils@=3.10%clang@=17.0.0-rocm5.7.1 build_system=autotools arch=linux-sles15-x86_64
module load diffutils/3.10-clang-17.0.0-rocm5.7.1-enb5vdz
# libsigsegv@=2.14%clang@=17.0.0-rocm5.7.1 build_system=autotools arch=linux-sles15-x86_64
module load libsigsegv/2.14-clang-17.0.0-rocm5.7.1-ruhgxlh
# m4@=1.4.19%clang@=17.0.0-rocm5.7.1+sigsegv build_system=autotools patches=9dc5fbd,bfdffa7 arch=linux-sles15-x86_64
module load m4/1.4.19-clang-17.0.0-rocm5.7.1-kt4rgm4
# autoconf@=2.72%clang@=17.0.0-rocm5.7.1 build_system=autotools arch=linux-sles15-x86_64
module load autoconf/2.72-clang-17.0.0-rocm5.7.1-3znqc65
# automake@=1.16.5%clang@=17.0.0-rocm5.7.1 build_system=autotools arch=linux-sles15-x86_64
module load automake/1.16.5-clang-17.0.0-rocm5.7.1-7zcyrz5
# findutils@=4.9.0%clang@=17.0.0-rocm5.7.1 build_system=autotools patches=440b954 arch=linux-sles15-x86_64
module load findutils/4.9.0-clang-17.0.0-rocm5.7.1-cqhq5yp
# libtool@=2.4.7%clang@=17.0.0-rocm5.7.1 build_system=autotools arch=linux-sles15-x86_64
module load libtool/2.4.7-clang-17.0.0-rocm5.7.1-dxuceh5
# gmp@=6.2.1%clang@=17.0.0-rocm5.7.1+cxx build_system=autotools libs=shared,static patches=69ad2e2 arch=linux-sles15-x86_64
module load gmp/6.2.1-clang-17.0.0-rocm5.7.1-fism3yk
# autoconf-archive@=2023.02.20%clang@=17.0.0-rocm5.7.1 build_system=autotools arch=linux-sles15-x86_64
module load autoconf-archive/2023.02.20-clang-17.0.0-rocm5.7.1-cu5ts4p
# bzip2@=1.0.8%clang@=17.0.0-rocm5.7.1~debug~pic+shared build_system=generic arch=linux-sles15-x86_64
module load bzip2/1.0.8-clang-17.0.0-rocm5.7.1-uwso2l7
# xz@=5.4.6%clang@=17.0.0-rocm5.7.1~pic build_system=autotools libs=shared,static arch=linux-sles15-x86_64
module load xz/5.4.6-clang-17.0.0-rocm5.7.1-767bzjj
# libxml2@=2.10.3%clang@=17.0.0-rocm5.7.1+pic~python+shared build_system=autotools arch=linux-sles15-x86_64
module load libxml2/2.10.3-clang-17.0.0-rocm5.7.1-gy576lq
# pigz@=2.8%clang@=17.0.0-rocm5.7.1 build_system=makefile arch=linux-sles15-x86_64
module load pigz/2.8-clang-17.0.0-rocm5.7.1-xmqkryq
# zstd@=1.5.5%clang@=17.0.0-rocm5.7.1+programs build_system=makefile compression=none libs=shared,static arch=linux-sles15-x86_64
module load zstd/1.5.5-clang-17.0.0-rocm5.7.1-pcjm6lw
# tar@=1.34%clang@=17.0.0-rocm5.7.1 build_system=autotools zip=pigz arch=linux-sles15-x86_64
module load tar/1.34-clang-17.0.0-rocm5.7.1-uojmtx3
# gettext@=0.22.4%clang@=17.0.0-rocm5.7.1+bzip2+curses+git~libunistring+libxml2+pic+shared+tar+xz build_system=autotools arch=linux-sles15-x86_64
module load gettext/0.22.4-clang-17.0.0-rocm5.7.1-oz3jmnq
# texinfo@=7.0.3%clang@=17.0.0-rocm5.7.1 build_system=autotools arch=linux-sles15-x86_64
module load texinfo/7.0.3-clang-17.0.0-rocm5.7.1-4j7pjhb
# mpfr@=4.2.1%clang@=17.0.0-rocm5.7.1 build_system=autotools libs=shared,static arch=linux-sles15-x86_64
module load mpfr/4.2.1-clang-17.0.0-rocm5.7.1-5jbmi3f
# suite-sparse@=5.13.0%clang@=17.0.0-rocm5.7.1~cuda~graphblas~openmp+pic build_system=generic arch=linux-sles15-x86_64
module load suite-sparse/5.13.0-clang-17.0.0-rocm5.7.1-cdejlpe
# umpire@=6.0.0%clang@=17.0.0-rocm5.7.1+c~cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make tests=none arch=linux-sles15-x86_64
module load umpire/6.0.0-clang-17.0.0-rocm5.7.1-s4nlsx5
# hiop@=develop%clang@=17.0.0-rocm5.7.1~cuda~deepchecking~ginkgo~ipo~jsrun+kron+mpi+raja+rocm~shared+sparse amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load hiop/develop-clang-17.0.0-rocm5.7.1-nmp3foh
# ipopt@=3.12.10%clang@=17.0.0-rocm5.7.1+coinhsl~debug~metis~mumps build_system=autotools arch=linux-sles15-x86_64
module load ipopt/3.12.10-clang-17.0.0-rocm5.7.1-bqpb4bi
# python@=3.11.5%clang@=17.0.0-rocm5.7.1+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=13fa8bf,b0615b2,ebdca64,f2fd060 arch=linux-sles15-x86_64
module load python/3.11.5-clang-17.0.0-rocm5.7.1-kmjxg2r
# petsc@=3.20.4%clang@=17.0.0-rocm5.7.1~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~sycl~tetgen~trilinos~valgrind~zoltan build_system=generic clanguage=C memalign=none arch=linux-sles15-x86_64
module load petsc/3.20.4-clang-17.0.0-rocm5.7.1-7b6bxwi
# exago@=develop%clang@=17.0.0-rocm5.7.1~cuda+hiop~ipo+ipopt+logging+mpi~python+raja+rocm amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
## module load exago/develop-clang-17.0.0-rocm5.7.1-vc5kuoe
