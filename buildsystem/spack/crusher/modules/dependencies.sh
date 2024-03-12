module use -a /lustre/orion/eng145/world-shared/spack-install/modules/linux-sles15-x86_64
# gmake@=4.4.1%clang@=16.0.0-rocm5.6.0-mixed~guile build_system=generic arch=linux-sles15-x86_64
module load gmake/4.4.1-clang-16.0.0-rocm5.6.0-mixed-k6al2qo
# pkgconf@=1.9.5%clang@=16.0.0-rocm5.6.0-mixed build_system=autotools arch=linux-sles15-x86_64
module load pkgconf/1.9.5-clang-16.0.0-rocm5.6.0-mixed-uouwfb7
# nghttp2@=1.57.0%clang@=16.0.0-rocm5.6.0-mixed build_system=autotools arch=linux-sles15-x86_64
module load nghttp2/1.57.0-clang-16.0.0-rocm5.6.0-mixed-4mp6nud
# ca-certificates-mozilla@=2023-05-30%clang@=16.0.0-rocm5.6.0-mixed build_system=generic arch=linux-sles15-x86_64
module load ca-certificates-mozilla/2023-05-30-clang-16.0.0-rocm5.6.0-mixed-b35brx3
# perl@=5.34.0%clang@=16.0.0-rocm5.6.0-mixed+cpanm+opcode+open+shared+threads build_system=generic arch=linux-sles15-x86_64
module load perl/5.34.0-clang-16.0.0-rocm5.6.0-mixed-kz4eken
# zlib-ng@=2.1.4%clang@=16.0.0-rocm5.6.0-mixed+compat+opt build_system=autotools arch=linux-sles15-x86_64
module load zlib-ng/2.1.4-clang-16.0.0-rocm5.6.0-mixed-3qkqsle
# openssl@=3.1.3%clang@=16.0.0-rocm5.6.0-mixed~docs+shared build_system=generic certs=mozilla arch=linux-sles15-x86_64
## module load openssl/3.1.3-clang-16.0.0-rocm5.6.0-mixed-z2b2n2v
# curl@=8.4.0%clang@=16.0.0-rocm5.6.0-mixed~gssapi~ldap~libidn2~librtmp~libssh~libssh2+nghttp2 build_system=autotools libs=shared,static tls=openssl arch=linux-sles15-x86_64
module load curl/8.4.0-clang-16.0.0-rocm5.6.0-mixed-f6bnt2g
# ncurses@=6.4%clang@=16.0.0-rocm5.6.0-mixed~symlinks+termlib abi=none build_system=autotools arch=linux-sles15-x86_64
module load ncurses/6.4-clang-16.0.0-rocm5.6.0-mixed-epb33eo
# cmake@=3.20.6%clang@=16.0.0-rocm5.6.0-mixed~doc+ncurses+ownlibs build_system=generic build_type=Release arch=linux-sles15-x86_64
module load cmake/3.20.6-clang-16.0.0-rocm5.6.0-mixed-yaow2oz
# blt@=0.4.1%clang@=16.0.0-rocm5.6.0-mixed build_system=generic arch=linux-sles15-x86_64
module load blt/0.4.1-clang-16.0.0-rocm5.6.0-mixed-7cfyib3
# hip@=5.6.0%clang@=16.0.0-rocm5.6.0-mixed~cuda+rocm build_system=cmake build_type=Release generator=make patches=aee7249,c2ee21c,e73e91b arch=linux-sles15-x86_64
module load hip/5.6.0-clang-16.0.0-rocm5.6.0-mixed-op3qkhw
# hsa-rocr-dev@=5.6.0%clang@=16.0.0-rocm5.6.0-mixed+image+shared build_system=cmake build_type=Release generator=make patches=9267179 arch=linux-sles15-x86_64
module load hsa-rocr-dev/5.6.0-clang-16.0.0-rocm5.6.0-mixed-ezukiit
# llvm-amdgpu@=5.6.0%clang@=16.0.0-rocm5.6.0-mixed~link_llvm_dylib~llvm_dylib~openmp+rocm-device-libs build_system=cmake build_type=Release generator=ninja patches=a08bbe1,b66529f,d35aec9 arch=linux-sles15-x86_64
module load llvm-amdgpu/5.6.0-clang-16.0.0-rocm5.6.0-mixed-hywg6dz
# camp@=0.2.3%clang@=16.0.0-rocm5.6.0-mixed~cuda~ipo~openmp+rocm~tests amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load camp/0.2.3-clang-16.0.0-rocm5.6.0-mixed-bpoapxo
# cray-mpich@=8.1.25%clang@=16.0.0-rocm5.6.0-mixed+wrappers build_system=generic arch=linux-sles15-x86_64
module load cray-mpich/8.1.25-clang-16.0.0-rocm5.6.0-mixed-7mbfpyl
# gmake@=4.4.1%gcc@=12.2.0-mixed~guile build_system=generic arch=linux-sles15-x86_64
module load gmake/4.4.1-gcc-12.2.0-mixed-j4cmkzq
# openblas@=0.3.20%gcc@=12.2.0-mixed~bignuma~consistent_fpcsr~ilp64+locking+pic+shared build_system=makefile patches=9f12903 symbol_suffix=none threads=none arch=linux-sles15-x86_64
module load openblas/0.3.20-gcc-12.2.0-mixed-tq6djj3
# coinhsl@=2019.05.21%gcc@=12.2.0-mixed+blas build_system=autotools arch=linux-sles15-x86_64
module load coinhsl/2019.05.21-gcc-12.2.0-mixed-tafgcu3
# hipblas@=5.6.0%clang@=16.0.0-rocm5.6.0-mixed~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load hipblas/5.6.0-clang-16.0.0-rocm5.6.0-mixed-lckx2a7
# hipsparse@=5.6.0%clang@=16.0.0-rocm5.6.0-mixed~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load hipsparse/5.6.0-clang-16.0.0-rocm5.6.0-mixed-v2eembx
# magma@=2.7.2%clang@=16.0.0-rocm5.6.0-mixed~cuda+fortran~ipo+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load magma/2.7.2-clang-16.0.0-rocm5.6.0-mixed-yfghydp
# metis@=5.1.0%clang@=16.0.0-rocm5.6.0-mixed~gdb~int64~ipo~real64+shared build_system=cmake build_type=Release generator=make patches=4991da9,93a7903 arch=linux-sles15-x86_64
module load metis/5.1.0-clang-16.0.0-rocm5.6.0-mixed-odzstni
# rocprim@=5.6.0%clang@=16.0.0-rocm5.6.0-mixed amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load rocprim/5.6.0-clang-16.0.0-rocm5.6.0-mixed-mbf63yj
# raja@=0.14.0%clang@=16.0.0-rocm5.6.0-mixed~cuda~examples~exercises~ipo~openmp~plugins+rocm+shared~tests amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load raja/0.14.0-clang-16.0.0-rocm5.6.0-mixed-z6v3hxh
# libiconv@=1.17%clang@=16.0.0-rocm5.6.0-mixed build_system=autotools libs=shared,static arch=linux-sles15-x86_64
module load libiconv/1.17-clang-16.0.0-rocm5.6.0-mixed-2aglobu
# diffutils@=3.9%clang@=16.0.0-rocm5.6.0-mixed build_system=autotools arch=linux-sles15-x86_64
module load diffutils/3.9-clang-16.0.0-rocm5.6.0-mixed-b7kmhwy
# libsigsegv@=2.14%clang@=16.0.0-rocm5.6.0-mixed build_system=autotools arch=linux-sles15-x86_64
module load libsigsegv/2.14-clang-16.0.0-rocm5.6.0-mixed-bdwxxrb
# m4@=1.4.19%clang@=16.0.0-rocm5.6.0-mixed+sigsegv build_system=autotools patches=9dc5fbd,bfdffa7 arch=linux-sles15-x86_64
module load m4/1.4.19-clang-16.0.0-rocm5.6.0-mixed-yehligm
# autoconf@=2.69%clang@=16.0.0-rocm5.6.0-mixed build_system=autotools patches=35c4492,7793209,a49dd5b arch=linux-sles15-x86_64
module load autoconf/2.69-clang-16.0.0-rocm5.6.0-mixed-dnxxoh2
# automake@=1.16.5%clang@=16.0.0-rocm5.6.0-mixed build_system=autotools arch=linux-sles15-x86_64
module load automake/1.16.5-clang-16.0.0-rocm5.6.0-mixed-3mm4fm5
# libtool@=2.4.7%clang@=16.0.0-rocm5.6.0-mixed build_system=autotools arch=linux-sles15-x86_64
module load libtool/2.4.7-clang-16.0.0-rocm5.6.0-mixed-23ehtda
# gmp@=6.2.1%clang@=16.0.0-rocm5.6.0-mixed+cxx build_system=autotools libs=shared,static patches=69ad2e2 arch=linux-sles15-x86_64
module load gmp/6.2.1-clang-16.0.0-rocm5.6.0-mixed-t5hlxtw
# autoconf-archive@=2023.02.20%clang@=16.0.0-rocm5.6.0-mixed build_system=autotools arch=linux-sles15-x86_64
module load autoconf-archive/2023.02.20-clang-16.0.0-rocm5.6.0-mixed-exate7t
# bzip2@=1.0.8%clang@=16.0.0-rocm5.6.0-mixed~debug~pic+shared build_system=generic arch=linux-sles15-x86_64
module load bzip2/1.0.8-clang-16.0.0-rocm5.6.0-mixed-jcvtab5
# xz@=5.4.1%clang@=16.0.0-rocm5.6.0-mixed~pic build_system=autotools libs=shared,static arch=linux-sles15-x86_64
module load xz/5.4.1-clang-16.0.0-rocm5.6.0-mixed-jny6a2y
# libxml2@=2.10.3%clang@=16.0.0-rocm5.6.0-mixed+pic~python+shared build_system=autotools arch=linux-sles15-x86_64
module load libxml2/2.10.3-clang-16.0.0-rocm5.6.0-mixed-bwzce2v
# pigz@=2.7%clang@=16.0.0-rocm5.6.0-mixed build_system=makefile arch=linux-sles15-x86_64
module load pigz/2.7-clang-16.0.0-rocm5.6.0-mixed-hndhspi
# zstd@=1.5.5%clang@=16.0.0-rocm5.6.0-mixed+programs build_system=makefile compression=none libs=shared,static arch=linux-sles15-x86_64
module load zstd/1.5.5-clang-16.0.0-rocm5.6.0-mixed-w7bcfca
# tar@=1.34%clang@=16.0.0-rocm5.6.0-mixed build_system=autotools zip=pigz arch=linux-sles15-x86_64
module load tar/1.34-clang-16.0.0-rocm5.6.0-mixed-yvw47t5
# gettext@=0.22.3%clang@=16.0.0-rocm5.6.0-mixed+bzip2+curses+git~libunistring+libxml2+pic+shared+tar+xz build_system=autotools arch=linux-sles15-x86_64
module load gettext/0.22.3-clang-16.0.0-rocm5.6.0-mixed-iliwbcy
# texinfo@=7.0.3%clang@=16.0.0-rocm5.6.0-mixed build_system=autotools arch=linux-sles15-x86_64
module load texinfo/7.0.3-clang-16.0.0-rocm5.6.0-mixed-gkcyjes
# mpfr@=4.2.0%clang@=16.0.0-rocm5.6.0-mixed build_system=autotools libs=shared,static arch=linux-sles15-x86_64
module load mpfr/4.2.0-clang-16.0.0-rocm5.6.0-mixed-j5g7ch2
# suite-sparse@=5.13.0%clang@=16.0.0-rocm5.6.0-mixed~cuda~graphblas~openmp+pic build_system=generic arch=linux-sles15-x86_64
module load suite-sparse/5.13.0-clang-16.0.0-rocm5.6.0-mixed-riysyew
# umpire@=6.0.0%clang@=16.0.0-rocm5.6.0-mixed+c~cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make tests=none arch=linux-sles15-x86_64
module load umpire/6.0.0-clang-16.0.0-rocm5.6.0-mixed-2li2g3r
# hiop@=develop%clang@=16.0.0-rocm5.6.0-mixed~cuda~deepchecking~ginkgo~ipo~jsrun+kron+mpi+raja+rocm~shared+sparse amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load hiop/develop-clang-16.0.0-rocm5.6.0-mixed-3ozi52g
# ipopt@=3.12.10%clang@=16.0.0-rocm5.6.0-mixed+coinhsl~debug~metis~mumps build_system=autotools arch=linux-sles15-x86_64
module load ipopt/3.12.10-clang-16.0.0-rocm5.6.0-mixed-rlgcopl
# python@=3.9.12%clang@=16.0.0-rocm5.6.0-mixed+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=0d98e93,4c24573,ebdca64,f2fd060 arch=linux-sles15-x86_64
module load python/3.9.12-clang-16.0.0-rocm5.6.0-mixed-mdt6uvm
# petsc@=3.20.1%clang@=16.0.0-rocm5.6.0-mixed~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~sycl~tetgen~trilinos~valgrind build_system=generic clanguage=C memalign=none arch=linux-sles15-x86_64
module load petsc/3.20.1-clang-16.0.0-rocm5.6.0-mixed-pa2upkm
# exago@=develop%clang@=16.0.0-rocm5.6.0-mixed~cuda+hiop~ipo+ipopt+logging+mpi~python+raja+rocm amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
## module load exago/develop-clang-16.0.0-rocm5.6.0-mixed-6xcrrvo
