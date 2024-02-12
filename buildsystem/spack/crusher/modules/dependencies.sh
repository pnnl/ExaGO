module use -a /lustre/orion/eng145/proj-shared/nkouk/spack-install/modules/linux-sles15-zen3
# gmake@=4.4.1%clang@=14.0.0-rocm5.2.0-mixed~guile build_system=generic arch=linux-sles15-zen3
module load gmake/4.4.1-clang-14.0.0-rocm5.2.0-mixed-7vcngsz
# pkgconf@=1.9.5%clang@=14.0.0-rocm5.2.0-mixed build_system=autotools arch=linux-sles15-zen3
module load pkgconf/1.9.5-clang-14.0.0-rocm5.2.0-mixed-4swj7eh
# nghttp2@=1.57.0%clang@=14.0.0-rocm5.2.0-mixed build_system=autotools arch=linux-sles15-zen3
module load nghttp2/1.57.0-clang-14.0.0-rocm5.2.0-mixed-iwxpikh
# ca-certificates-mozilla@=2023-05-30%clang@=14.0.0-rocm5.2.0-mixed build_system=generic arch=linux-sles15-zen3
module load ca-certificates-mozilla/2023-05-30-clang-14.0.0-rocm5.2.0-mixed-q3hsdil
# perl@=5.34.0%clang@=14.0.0-rocm5.2.0-mixed+cpanm+opcode+open+shared+threads build_system=generic arch=linux-sles15-zen3
module load perl/5.34.0-clang-14.0.0-rocm5.2.0-mixed-ex2kxdz
# zlib-ng@=2.1.5%clang@=14.0.0-rocm5.2.0-mixed+compat+opt build_system=autotools arch=linux-sles15-zen3
module load zlib-ng/2.1.5-clang-14.0.0-rocm5.2.0-mixed-nhxuhp3
# openssl@=3.1.3%clang@=14.0.0-rocm5.2.0-mixed~docs+shared build_system=generic certs=mozilla arch=linux-sles15-zen3
## module load openssl/3.1.3-clang-14.0.0-rocm5.2.0-mixed-y5566cf
# curl@=8.4.0%clang@=14.0.0-rocm5.2.0-mixed~gssapi~ldap~libidn2~librtmp~libssh~libssh2+nghttp2 build_system=autotools libs=shared,static tls=openssl arch=linux-sles15-zen3
module load curl/8.4.0-clang-14.0.0-rocm5.2.0-mixed-2yz3o3v
# ncurses@=6.4%clang@=14.0.0-rocm5.2.0-mixed~symlinks+termlib abi=none build_system=autotools arch=linux-sles15-zen3
module load ncurses/6.4-clang-14.0.0-rocm5.2.0-mixed-vy6yuss
# cmake@=3.20.6%clang@=14.0.0-rocm5.2.0-mixed~doc+ncurses+ownlibs build_system=generic build_type=Release arch=linux-sles15-zen3
module load cmake/3.20.6-clang-14.0.0-rocm5.2.0-mixed-lilujjv
# blt@=0.4.1%clang@=14.0.0-rocm5.2.0-mixed build_system=generic arch=linux-sles15-zen3
module load blt/0.4.1-clang-14.0.0-rocm5.2.0-mixed-pihruyq
# hip@=5.2.0%clang@=14.0.0-rocm5.2.0-mixed~cuda+rocm build_system=cmake build_type=Release generator=make patches=959d1fe,c2ee21c arch=linux-sles15-zen3
module load hip/5.2.0-clang-14.0.0-rocm5.2.0-mixed-6ftqihk
# hsa-rocr-dev@=5.2.0%clang@=14.0.0-rocm5.2.0-mixed~asan+image+shared build_system=cmake build_type=Release generator=make patches=9267179 arch=linux-sles15-zen3
module load hsa-rocr-dev/5.2.0-clang-14.0.0-rocm5.2.0-mixed-p6fvpym
# llvm-amdgpu@=5.2.0%clang@=14.0.0-rocm5.2.0-mixed~link_llvm_dylib~llvm_dylib~openmp+rocm-device-libs build_system=cmake build_type=Release generator=ninja patches=a08bbe1 arch=linux-sles15-zen3
module load llvm-amdgpu/5.2.0-clang-14.0.0-rocm5.2.0-mixed-7g3i4ov
# camp@=0.2.3%clang@=14.0.0-rocm5.2.0-mixed~cuda~ipo~openmp+rocm~tests amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load camp/0.2.3-clang-14.0.0-rocm5.2.0-mixed-5hqx2ao
# cray-mpich@=8.1.25%clang@=14.0.0-rocm5.2.0-mixed+wrappers build_system=generic arch=linux-sles15-zen3
module load cray-mpich/8.1.25-clang-14.0.0-rocm5.2.0-mixed-qd6s47h
# hipblas@=5.2.0%clang@=14.0.0-rocm5.2.0-mixed~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load hipblas/5.2.0-clang-14.0.0-rocm5.2.0-mixed-pz2pdk4
# hipfft@=5.2.0%clang@=14.0.0-rocm5.2.0-mixed~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load hipfft/5.2.0-clang-14.0.0-rocm5.2.0-mixed-pje75d4
# hiprand@=5.2.0%clang@=14.0.0-rocm5.2.0-mixed~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load hiprand/5.2.0-clang-14.0.0-rocm5.2.0-mixed-224knmx
# hipsparse@=5.2.0%clang@=14.0.0-rocm5.2.0-mixed~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make patches=c447537 arch=linux-sles15-zen3
module load hipsparse/5.2.0-clang-14.0.0-rocm5.2.0-mixed-vqshjsa
# rocprim@=5.2.0%clang@=14.0.0-rocm5.2.0-mixed amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load rocprim/5.2.0-clang-14.0.0-rocm5.2.0-mixed-6jqmhkg
# rocrand@=5.2.0%clang@=14.0.0-rocm5.2.0-mixed+hiprand amdgpu_target=auto build_system=cmake build_type=Release generator=make patches=a35e689 arch=linux-sles15-zen3
module load rocrand/5.2.0-clang-14.0.0-rocm5.2.0-mixed-nwexus7
# rocthrust@=5.2.0%clang@=14.0.0-rocm5.2.0-mixed amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load rocthrust/5.2.0-clang-14.0.0-rocm5.2.0-mixed-3upg7wr
# ginkgo@=1.5.0.glu_experimental%clang@=14.0.0-rocm5.2.0-mixed~cuda~develtools~full_optimizations~hwloc~ipo~mpi~openmp+rocm+shared~sycl amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make patches=ba0956e arch=linux-sles15-zen3
module load ginkgo/1.5.0.glu_experimental-clang-14.0.0-rocm5.2.0-mixed-6bembpq
# gmake@=4.4.1%gcc@=12.2.0-mixed~guile build_system=generic arch=linux-sles15-zen3
module load gmake/4.4.1-gcc-12.2.0-mixed-svg2oxk
# openblas@=0.3.20%gcc@=12.2.0-mixed~bignuma~consistent_fpcsr~ilp64+locking+pic+shared build_system=makefile patches=9f12903 symbol_suffix=none threads=none arch=linux-sles15-zen3
module load openblas/0.3.20-gcc-12.2.0-mixed-lkwa4l7
# coinhsl@=2019.05.21%gcc@=12.2.0-mixed+blas build_system=autotools arch=linux-sles15-zen3
module load coinhsl/2019.05.21-gcc-12.2.0-mixed-4lf6df7
# magma@=2.6.2%clang@=14.0.0-rocm5.2.0-mixed~cuda+fortran~ipo+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load magma/2.6.2-clang-14.0.0-rocm5.2.0-mixed-5v4sbqv
# metis@=5.1.0%clang@=14.0.0-rocm5.2.0-mixed~gdb~int64~ipo~real64+shared build_system=cmake build_type=Release generator=make patches=4991da9,93a7903 arch=linux-sles15-zen3
module load metis/5.1.0-clang-14.0.0-rocm5.2.0-mixed-3loietj
# raja@=0.14.0%clang@=14.0.0-rocm5.2.0-mixed~cuda~examples~exercises~ipo~openmp~plugins+rocm+shared~tests amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load raja/0.14.0-clang-14.0.0-rocm5.2.0-mixed-plistg5
# libiconv@=1.17%clang@=14.0.0-rocm5.2.0-mixed build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load libiconv/1.17-clang-14.0.0-rocm5.2.0-mixed-gryio4w
# diffutils@=3.9%clang@=14.0.0-rocm5.2.0-mixed build_system=autotools arch=linux-sles15-zen3
module load diffutils/3.9-clang-14.0.0-rocm5.2.0-mixed-yz7kuqz
# libsigsegv@=2.14%clang@=14.0.0-rocm5.2.0-mixed build_system=autotools arch=linux-sles15-zen3
module load libsigsegv/2.14-clang-14.0.0-rocm5.2.0-mixed-mjraljk
# m4@=1.4.19%clang@=14.0.0-rocm5.2.0-mixed+sigsegv build_system=autotools patches=9dc5fbd,bfdffa7 arch=linux-sles15-zen3
module load m4/1.4.19-clang-14.0.0-rocm5.2.0-mixed-xrlual5
# autoconf@=2.69%clang@=14.0.0-rocm5.2.0-mixed build_system=autotools patches=35c4492,7793209,a49dd5b arch=linux-sles15-zen3
module load autoconf/2.69-clang-14.0.0-rocm5.2.0-mixed-4q2ohdh
# automake@=1.16.5%clang@=14.0.0-rocm5.2.0-mixed build_system=autotools arch=linux-sles15-zen3
module load automake/1.16.5-clang-14.0.0-rocm5.2.0-mixed-vpi6wmc
# libtool@=2.4.7%clang@=14.0.0-rocm5.2.0-mixed build_system=autotools arch=linux-sles15-zen3
module load libtool/2.4.7-clang-14.0.0-rocm5.2.0-mixed-lm5hq6p
# gmp@=6.2.1%clang@=14.0.0-rocm5.2.0-mixed+cxx build_system=autotools libs=shared,static patches=69ad2e2 arch=linux-sles15-zen3
module load gmp/6.2.1-clang-14.0.0-rocm5.2.0-mixed-kvk75vv
# autoconf-archive@=2023.02.20%clang@=14.0.0-rocm5.2.0-mixed build_system=autotools arch=linux-sles15-zen3
module load autoconf-archive/2023.02.20-clang-14.0.0-rocm5.2.0-mixed-2gml5wh
# bzip2@=1.0.8%clang@=14.0.0-rocm5.2.0-mixed~debug~pic+shared build_system=generic arch=linux-sles15-zen3
module load bzip2/1.0.8-clang-14.0.0-rocm5.2.0-mixed-e22jeoj
# xz@=5.4.1%clang@=14.0.0-rocm5.2.0-mixed~pic build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load xz/5.4.1-clang-14.0.0-rocm5.2.0-mixed-il6uern
# libxml2@=2.10.3%clang@=14.0.0-rocm5.2.0-mixed+pic~python+shared build_system=autotools arch=linux-sles15-zen3
module load libxml2/2.10.3-clang-14.0.0-rocm5.2.0-mixed-rs5hdf4
# pigz@=2.7%clang@=14.0.0-rocm5.2.0-mixed build_system=makefile arch=linux-sles15-zen3
module load pigz/2.7-clang-14.0.0-rocm5.2.0-mixed-ghuoaqc
# zstd@=1.5.5%clang@=14.0.0-rocm5.2.0-mixed+programs build_system=makefile compression=none libs=shared,static arch=linux-sles15-zen3
module load zstd/1.5.5-clang-14.0.0-rocm5.2.0-mixed-q7yzxzw
# tar@=1.34%clang@=14.0.0-rocm5.2.0-mixed build_system=autotools zip=pigz arch=linux-sles15-zen3
module load tar/1.34-clang-14.0.0-rocm5.2.0-mixed-xjye57q
# gettext@=0.22.4%clang@=14.0.0-rocm5.2.0-mixed+bzip2+curses+git~libunistring+libxml2+pic+shared+tar+xz build_system=autotools arch=linux-sles15-zen3
module load gettext/0.22.4-clang-14.0.0-rocm5.2.0-mixed-mdy7wdc
# texinfo@=7.0.3%clang@=14.0.0-rocm5.2.0-mixed build_system=autotools arch=linux-sles15-zen3
module load texinfo/7.0.3-clang-14.0.0-rocm5.2.0-mixed-zpog5zs
# mpfr@=4.2.0%clang@=14.0.0-rocm5.2.0-mixed build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load mpfr/4.2.0-clang-14.0.0-rocm5.2.0-mixed-k6y6ufe
# suite-sparse@=5.13.0%clang@=14.0.0-rocm5.2.0-mixed~cuda~graphblas~openmp+pic build_system=generic arch=linux-sles15-zen3
module load suite-sparse/5.13.0-clang-14.0.0-rocm5.2.0-mixed-vfmchqd
# umpire@=6.0.0%clang@=14.0.0-rocm5.2.0-mixed+c~cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make tests=none arch=linux-sles15-zen3
module load umpire/6.0.0-clang-14.0.0-rocm5.2.0-mixed-q3vler4
# hiop@=develop%clang@=14.0.0-rocm5.2.0-mixed~cuda~deepchecking+ginkgo~ipo~jsrun+kron+mpi+raja+rocm~shared+sparse amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load hiop/develop-clang-14.0.0-rocm5.2.0-mixed-jxs7pn7
# ipopt@=3.12.10%clang@=14.0.0-rocm5.2.0-mixed+coinhsl~debug~metis~mumps build_system=autotools arch=linux-sles15-zen3
module load ipopt/3.12.10-clang-14.0.0-rocm5.2.0-mixed-di6sy6h
# python@=3.9.12%clang@=14.0.0-rocm5.2.0-mixed+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=0d98e93,4c24573,ebdca64,f2fd060 arch=linux-sles15-zen3
module load python/3.9.12-clang-14.0.0-rocm5.2.0-mixed-rmdxnw5
# petsc@=3.20.2%clang@=14.0.0-rocm5.2.0-mixed~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~sycl~tetgen~trilinos~valgrind~zoltan build_system=generic clanguage=C memalign=none arch=linux-sles15-zen3
module load petsc/3.20.2-clang-14.0.0-rocm5.2.0-mixed-lnhejpn
# exago@=develop%clang@=14.0.0-rocm5.2.0-mixed~cuda+hiop~ipo+ipopt+logging+mpi~python+raja+rocm amdgpu_target=gfx90a build_system=cmake build_type=Release dev_path=/lustre/orion/scratch/nkouk/eng145/ExaGO generator=make arch=linux-sles15-zen3
## module load exago/develop-clang-14.0.0-rocm5.2.0-mixed-stynmti
