module use -a /lustre/orion/stf006/world-shared/nkouk/exago/spack-install/modules/linux-sles15-zen3
# cmake@=3.27.9%rocmcc@=6.3.1~doc+ncurses+ownlibs~qtgui build_system=generic build_type=Release patches=dbc3892 arch=linux-sles15-zen3
module load cmake/3.27.9-rocmcc-6.3.1-glbjtfr
# glibc@=2.31%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-zen3
module load glibc/2.31-rocmcc-6.3.1-gjx74zv
# blt@=0.4.1%rocmcc@=6.3.1 build_system=generic arch=linux-sles15-zen3
module load blt/0.4.1-rocmcc-6.3.1-patyhju
# gmake@=4.4.1%rocmcc@=6.3.1~guile build_system=generic arch=linux-sles15-zen3
module load gmake/4.4.1-rocmcc-6.3.1-ugara3u
# hip@=6.3.1%rocmcc@=6.3.1~asan~cuda+rocm build_system=cmake build_type=Release generator=make patches=1f65dfe arch=linux-sles15-zen3
module load hip/6.3.1-rocmcc-6.3.1-je2yq3g
# hsa-rocr-dev@=6.3.1%rocmcc@=6.3.1~asan+image+shared build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load hsa-rocr-dev/6.3.1-rocmcc-6.3.1-b55uik5
# llvm-amdgpu@=6.3.1%rocmcc@=6.3.1~link_llvm_dylib~llvm_dylib+rocm-device-libs build_system=cmake build_type=Release generator=ninja patches=b4774ca arch=linux-sles15-zen3
module load llvm-amdgpu/6.3.1-rocmcc-6.3.1-7ltelvs
# camp@=0.2.3%rocmcc@=6.3.1~cuda~ipo~omptarget~openmp+rocm~sycl~tests amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make patches=cb9e25b,f854571 arch=linux-sles15-zen3
module load camp/0.2.3-rocmcc-6.3.1-7uwcpvu
# cray-mpich@=8.1.28%rocmcc@=6.3.1+wrappers build_system=generic arch=linux-sles15-zen3
module load cray-mpich/8.1.28-rocmcc-6.3.1-ohlzfme
# gcc-runtime@=12.3%gcc@=12.3 build_system=generic arch=linux-sles15-zen3
module load gcc-runtime/12.3-gcc-12.3-tmqowzy
# gmake@=4.4.1%gcc@=12.3~guile build_system=generic arch=linux-sles15-zen3
module load gmake/4.4.1-gcc-12.3-i7fdisk
# berkeley-db@=18.1.40%rocmcc@=6.3.1+cxx~docs+stl build_system=autotools patches=26090f4,b231fcc arch=linux-sles15-zen3
module load berkeley-db/18.1.40-rocmcc-6.3.1-snty35e
# libiconv@=1.17%rocmcc@=6.3.1 build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load libiconv/1.17-rocmcc-6.3.1-fu3sqe6
# diffutils@=3.10%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-zen3
module load diffutils/3.10-rocmcc-6.3.1-c4sfpr4
# bzip2@=1.0.8%rocmcc@=6.3.1~debug~pic+shared build_system=generic arch=linux-sles15-zen3
module load bzip2/1.0.8-rocmcc-6.3.1-ghzkzol
# pkgconf@=2.3.0%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-zen3
module load pkgconf/2.3.0-rocmcc-6.3.1-jxagujx
# ncurses@=6.5%rocmcc@=6.3.1~symlinks+termlib abi=none build_system=autotools patches=7a351bc arch=linux-sles15-zen3
module load ncurses/6.5-rocmcc-6.3.1-i2h3up6
# readline@=8.2%rocmcc@=6.3.1 build_system=autotools patches=1ea4349,24f587b,3d9885e,5911a5b,622ba38,6c8adf8,758e2ec,79572ee,a177edc,bbf97f1,c7b45ff,e0013d9,e065038 arch=linux-sles15-zen3
module load readline/8.2-rocmcc-6.3.1-snonwhj
# gdbm@=1.23%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-zen3
module load gdbm/1.23-rocmcc-6.3.1-o5dm5k3
# zlib-ng@=2.2.3%rocmcc@=6.3.1+compat+new_strategies+opt+pic+shared build_system=autotools arch=linux-sles15-zen3
module load zlib-ng/2.2.3-rocmcc-6.3.1-qagngsw
# perl@=5.40.0%rocmcc@=6.3.1+cpanm+opcode+open+shared+threads build_system=generic arch=linux-sles15-zen3
module load perl/5.40.0-rocmcc-6.3.1-zbot3g3
# openblas@=0.3.20%gcc@=12.3~bignuma~consistent_fpcsr+dynamic_dispatch~ilp64+locking+pic+shared build_system=makefile patches=9f12903 symbol_suffix=none threads=none arch=linux-sles15-zen3
module load openblas/0.3.20-gcc-12.3-5lo7g62
# coinhsl@=2019.05.21%gcc@=12.3+blas build_system=autotools arch=linux-sles15-zen3
module load coinhsl/2019.05.21-gcc-12.3-2jugtcw
# hipblas@=6.3.1%rocmcc@=6.3.1~asan~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make patches=8d71578,b05b34b arch=linux-sles15-zen3
module load hipblas/6.3.1-rocmcc-6.3.1-gbtlgcy
# hiprand@=6.3.1%rocmcc@=6.3.1~asan~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load hiprand/6.3.1-rocmcc-6.3.1-fz3jffm
# hipsparse@=6.3.1%rocmcc@=6.3.1~asan~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load hipsparse/6.3.1-rocmcc-6.3.1-lwhon6d
# rocm-core@=6.3.1%rocmcc@=6.3.1~asan build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load rocm-core/6.3.1-rocmcc-6.3.1-3fytkbd
# magma@=2.8.0%rocmcc@=6.3.1~cuda+fortran~ipo+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load magma/2.8.0-rocmcc-6.3.1-gk7atiu
# metis@=5.1.0%rocmcc@=6.3.1~gdb~int64~no_warning~real64+shared build_system=cmake build_type=Release generator=make patches=4991da9,93a7903 arch=linux-sles15-zen3
module load metis/5.1.0-rocmcc-6.3.1-eqv3lql
# rocprim@=6.3.1%rocmcc@=6.3.1~asan amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load rocprim/6.3.1-rocmcc-6.3.1-r2wbqsu
# raja@=0.14.0%rocmcc@=6.3.1~cuda~desul~examples~exercises~ipo~omptarget~omptask~openmp~plugins+rocm~run-all-tests+shared~sycl~tests~vectorization amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load raja/0.14.0-rocmcc-6.3.1-goyrfcj
# libsigsegv@=2.14%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-zen3
module load libsigsegv/2.14-rocmcc-6.3.1-wmk7xpl
# m4@=1.4.19%rocmcc@=6.3.1+sigsegv build_system=autotools patches=9dc5fbd,bfdffa7 arch=linux-sles15-zen3
module load m4/1.4.19-rocmcc-6.3.1-vlezrbi
# autoconf@=2.72%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-zen3
module load autoconf/2.72-rocmcc-6.3.1-wvqalwi
# automake@=1.16.5%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-zen3
module load automake/1.16.5-rocmcc-6.3.1-tymrcpv
# xz@=5.4.6%rocmcc@=6.3.1~pic build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load xz/5.4.6-rocmcc-6.3.1-uam3imj
# libxml2@=2.13.5%rocmcc@=6.3.1~http+pic~python+shared build_system=autotools arch=linux-sles15-zen3
module load libxml2/2.13.5-rocmcc-6.3.1-42swsdo
# pigz@=2.8%rocmcc@=6.3.1 build_system=makefile arch=linux-sles15-zen3
module load pigz/2.8-rocmcc-6.3.1-vhc4ve7
# zstd@=1.5.6%rocmcc@=6.3.1+programs build_system=makefile compression=none libs=shared,static arch=linux-sles15-zen3
module load zstd/1.5.6-rocmcc-6.3.1-mvni44r
# tar@=1.35%rocmcc@=6.3.1 build_system=autotools zip=pigz arch=linux-sles15-zen3
module load tar/1.35-rocmcc-6.3.1-5hmruwv
# gettext@=0.23.1%rocmcc@=6.3.1+bzip2+curses+git~libunistring+libxml2+pic+shared+tar+xz build_system=autotools arch=linux-sles15-zen3
module load gettext/0.23.1-rocmcc-6.3.1-olazbjn
# findutils@=4.10.0%rocmcc@=6.3.1 build_system=autotools patches=440b954 arch=linux-sles15-zen3
module load findutils/4.10.0-rocmcc-6.3.1-xfx74g5
# libtool@=2.4.7%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-zen3
module load libtool/2.4.7-rocmcc-6.3.1-vuyual7
# gmp@=6.3.0%rocmcc@=6.3.1+cxx build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load gmp/6.3.0-rocmcc-6.3.1-wn2ocfd
# autoconf-archive@=2023.02.20%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-zen3
module load autoconf-archive/2023.02.20-rocmcc-6.3.1-cr25bvg
# texinfo@=7.1%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-zen3
module load texinfo/7.1-rocmcc-6.3.1-a7zkg2r
# mpfr@=4.2.1%rocmcc@=6.3.1 build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load mpfr/4.2.1-rocmcc-6.3.1-4zfn7vk
# suite-sparse@=7.8.3%rocmcc@=6.3.1~cuda~graphblas~openmp+pic build_system=generic arch=linux-sles15-zen3
module load suite-sparse/7.8.3-rocmcc-6.3.1-ldbpuyy
# umpire@=6.0.0%rocmcc@=6.3.1~asan~backtrace+c~cuda~dev_benchmarks~device_alloc~deviceconst~examples+fmt_header_only~fortran~ipc_shmem~ipo~mpi~numa~omptarget~openmp+rocm~sanitizer_tests+shared~sqlite_experimental~tools~werror amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make tests=none arch=linux-sles15-zen3
module load umpire/6.0.0-rocmcc-6.3.1-4g4w7pb
# hiop@=develop%rocmcc@=6.3.1~cuda~deepchecking~ginkgo~ipo~jsrun+kron+mpi+raja+rocm~shared+sparse amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make patches=bb62ae1 arch=linux-sles15-zen3
module load hiop/develop-rocmcc-6.3.1-kexanvc
# ipopt@=3.12.10%rocmcc@=6.3.1+coinhsl~debug~java~metis~mumps build_system=autotools arch=linux-sles15-zen3
module load ipopt/3.12.10-rocmcc-6.3.1-gnbtqod
# python@=3.11.5%rocmcc@=6.3.1+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=13fa8bf,b0615b2,ebdca64,f2fd060 arch=linux-sles15-zen3
module load python/3.11.5-rocmcc-6.3.1-a47nn7m
# petsc@=3.22.2%rocmcc@=6.3.1~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~sycl~tetgen~trilinos~valgrind~zoltan build_system=generic clanguage=C memalign=none arch=linux-sles15-zen3
module load petsc/3.22.2-rocmcc-6.3.1-gulnlvd
# exago@=develop%rocmcc@=6.3.1~cuda+hiop~ipo+ipopt+logging+mpi~python+raja+rocm amdgpu_target=gfx90a build_system=cmake build_type=Release dev_path=/lustre/orion/scratch/nkouk/stf006/Codes/ExaGO generator=make arch=linux-sles15-zen3
## module load exago/develop-rocmcc-6.3.1-5tcfzm4
