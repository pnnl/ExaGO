module use -a /lustre/orion/eng145/world-shared/spack-install/modules/linux-sles15-x86_64
# cmake@=3.27.9%rocmcc@=6.3.1~doc+ncurses+ownlibs~qtgui build_system=generic build_type=Release patches=dbc3892 arch=linux-sles15-x86_64
module load cmake/3.27.9-rocmcc-6.3.1-se2svx5
# glibc@=2.31%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-x86_64
module load glibc/2.31-rocmcc-6.3.1-lrviqsu
# blt@=0.4.1%rocmcc@=6.3.1 build_system=generic arch=linux-sles15-x86_64
module load blt/0.4.1-rocmcc-6.3.1-vazvlwp
# gmake@=4.4.1%rocmcc@=6.3.1~guile build_system=generic arch=linux-sles15-x86_64
module load gmake/4.4.1-rocmcc-6.3.1-qrwqcmf
# hip@=6.3.1%rocmcc@=6.3.1~asan~cuda+rocm build_system=cmake build_type=Release generator=make patches=1f65dfe arch=linux-sles15-x86_64
module load hip/6.3.1-rocmcc-6.3.1-adfbkqz
# hsa-rocr-dev@=6.3.1%rocmcc@=6.3.1~asan+image+shared build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load hsa-rocr-dev/6.3.1-rocmcc-6.3.1-ssifoxs
# llvm-amdgpu@=6.3.1%rocmcc@=6.3.1~link_llvm_dylib~llvm_dylib+rocm-device-libs build_system=cmake build_type=Release generator=ninja patches=b4774ca arch=linux-sles15-x86_64
module load llvm-amdgpu/6.3.1-rocmcc-6.3.1-3iotzh4
# camp@=0.2.3%rocmcc@=6.3.1~cuda~ipo~omptarget~openmp+rocm~sycl~tests amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make patches=cb9e25b,f854571 arch=linux-sles15-x86_64
module load camp/0.2.3-rocmcc-6.3.1-t4nudcc
# cray-mpich@=8.1.28%rocmcc@=6.3.1+wrappers build_system=generic arch=linux-sles15-x86_64
module load cray-mpich/8.1.28-rocmcc-6.3.1-om4wzu5
# gcc-runtime@=12.3%gcc@=12.3 build_system=generic arch=linux-sles15-x86_64
module load gcc-runtime/12.3-gcc-12.3-qdgdadl
# gmake@=4.4.1%gcc@=12.3~guile build_system=generic arch=linux-sles15-x86_64
module load gmake/4.4.1-gcc-12.3-kfvip5a
# berkeley-db@=18.1.40%rocmcc@=6.3.1+cxx~docs+stl build_system=autotools patches=26090f4,b231fcc arch=linux-sles15-x86_64
module load berkeley-db/18.1.40-rocmcc-6.3.1-bgxurz7
# libiconv@=1.17%rocmcc@=6.3.1 build_system=autotools libs=shared,static arch=linux-sles15-x86_64
module load libiconv/1.17-rocmcc-6.3.1-ku73cc7
# diffutils@=3.10%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-x86_64
module load diffutils/3.10-rocmcc-6.3.1-zo2edcf
# bzip2@=1.0.8%rocmcc@=6.3.1~debug~pic+shared build_system=generic arch=linux-sles15-x86_64
module load bzip2/1.0.8-rocmcc-6.3.1-oruwcbe
# pkgconf@=2.3.0%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-x86_64
module load pkgconf/2.3.0-rocmcc-6.3.1-nahfiwg
# ncurses@=6.5%rocmcc@=6.3.1~symlinks+termlib abi=none build_system=autotools patches=7a351bc arch=linux-sles15-x86_64
module load ncurses/6.5-rocmcc-6.3.1-rarrrlo
# readline@=8.2%rocmcc@=6.3.1 build_system=autotools patches=1ea4349,24f587b,3d9885e,5911a5b,622ba38,6c8adf8,758e2ec,79572ee,a177edc,bbf97f1,c7b45ff,e0013d9,e065038 arch=linux-sles15-x86_64
module load readline/8.2-rocmcc-6.3.1-i4zgaib
# gdbm@=1.23%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-x86_64
module load gdbm/1.23-rocmcc-6.3.1-4j7bbmo
# zlib-ng@=2.2.3%rocmcc@=6.3.1+compat+new_strategies+opt+pic+shared build_system=autotools arch=linux-sles15-x86_64
module load zlib-ng/2.2.3-rocmcc-6.3.1-nykmkcx
# perl@=5.40.0%rocmcc@=6.3.1+cpanm+opcode+open+shared+threads build_system=generic arch=linux-sles15-x86_64
module load perl/5.40.0-rocmcc-6.3.1-3fcevvu
# openblas@=0.3.20%gcc@=12.3~bignuma~consistent_fpcsr+dynamic_dispatch~ilp64+locking+pic+shared build_system=makefile patches=9f12903 symbol_suffix=none threads=none arch=linux-sles15-x86_64
module load openblas/0.3.20-gcc-12.3-yg3jpuq
# coinhsl@=2019.05.21%gcc@=12.3+blas build_system=autotools arch=linux-sles15-x86_64
module load coinhsl/2019.05.21-gcc-12.3-wxscvfx
# hipblas@=6.3.1%rocmcc@=6.3.1~asan~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make patches=8d71578,b05b34b arch=linux-sles15-x86_64
module load hipblas/6.3.1-rocmcc-6.3.1-ax2xpo6
# hiprand@=6.3.1%rocmcc@=6.3.1~asan~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load hiprand/6.3.1-rocmcc-6.3.1-s3e37yw
# hipsparse@=6.3.1%rocmcc@=6.3.1~asan~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load hipsparse/6.3.1-rocmcc-6.3.1-lck4zcu
# rocm-core@=6.3.1%rocmcc@=6.3.1~asan build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load rocm-core/6.3.1-rocmcc-6.3.1-z2q67i7
# magma@=2.8.0%rocmcc@=6.3.1~cuda+fortran~ipo+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load magma/2.8.0-rocmcc-6.3.1-6hgzai5
# metis@=5.1.0%rocmcc@=6.3.1~gdb~int64~no_warning~real64+shared build_system=cmake build_type=Release generator=make patches=4991da9,93a7903 arch=linux-sles15-x86_64
module load metis/5.1.0-rocmcc-6.3.1-tyhp2lo
# rocprim@=6.3.1%rocmcc@=6.3.1~asan amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load rocprim/6.3.1-rocmcc-6.3.1-fllulif
# raja@=0.14.0%rocmcc@=6.3.1~cuda~desul~examples~exercises~ipo~omptarget~omptask~openmp~plugins+rocm~run-all-tests+shared~sycl~tests~vectorization amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-x86_64
module load raja/0.14.0-rocmcc-6.3.1-u37hjzf
# libsigsegv@=2.14%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-x86_64
module load libsigsegv/2.14-rocmcc-6.3.1-6xlcqf3
# m4@=1.4.19%rocmcc@=6.3.1+sigsegv build_system=autotools patches=9dc5fbd,bfdffa7 arch=linux-sles15-x86_64
module load m4/1.4.19-rocmcc-6.3.1-4rovv2p
# autoconf@=2.72%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-x86_64
module load autoconf/2.72-rocmcc-6.3.1-h6icbph
# automake@=1.16.5%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-x86_64
module load automake/1.16.5-rocmcc-6.3.1-3jxtpcx
# xz@=5.4.6%rocmcc@=6.3.1~pic build_system=autotools libs=shared,static arch=linux-sles15-x86_64
module load xz/5.4.6-rocmcc-6.3.1-ph3kr5g
# libxml2@=2.13.5%rocmcc@=6.3.1~http+pic~python+shared build_system=autotools arch=linux-sles15-x86_64
module load libxml2/2.13.5-rocmcc-6.3.1-zmwsnxn
# pigz@=2.8%rocmcc@=6.3.1 build_system=makefile arch=linux-sles15-x86_64
module load pigz/2.8-rocmcc-6.3.1-7c4aur3
# zstd@=1.5.6%rocmcc@=6.3.1+programs build_system=makefile compression=none libs=shared,static arch=linux-sles15-x86_64
module load zstd/1.5.6-rocmcc-6.3.1-77qkqle
# tar@=1.35%rocmcc@=6.3.1 build_system=autotools zip=pigz arch=linux-sles15-x86_64
module load tar/1.35-rocmcc-6.3.1-u25yaw6
# gettext@=0.23.1%rocmcc@=6.3.1+bzip2+curses+git~libunistring+libxml2+pic+shared+tar+xz build_system=autotools arch=linux-sles15-x86_64
module load gettext/0.23.1-rocmcc-6.3.1-he7pw5h
# findutils@=4.10.0%rocmcc@=6.3.1 build_system=autotools patches=440b954 arch=linux-sles15-x86_64
module load findutils/4.10.0-rocmcc-6.3.1-akq3qmr
# libtool@=2.4.7%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-x86_64
module load libtool/2.4.7-rocmcc-6.3.1-shv6b7n
# gmp@=6.3.0%rocmcc@=6.3.1+cxx build_system=autotools libs=shared,static arch=linux-sles15-x86_64
module load gmp/6.3.0-rocmcc-6.3.1-ugcvc4h
# autoconf-archive@=2023.02.20%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-x86_64
module load autoconf-archive/2023.02.20-rocmcc-6.3.1-gkm3cna
# texinfo@=7.1%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-x86_64
module load texinfo/7.1-rocmcc-6.3.1-rvpweq4
# mpfr@=4.2.1%rocmcc@=6.3.1 build_system=autotools libs=shared,static arch=linux-sles15-x86_64
module load mpfr/4.2.1-rocmcc-6.3.1-v36rzqo
# suite-sparse@=7.8.3%rocmcc@=6.3.1~cuda~graphblas~openmp+pic build_system=generic arch=linux-sles15-x86_64
module load suite-sparse/7.8.3-rocmcc-6.3.1-lbosvok
# umpire@=6.0.0%rocmcc@=6.3.1~asan~backtrace+c~cuda~dev_benchmarks~device_alloc~deviceconst~examples+fmt_header_only~fortran~ipc_shmem~ipo~mpi~numa~omptarget~openmp+rocm~sanitizer_tests+shared~sqlite_experimental~tools~werror amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make tests=none arch=linux-sles15-x86_64
module load umpire/6.0.0-rocmcc-6.3.1-gi6w3oz
# hiop@=develop%rocmcc@=6.3.1~cuda~deepchecking~ginkgo~ipo~jsrun+kron+mpi+raja+rocm~shared+sparse amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make patches=bb62ae1 arch=linux-sles15-x86_64
module load hiop/develop-rocmcc-6.3.1-ci7ymop
# ipopt@=3.12.10%rocmcc@=6.3.1+coinhsl~debug~java~metis~mumps build_system=autotools arch=linux-sles15-x86_64
module load ipopt/3.12.10-rocmcc-6.3.1-fctobiz
# python@=3.11.5%rocmcc@=6.3.1+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=13fa8bf,b0615b2,ebdca64,f2fd060 arch=linux-sles15-x86_64
module load python/3.11.5-rocmcc-6.3.1-gzzj7bl
# petsc@=3.22.2%rocmcc@=6.3.1~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~sycl~tetgen~trilinos~valgrind~zoltan build_system=generic clanguage=C memalign=none arch=linux-sles15-x86_64
module load petsc/3.22.2-rocmcc-6.3.1-mr3e5zv
# exago@=develop%rocmcc@=6.3.1~cuda+hiop~ipo+ipopt+logging+mpi~python+raja+rocm amdgpu_target=gfx90a build_system=cmake build_type=Release dev_path=/lustre/orion/scratch/nkouk/stf006/Codes/ExaGO generator=make arch=linux-sles15-x86_64
## module load exago/develop-rocmcc-6.3.1-qs2bkd7
