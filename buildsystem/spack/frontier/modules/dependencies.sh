module use -a /lustre/orion/stf006/world-shared/nkouk/exago-05-2025/spack-install/modules/linux-sles15-zen3
# cmake@=3.30.5%rocmcc@=6.3.1~doc+ncurses+ownlibs~qtgui build_system=generic build_type=Release patches=dbc3892 arch=linux-sles15-zen3
module load cmake/3.30.5-rocmcc-6.3.1-hv4mzvr
# glibc@=2.38%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-zen3
module load glibc/2.38-rocmcc-6.3.1-3hkbcwm
# blt@=0.4.1%rocmcc@=6.3.1 build_system=generic arch=linux-sles15-zen3
module load blt/0.4.1-rocmcc-6.3.1-ru66g72
# gmake@=4.4.1%rocmcc@=6.3.1~guile build_system=generic arch=linux-sles15-zen3
module load gmake/4.4.1-rocmcc-6.3.1-jqumzoc
# hip@=6.3.1%rocmcc@=6.3.1~asan~cuda+rocm build_system=cmake build_type=Release generator=make patches=1f65dfe arch=linux-sles15-zen3
module load hip/6.3.1-rocmcc-6.3.1-ar4rzg2
# hsa-rocr-dev@=6.3.1%rocmcc@=6.3.1~asan+image+shared build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load hsa-rocr-dev/6.3.1-rocmcc-6.3.1-wwlghld
# llvm-amdgpu@=6.3.1%rocmcc@=6.3.1~link_llvm_dylib~llvm_dylib+rocm-device-libs build_system=cmake build_type=Release generator=ninja patches=b4774ca arch=linux-sles15-zen3
module load llvm-amdgpu/6.3.1-rocmcc-6.3.1-xcjzoi6
# camp@=0.2.3%rocmcc@=6.3.1~cuda~ipo~omptarget~openmp+rocm~sycl~tests amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make patches=cb9e25b,f854571 arch=linux-sles15-zen3
module load camp/0.2.3-rocmcc-6.3.1-nwyc533
# cray-mpich@=8.1.31%rocmcc@=6.3.1+wrappers build_system=generic arch=linux-sles15-zen3
module load cray-mpich/8.1.31-rocmcc-6.3.1-bvy2lpv
# gcc-runtime@=14.2%gcc@=14.2 build_system=generic arch=linux-sles15-zen3
module load gcc-runtime/14.2-gcc-14.2-yokdhjr
# python@=3.11.7%gcc@=14.2+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=13fa8bf,b0615b2,ebdca64,f2fd060 arch=linux-sles15-zen3
module load python/3.11.7-gcc-14.2-wxklcnc
# gmake@=4.4.1%gcc@=14.2~guile build_system=generic arch=linux-sles15-zen3
module load gmake/4.4.1-gcc-14.2-367srkr
# re2c@=3.1%gcc@=14.2 build_system=autotools arch=linux-sles15-zen3
module load re2c/3.1-gcc-14.2-zzzjdxf
# ninja@=1.12.1%gcc@=14.2+re2c build_system=generic arch=linux-sles15-zen3
module load ninja/1.12.1-gcc-14.2-ru7pywf
# python-venv@=1.0%gcc@=14.2 build_system=generic arch=linux-sles15-zen3
module load python-venv/1.0-gcc-14.2-bkaefwo
# py-pip@=23.1.2%gcc@=14.2 build_system=generic arch=linux-sles15-zen3
module load py-pip/23.1.2-gcc-14.2-np2dct5
# py-setuptools@=69.2.0%gcc@=14.2 build_system=generic arch=linux-sles15-zen3
module load py-setuptools/69.2.0-gcc-14.2-ub2pryb
# py-wheel@=0.41.2%gcc@=14.2 build_system=generic arch=linux-sles15-zen3
module load py-wheel/0.41.2-gcc-14.2-xonikb7
# meson@=1.5.1%gcc@=14.2 build_system=python_pip patches=0f0b1bd arch=linux-sles15-zen3
module load meson/1.5.1-gcc-14.2-ddeivby
# metis@=5.1.0%rocmcc@=6.3.1~gdb~int64~ipo~real64+shared build_system=cmake build_type=Release generator=make patches=4991da9,93a7903 arch=linux-sles15-zen3
module load metis/5.1.0-rocmcc-6.3.1-xxxrqja
# berkeley-db@=18.1.40%rocmcc@=6.3.1+cxx~docs+stl build_system=autotools patches=26090f4,b231fcc arch=linux-sles15-zen3
module load berkeley-db/18.1.40-rocmcc-6.3.1-4c4xyjg
# libiconv@=1.17%rocmcc@=6.3.1 build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load libiconv/1.17-rocmcc-6.3.1-m52kzif
# diffutils@=3.10%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-zen3
module load diffutils/3.10-rocmcc-6.3.1-obnafli
# bzip2@=1.0.8%rocmcc@=6.3.1~debug~pic+shared build_system=generic arch=linux-sles15-zen3
module load bzip2/1.0.8-rocmcc-6.3.1-6opfml5
# pkgconf@=2.2.0%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-zen3
module load pkgconf/2.2.0-rocmcc-6.3.1-6qyuzby
# ncurses@=6.5%rocmcc@=6.3.1~symlinks+termlib abi=none build_system=autotools patches=7a351bc arch=linux-sles15-zen3
module load ncurses/6.5-rocmcc-6.3.1-xo7cmdh
# readline@=8.2%rocmcc@=6.3.1 build_system=autotools patches=bbf97f1 arch=linux-sles15-zen3
module load readline/8.2-rocmcc-6.3.1-jfidt3c
# gdbm@=1.23%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-zen3
module load gdbm/1.23-rocmcc-6.3.1-n2wrrdx
# zlib-ng@=2.2.1%rocmcc@=6.3.1+compat+new_strategies+opt+pic+shared build_system=autotools arch=linux-sles15-zen3
module load zlib-ng/2.2.1-rocmcc-6.3.1-ozhgif5
# perl@=5.40.0%rocmcc@=6.3.1+cpanm+opcode+open+shared+threads build_system=generic arch=linux-sles15-zen3
module load perl/5.40.0-rocmcc-6.3.1-yuniufu
# openblas@=0.3.20%gcc@=14.2~bignuma~consistent_fpcsr+dynamic_dispatch~ilp64+locking+pic+shared build_system=makefile patches=9f12903 symbol_suffix=none threads=none arch=linux-sles15-zen3
module load openblas/0.3.20-gcc-14.2-dqkevi2
# coinhsl@=2024.05.15%gcc@=14.2+metis~strip build_system=meson buildtype=release default_library=shared arch=linux-sles15-zen3
module load coinhsl/2024.05.15-gcc-14.2-3oirf5f
# hipblas@=6.3.1%rocmcc@=6.3.1~asan~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make patches=b05b34b arch=linux-sles15-zen3
module load hipblas/6.3.1-rocmcc-6.3.1-wxwuiie
# hiprand@=6.3.1%rocmcc@=6.3.1~asan~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load hiprand/6.3.1-rocmcc-6.3.1-67er6wz
# hipsparse@=6.3.1%rocmcc@=6.3.1~asan~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load hipsparse/6.3.1-rocmcc-6.3.1-k7565cr
# rocm-core@=6.3.1%rocmcc@=6.3.1~asan build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load rocm-core/6.3.1-rocmcc-6.3.1-3fytkbd
# magma@=2.8.0%rocmcc@=6.3.1~cuda+fortran~ipo+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load magma/2.8.0-rocmcc-6.3.1-bkz6mcz
# rocprim@=6.3.1%rocmcc@=6.3.1~asan amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load rocprim/6.3.1-rocmcc-6.3.1-r2wbqsu
# raja@=0.14.0%rocmcc@=6.3.1~cuda~desul~examples~exercises~ipo~omptarget~omptask~openmp~plugins+rocm~run-all-tests+shared~sycl~tests~vectorization amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load raja/0.14.0-rocmcc-6.3.1-baohou3
# libsigsegv@=2.14%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-zen3
module load libsigsegv/2.14-rocmcc-6.3.1-h2zoyuo
# m4@=1.4.19%rocmcc@=6.3.1+sigsegv build_system=autotools patches=9dc5fbd,bfdffa7 arch=linux-sles15-zen3
module load m4/1.4.19-rocmcc-6.3.1-6jkfljq
# autoconf@=2.72%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-zen3
module load autoconf/2.72-rocmcc-6.3.1-7g3s7zf
# automake@=1.16.5%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-zen3
module load automake/1.16.5-rocmcc-6.3.1-tnehnbg
# findutils@=4.9.0%rocmcc@=6.3.1 build_system=autotools patches=440b954 arch=linux-sles15-zen3
module load findutils/4.9.0-rocmcc-6.3.1-5rp3cmp
# libtool@=2.4.7%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-zen3
module load libtool/2.4.7-rocmcc-6.3.1-cikm7tr
# gmp@=6.3.0%rocmcc@=6.3.1+cxx build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load gmp/6.3.0-rocmcc-6.3.1-t4wddaf
# autoconf-archive@=2023.02.20%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-zen3
module load autoconf-archive/2023.02.20-rocmcc-6.3.1-ejmrdvy
# xz@=5.4.6%rocmcc@=6.3.1~pic build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load xz/5.4.6-rocmcc-6.3.1-tbfzi5n
# libxml2@=2.13.4%rocmcc@=6.3.1+pic~python+shared build_system=autotools arch=linux-sles15-zen3
module load libxml2/2.13.4-rocmcc-6.3.1-okrmncg
# pigz@=2.8%rocmcc@=6.3.1 build_system=makefile arch=linux-sles15-zen3
module load pigz/2.8-rocmcc-6.3.1-lywcdd5
# zstd@=1.5.6%rocmcc@=6.3.1+programs build_system=makefile compression=none libs=shared,static arch=linux-sles15-zen3
module load zstd/1.5.6-rocmcc-6.3.1-5qaccww
# tar@=1.34%rocmcc@=6.3.1 build_system=autotools zip=pigz arch=linux-sles15-zen3
module load tar/1.34-rocmcc-6.3.1-gsdfagp
# gettext@=0.22.5%rocmcc@=6.3.1+bzip2+curses+git~libunistring+libxml2+pic+shared+tar+xz build_system=autotools arch=linux-sles15-zen3
module load gettext/0.22.5-rocmcc-6.3.1-khcpe53
# texinfo@=7.1%rocmcc@=6.3.1 build_system=autotools arch=linux-sles15-zen3
module load texinfo/7.1-rocmcc-6.3.1-5fjqqw6
# mpfr@=4.2.1%rocmcc@=6.3.1 build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load mpfr/4.2.1-rocmcc-6.3.1-dbj45ac
# suite-sparse@=7.7.0%rocmcc@=6.3.1~cuda~graphblas~openmp+pic build_system=generic arch=linux-sles15-zen3
module load suite-sparse/7.7.0-rocmcc-6.3.1-iq7sq7v
# umpire@=6.0.0%rocmcc@=6.3.1~asan~backtrace+c~cuda~dev_benchmarks~device_alloc~deviceconst~examples+fmt_header_only~fortran~ipc_shmem~ipo~mpi~numa~omptarget~openmp+rocm~sanitizer_tests+shared~sqlite_experimental~tools~werror amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make tests=none arch=linux-sles15-zen3
module load umpire/6.0.0-rocmcc-6.3.1-3zjinbv
# hiop@=1.1.0%rocmcc@=6.3.1~cuda~deepchecking~ginkgo~ipo~jsrun+kron+mpi+raja+rocm~shared+sparse amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make patches=bb62ae1 arch=linux-sles15-zen3
module load hiop/1.1.0-rocmcc-6.3.1-iqwppvk
# ipopt@=3.14.14%rocmcc@=6.3.1+coinhsl~debug~java~metis~mumps build_system=autotools arch=linux-sles15-zen3
module load ipopt/3.14.14-rocmcc-6.3.1-djalso6
# parmetis@=4.0.3%rocmcc@=6.3.1~gdb~int64~ipo+shared build_system=cmake build_type=Release generator=make patches=4f89253,50ed208,704b84f arch=linux-sles15-zen3
module load parmetis/4.0.3-rocmcc-6.3.1-amlsqtx
# petsc@=3.22.1%rocmcc@=6.3.1~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind+metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~sycl~tetgen~trilinos~valgrind~zoltan build_system=generic clanguage=C memalign=none arch=linux-sles15-zen3
module load petsc/3.22.1-rocmcc-6.3.1-krpsfil
# exago@=develop%rocmcc@=6.3.1~cuda+hiop~ipo+ipopt+logging+mpi~python+raja+rocm amdgpu_target=gfx90a build_system=cmake build_type=Release dev_path=/lustre/orion/scratch/nkouk/stf006/Codes/ExaGO generator=make arch=linux-sles15-zen3
## module load exago/develop-rocmcc-6.3.1-rn7oc3p
