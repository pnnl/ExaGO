module use -a /lustre/orion/stf006/world-shared/nkouk/exago-02-2026/spack-install/modules/linux-sles15-zen3
## excluded or missing from upstream: cmake@=3.30.5~doc+ncurses+ownlibs~qtgui build_system=generic build_type=Release platform=linux os=sles15 target=x86_64
# compiler-wrapper@=1.0 build_system=generic platform=linux os=sles15 target=zen3
module load compiler-wrapper/1.0-none-none-rfhw7ey
## excluded or missing from upstream: gcc@=14.2.0~binutils+bootstrap~graphite~mold~nvptx~piclibs~profiled~strip build_system=autotools build_type=RelWithDebInfo languages:='c,c++,fortran' platform=linux os=sles15 target=x86_64
## excluded or missing from upstream: glibc@=2.38 build_system=autotools platform=linux os=sles15 target=x86_64
# gcc-runtime@=14.2.0 build_system=generic platform=linux os=sles15 target=zen3
module load gcc-runtime/14.2.0-none-none-ddc453r
# blt@=0.7.1 build_system=generic platform=linux os=sles15 target=zen3
module load blt/0.7.1-gcc-14.2.0-ckjjt7m
# gmake@=4.4.1~guile build_system=generic platform=linux os=sles15 target=zen3
module load gmake/4.4.1-gcc-14.2.0-kqtm7j7
## excluded or missing from upstream: hip@=6.3.1~asan~cuda+rocm build_system=cmake build_type=Release generator=make platform=linux os=sles15 target=x86_64
## excluded or missing from upstream: hsa-rocr-dev@=6.3.1~asan+image+shared build_system=cmake build_type=Release generator=make platform=linux os=sles15 target=x86_64
## excluded or missing from upstream: llvm-amdgpu@=6.3.1~link_llvm_dylib~llvm_dylib+rocm-device-libs build_system=cmake build_type=Release generator=ninja languages:='c,c++' platform=linux os=sles15 target=x86_64
# camp@=2024.07.0~cuda~ipo~omptarget~openmp+rocm~sycl~tests amdgpu_target:=gfx90a build_system=cmake build_type=Release commit=0f07de4240c42e0b38a8d872a20440cb4b33d9f5 generator=make platform=linux os=sles15 target=zen3
module load camp/2024.07.0-gcc-14.2.0-6zkynai
## excluded or missing from upstream: cray-mpich@=8.1.31~cuda~rocm+wrappers build_system=generic platform=linux os=sles15 target=x86_64
# fmt@=11.0.2~ipo+pic~shared build_system=cmake build_type=Release cxxstd=11 generator=make platform=linux os=sles15 target=zen3
module load fmt/11.0.2-gcc-14.2.0-hmdf2um
## excluded or missing from upstream: python@=3.11.7+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic platform=linux os=sles15 target=x86_64
# re2c@=3.1 build_system=autotools platform=linux os=sles15 target=zen3
module load re2c/3.1-gcc-14.2.0-kjfvtlu
# ninja@=1.13.0+re2c build_system=generic platform=linux os=sles15 target=zen3
module load ninja/1.13.0-gcc-14.2.0-3xtgf4s
# python-venv@=1.0 build_system=generic platform=linux os=sles15 target=zen3
module load python-venv/1.0-none-none-w2bm77w
# py-pip@=25.1.1 build_system=generic platform=linux os=sles15 target=zen3
module load py-pip/25.1.1-none-none-qeddnbh
# py-setuptools@=80.9.0 build_system=generic platform=linux os=sles15 target=zen3
module load py-setuptools/80.9.0-none-none-mjws5td
# py-wheel@=0.45.1 build_system=generic platform=linux os=sles15 target=zen3
module load py-wheel/0.45.1-none-none-mrtxpal
# meson@=1.8.5 build_system=python_pip patches:=0f0b1bd platform=linux os=sles15 target=zen3
module load meson/1.8.5-none-none-3oz75tc
# metis@=5.1.0~gdb~int64~ipo~no_warning~real64+shared build_system=cmake build_type=Release generator=make patches:=4991da9,93a7903,b1225da platform=linux os=sles15 target=zen3
module load metis/5.1.0-gcc-14.2.0-yujvyqz
# berkeley-db@=18.1.40+cxx~docs+stl build_system=autotools patches:=26090f4,b231fcc platform=linux os=sles15 target=zen3
module load berkeley-db/18.1.40-gcc-14.2.0-hs3yf5g
# libiconv@=1.18 build_system=autotools libs:=shared,static platform=linux os=sles15 target=zen3
module load libiconv/1.18-gcc-14.2.0-4cxmjca
# diffutils@=3.12 build_system=autotools platform=linux os=sles15 target=zen3
module load diffutils/3.12-gcc-14.2.0-2jwfpug
# bzip2@=1.0.8~debug~pic+shared build_system=generic platform=linux os=sles15 target=zen3
module load bzip2/1.0.8-gcc-14.2.0-bxs2wih
# pkgconf@=2.5.1 build_system=autotools platform=linux os=sles15 target=zen3
module load pkgconf/2.5.1-gcc-14.2.0-4bsp2cj
# ncurses@=6.5-20250705~symlinks+termlib abi=none build_system=autotools patches:=7a351bc platform=linux os=sles15 target=zen3
module load ncurses/6.5-20250705-gcc-14.2.0-5k2higm
# readline@=8.3 build_system=autotools patches:=21f0a03 platform=linux os=sles15 target=zen3
module load readline/8.3-gcc-14.2.0-66mb24c
# gdbm@=1.25 build_system=autotools platform=linux os=sles15 target=zen3
module load gdbm/1.25-gcc-14.2.0-w67thwv
# zlib-ng@=2.2.4+compat+new_strategies+opt+pic+shared build_system=autotools platform=linux os=sles15 target=zen3
module load zlib-ng/2.2.4-gcc-14.2.0-mm3wqqs
# perl@=5.42.0+cpanm+opcode+open+shared+threads build_system=generic platform=linux os=sles15 target=zen3
module load perl/5.42.0-gcc-14.2.0-g7cf6rz
# openblas@=0.3.20~bignuma~consistent_fpcsr+dynamic_dispatch~ilp64+locking+pic+shared build_system=makefile patches:=9f12903 symbol_suffix=none threads=none platform=linux os=sles15 target=zen3
module load openblas/0.3.20-gcc-14.2.0-holyhf3
# coinhsl@=2024.05.15+metis~strip build_system=meson buildtype=release default_library:=shared platform=linux os=sles15 target=zen3
module load coinhsl/2024.05.15-gcc-14.2.0-ke5biep
## excluded or missing from upstream: hipblas@=6.3.1~asan~cuda+rocm amdgpu_target:=auto build_system=cmake build_type=Release generator=make platform=linux os=sles15 target=x86_64
## excluded or missing from upstream: hiprand@=6.3.1~asan~cuda+rocm amdgpu_target:=auto build_system=cmake build_type=Release generator=make platform=linux os=sles15 target=x86_64
## excluded or missing from upstream: hipsparse@=6.3.1~asan~cuda+rocm amdgpu_target:=auto build_system=cmake build_type=Release generator=make platform=linux os=sles15 target=x86_64
## excluded or missing from upstream: rocm-core@=6.3.1~asan build_system=cmake build_type=Release generator=make platform=linux os=sles15 target=x86_64
# magma@=2.8.0~cuda+fortran~ipo+rocm+shared amdgpu_target:=gfx90a build_system=cmake build_type=Release generator=make platform=linux os=sles15 target=zen3
module load magma/2.8.0-gcc-14.2.0-ikok5nn
## excluded or missing from upstream: rocprim@=6.3.1~asan amdgpu_target:=auto build_system=cmake build_type=Release generator=make platform=linux os=sles15 target=x86_64
# raja@=2024.07.0~cuda~desul~examples~exercises~gpu-profiling~ipo~lowopttest~omptarget~omptask~openmp~plugins+rocm~run-all-tests~shared~sycl~tests~vectorization amdgpu_target:=gfx90a build_system=cmake build_type=Release commit=4d7fcba55ebc7cb972b7cc9f6778b48e43792ea1 generator=make platform=linux os=sles15 target=zen3
module load raja/2024.07.0-gcc-14.2.0-3e7fwpj
# libsigsegv@=2.14 build_system=autotools platform=linux os=sles15 target=zen3
module load libsigsegv/2.14-gcc-14.2.0-j5nwpxr
# m4@=1.4.20+sigsegv build_system=autotools platform=linux os=sles15 target=zen3
module load m4/1.4.20-gcc-14.2.0-eyfod7z
# autoconf@=2.72 build_system=autotools platform=linux os=sles15 target=zen3
module load autoconf/2.72-none-none-e32o2ms
# automake@=1.16.5 build_system=autotools platform=linux os=sles15 target=zen3
module load automake/1.16.5-gcc-14.2.0-qupgkqr
# xz@=5.6.3~pic build_system=autotools libs:=shared,static platform=linux os=sles15 target=zen3
module load xz/5.6.3-gcc-14.2.0-n3iixvg
# libxml2@=2.13.5~http+pic~python+shared build_system=autotools platform=linux os=sles15 target=zen3
module load libxml2/2.13.5-gcc-14.2.0-kkshkbr
# pigz@=2.8 build_system=makefile platform=linux os=sles15 target=zen3
module load pigz/2.8-gcc-14.2.0-d4ifmkj
# zstd@=1.5.7+programs build_system=makefile compression:=none libs:=shared,static platform=linux os=sles15 target=zen3
module load zstd/1.5.7-gcc-14.2.0-r3pky4w
# tar@=1.35 build_system=autotools zip=pigz platform=linux os=sles15 target=zen3
module load tar/1.35-gcc-14.2.0-ae4witc
# gettext@=0.23.1+bzip2+curses+git~libunistring+libxml2+pic+shared+tar+xz build_system=autotools platform=linux os=sles15 target=zen3
module load gettext/0.23.1-gcc-14.2.0-g4d6jeh
# findutils@=4.10.0 build_system=autotools patches:=440b954 platform=linux os=sles15 target=zen3
module load findutils/4.10.0-gcc-14.2.0-r2qqgvh
# libtool@=2.4.7 build_system=autotools platform=linux os=sles15 target=zen3
module load libtool/2.4.7-gcc-14.2.0-uiult5c
# gmp@=6.3.0+cxx build_system=autotools libs:=shared,static platform=linux os=sles15 target=zen3
module load gmp/6.3.0-gcc-14.2.0-hf6zhfr
# autoconf-archive@=2023.02.20 build_system=autotools platform=linux os=sles15 target=zen3
module load autoconf-archive/2023.02.20-none-none-dw7lrnr
# texinfo@=7.2~xs build_system=autotools platform=linux os=sles15 target=zen3
module load texinfo/7.2-gcc-14.2.0-wtmf5g5
# mpfr@=4.2.1 build_system=autotools libs:=shared,static patches:=3ec29a6 platform=linux os=sles15 target=zen3
module load mpfr/4.2.1-gcc-14.2.0-p4r33ds
# suite-sparse@=7.8.3~cuda~graphblas~openmp+pic build_system=generic platform=linux os=sles15 target=zen3
module load suite-sparse/7.8.3-gcc-14.2.0-3f2u2ru
# umpire@=2024.07.0~asan~backtrace+c~cuda~dev_benchmarks~device_alloc~deviceconst~examples+fmt_header_only~fortran~ipc_shmem~ipo~mpi~mpi3_shmem~numa~omptarget~openmp+rocm~sanitizer_tests+shared~sqlite_experimental~tools~werror amdgpu_target:=gfx90a build_system=cmake build_type=Release commit=abd729f40064175e999a83d11d6b073dac4c01d2 generator=make tests=none platform=linux os=sles15 target=zen3
module load umpire/2024.07.0-gcc-14.2.0-tzf24hf
# hiop@=1.1.1~axom~cuda~deepchecking~ginkgo~ipo~jsrun+kron+mpi+raja+rocm~shared+sparse amdgpu_target:=gfx90a build_system=cmake build_type=Release commit=d8762e05150b2040a27f69d8bf6603f22190a869 generator=make platform=linux os=sles15 target=zen3
module load hiop/1.1.1-gcc-14.2.0-beu4x7r
# ipopt@=3.14.14+coinhsl~debug~java~metis~mumps build_system=autotools platform=linux os=sles15 target=zen3
module load ipopt/3.14.14-gcc-14.2.0-gj4qlx4
# parmetis@=4.0.3~gdb~int64~ipo+shared build_system=cmake build_type=Release generator=make patches:=4f89253,50ed208,704b84f platform=linux os=sles15 target=zen3
module load parmetis/4.0.3-gcc-14.2.0-xatjoum
# petsc@=3.24.1~X~batch~cgns~complex~cuda~debug+double+examples~exodusii~fftw+fortran+fortran-bindings~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind+metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~sycl~tetgen~trilinos~valgrind~zoltan build_system=generic clanguage=C memalign=none patches:=fa5ef56 platform=linux os=sles15 target=zen3
module load petsc/3.24.1-gcc-14.2.0-4izpiqe
# spdlog@=1.15.0~ipo+shared build_system=cmake build_type=Release generator=make patches:=5ed92f4,fd4cbb1,fdc325d platform=linux os=sles15 target=zen3
module load spdlog/1.15.0-gcc-14.2.0-56k2ubv
# exago@=develop~cuda+hiop~ipo+ipopt+logging+mpi~python+raja+rocm amdgpu_target:=gfx90a build_system=cmake build_type=Release dev_path=/lustre/orion/scratch/nkouk/stf006/Tmp/ExaGO generator=make platform=linux os=sles15 target=zen3
## module load exago/develop-gcc-14.2.0-lj7562k
