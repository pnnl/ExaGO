module use -a /gpfs/alpine2/stf006/world-shared/nkouk/exago-spack-install/summit-modules/linux-rhel8-power9le
# cmake@=3.23.2%gcc@=10.2.0~doc+ncurses+ownlibs build_system=generic build_type=Release patches=dbc3892 arch=linux-rhel8-power9le
module load cmake/3.23.2-gcc-10.2.0-5qroxa4
# glibc@=2.28%gcc@=10.2.0 build_system=autotools arch=linux-rhel8-power9le
module load glibc/2.28-gcc-10.2.0-gla3q5m
# gcc-runtime@=10.2.0%gcc@=10.2.0 build_system=generic arch=linux-rhel8-power9le
module load gcc-runtime/10.2.0-gcc-10.2.0-swvql6q
# blt@=0.4.1%gcc@=10.2.0 build_system=generic arch=linux-rhel8-power9le
module load blt/0.4.1-gcc-10.2.0-au5462o
# cub@=2.1.0%gcc@=10.2.0 build_system=generic arch=linux-rhel8-power9le
module load cub/2.1.0-gcc-10.2.0-cntbldc
# gmake@=4.4.1%gcc@=10.2.0~guile build_system=generic arch=linux-rhel8-power9le
module load gmake/4.4.1-gcc-10.2.0-dnqn5ag
# gnuconfig@=2022-09-17%gcc@=10.2.0 build_system=generic arch=linux-rhel8-power9le
module load gnuconfig/2022-09-17-gcc-10.2.0-qyktplu
# libiconv@=1.17%gcc@=10.2.0 build_system=autotools libs=shared,static arch=linux-rhel8-power9le
module load libiconv/1.17-gcc-10.2.0-rjjfi2f
# pkgconf@=2.2.0%gcc@=10.2.0 build_system=autotools arch=linux-rhel8-power9le
module load pkgconf/2.2.0-gcc-10.2.0-6gw2wxr
# xz@=5.4.6%gcc@=10.2.0~pic build_system=autotools libs=shared,static arch=linux-rhel8-power9le
module load xz/5.4.6-gcc-10.2.0-zb5x4it
# zlib-ng@=2.2.1%gcc@=10.2.0+compat+new_strategies+opt+pic+shared build_system=autotools arch=linux-rhel8-power9le
module load zlib-ng/2.2.1-gcc-10.2.0-7lkqxgk
# libxml2@=2.10.3%gcc@=10.2.0+pic~python+shared build_system=autotools arch=linux-rhel8-power9le
module load libxml2/2.10.3-gcc-10.2.0-rufygxt
# cuda@=11.4.2%gcc@=10.2.0~allow-unsupported-compilers~dev build_system=generic arch=linux-rhel8-power9le
module load cuda/11.4.2-gcc-10.2.0-uqlz4r6
# camp@=0.2.3%gcc@=10.2.0+cuda~ipo~openmp~rocm~tests build_system=cmake build_type=Release cuda_arch=70 generator=make patches=cb9e25b arch=linux-rhel8-power9le
module load camp/0.2.3-gcc-10.2.0-p2kfyr7
# ginkgo@=1.5.0.glu_experimental%gcc@=10.2.0+cuda~develtools~full_optimizations~hwloc~ipo~mpi+openmp~rocm+shared~sycl build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-rhel8-power9le
module load ginkgo/1.5.0.glu_experimental-gcc-10.2.0-6ql6qmh
# openblas@=0.3.17%gcc@=10.2.0~bignuma~consistent_fpcsr+dynamic_dispatch~ilp64+locking+pic+shared build_system=makefile symbol_suffix=none threads=none arch=linux-rhel8-power9le
module load openblas/0.3.17-gcc-10.2.0-wuwjnlg
# coinhsl@=2019.05.21%gcc@=10.2.0+blas build_system=autotools arch=linux-rhel8-power9le
module load coinhsl/2019.05.21-gcc-10.2.0-yrlo3sp
# magma@=2.7.2%gcc@=10.2.0+cuda+fortran~ipo~rocm+shared build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-rhel8-power9le
module load magma/2.7.2-gcc-10.2.0-xgh4g6c
# metis@=5.1.0%gcc@=10.2.0~gdb~int64~real64+shared build_system=cmake build_type=Release generator=make patches=4991da9,93a7903,b1225da arch=linux-rhel8-power9le
module load metis/5.1.0-gcc-10.2.0-hciojcu
# raja@=0.14.0%gcc@=10.2.0+cuda~desul+examples+exercises~ipo~omptask~openmp~plugins~rocm~run-all-tests+shared~tests~vectorization build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-rhel8-power9le
module load raja/0.14.0-gcc-10.2.0-gmsy36e
# spectrum-mpi@=10.4.0.3-20210112%gcc@=10.2.0 build_system=bundle arch=linux-rhel8-power9le
module load spectrum-mpi/10.4.0.3-20210112-gcc-10.2.0-iqladxy
# diffutils@=3.10%gcc@=10.2.0 build_system=autotools arch=linux-rhel8-power9le
module load diffutils/3.10-gcc-10.2.0-n3hree5
# libsigsegv@=2.14%gcc@=10.2.0 build_system=autotools arch=linux-rhel8-power9le
module load libsigsegv/2.14-gcc-10.2.0-uml73jg
# m4@=1.4.19%gcc@=10.2.0+sigsegv build_system=autotools patches=9dc5fbd,bfdffa7 arch=linux-rhel8-power9le
module load m4/1.4.19-gcc-10.2.0-jedmjig
# perl@=5.30.1%gcc@=10.2.0+cpanm+opcode+open+shared+threads build_system=generic arch=linux-rhel8-power9le
module load perl/5.30.1-gcc-10.2.0-durrgca
# autoconf@=2.72%gcc@=10.2.0 build_system=autotools arch=linux-rhel8-power9le
module load autoconf/2.72-gcc-10.2.0-nho5qoo
# automake@=1.16.5%gcc@=10.2.0 build_system=autotools arch=linux-rhel8-power9le
module load automake/1.16.5-gcc-10.2.0-7gn6zit
# findutils@=4.9.0%gcc@=10.2.0 build_system=autotools patches=440b954 arch=linux-rhel8-power9le
module load findutils/4.9.0-gcc-10.2.0-tearl3v
# libtool@=2.4.7%gcc@=10.2.0 build_system=autotools arch=linux-rhel8-power9le
module load libtool/2.4.7-gcc-10.2.0-hnmnqdp
# gmp@=6.3.0%gcc@=10.2.0+cxx build_system=autotools libs=shared,static arch=linux-rhel8-power9le
module load gmp/6.3.0-gcc-10.2.0-p3siics
# autoconf-archive@=2023.02.20%gcc@=10.2.0 build_system=autotools arch=linux-rhel8-power9le
module load autoconf-archive/2023.02.20-gcc-10.2.0-kfdcf4y
# bzip2@=1.0.8%gcc@=10.2.0~debug~pic+shared build_system=generic arch=linux-rhel8-power9le
module load bzip2/1.0.8-gcc-10.2.0-wddlbky
# ncurses@=6.5%gcc@=10.2.0~symlinks+termlib abi=none build_system=autotools patches=7a351bc arch=linux-rhel8-power9le
module load ncurses/6.5-gcc-10.2.0-3bkq5f4
# pigz@=2.8%gcc@=10.2.0 build_system=makefile arch=linux-rhel8-power9le
module load pigz/2.8-gcc-10.2.0-aejawbc
# zstd@=1.5.6%gcc@=10.2.0+programs build_system=makefile compression=none libs=shared,static arch=linux-rhel8-power9le
module load zstd/1.5.6-gcc-10.2.0-oyxbdrb
# tar@=1.34%gcc@=10.2.0 build_system=autotools zip=pigz arch=linux-rhel8-power9le
module load tar/1.34-gcc-10.2.0-3bm2duw
# gettext@=0.22.5%gcc@=10.2.0+bzip2+curses+git~libunistring+libxml2+pic+shared+tar+xz build_system=autotools arch=linux-rhel8-power9le
module load gettext/0.22.5-gcc-10.2.0-rmhxj6q
# texinfo@=7.1%gcc@=10.2.0 build_system=autotools arch=linux-rhel8-power9le
module load texinfo/7.1-gcc-10.2.0-taix7z2
# mpfr@=4.2.1%gcc@=10.2.0 build_system=autotools libs=shared,static arch=linux-rhel8-power9le
module load mpfr/4.2.1-gcc-10.2.0-2aeysfy
# suite-sparse@=7.7.0%gcc@=10.2.0~cuda~graphblas~openmp+pic build_system=generic arch=linux-rhel8-power9le
module load suite-sparse/7.7.0-gcc-10.2.0-f5odelu
# umpire@=6.0.0%gcc@=10.2.0~asan~backtrace~c+cuda~dev_benchmarks~device_alloc~deviceconst~examples~fortran~ipc_shmem~ipo~mpi~numa~openmp~openmp_target~rocm~sanitizer_tests~shared~sqlite_experimental~tools~werror build_system=cmake build_type=Release cuda_arch=70 generator=make tests=none arch=linux-rhel8-power9le
module load umpire/6.0.0-gcc-10.2.0-kcwtinj
# hiop@=develop%gcc@=10.2.0+cuda+cusolver_lu~deepchecking+ginkgo~ipo~jsrun+kron+mpi+raja~rocm~shared+sparse build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-rhel8-power9le
module load hiop/develop-gcc-10.2.0-eemr776
# ipopt@=3.12.10%gcc@=10.2.0+coinhsl~debug+metis~mumps build_system=autotools arch=linux-rhel8-power9le
module load ipopt/3.12.10-gcc-10.2.0-o7rvvcf
# hdf5@=1.14.3%gcc@=10.2.0~cxx~fortran~hl~ipo~java~map+mpi+shared~subfiling~szip~threadsafe+tools api=default build_system=cmake build_type=Release generator=make patches=82088c8 arch=linux-rhel8-power9le
module load hdf5/1.14.3-gcc-10.2.0-h376yw2
# hypre@=2.31.0%gcc@=10.2.0~caliper~complex~cublas~cuda~debug+fortran~gptune~gpu-aware-mpi~int64~internal-superlu~magma~mixedint+mpi~openmp~rocblas~rocm+shared~superlu-dist~sycl~umpire~unified-memory build_system=autotools precision=double arch=linux-rhel8-power9le
module load hypre/2.31.0-gcc-10.2.0-3ijm4v7
# parmetis@=4.0.3%gcc@=10.2.0~gdb~int64~ipo+shared build_system=cmake build_type=Release generator=make patches=4f89253,50ed208,704b84f arch=linux-rhel8-power9le
module load parmetis/4.0.3-gcc-10.2.0-ngzr2hq
# python@=3.8.10%gcc@=10.2.0+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=0d98e93,4c24573,ebdca64,f2fd060 arch=linux-rhel8-power9le
module load python/3.8.10-gcc-10.2.0-ajrbpaz
# superlu-dist@=8.2.1%gcc@=10.2.0~cuda~int64~ipo~openmp+parmetis~rocm+shared build_system=cmake build_type=Release generator=make arch=linux-rhel8-power9le
module load superlu-dist/8.2.1-gcc-10.2.0-da3dc3r
# petsc@=3.21.5%gcc@=10.2.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib+hdf5~hpddm~hwloc+hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind+metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse+superlu-dist~sycl~tetgen~trilinos~valgrind~zoltan build_system=generic clanguage=C memalign=none arch=linux-rhel8-power9le
module load petsc/3.21.5-gcc-10.2.0-cqs7geh
# python-venv@=1.0%gcc@=10.2.0 build_system=generic arch=linux-rhel8-power9le
module load python-venv/1.0-gcc-10.2.0-luagd3n
# py-pip@=23.1.2%gcc@=10.2.0 build_system=generic arch=linux-rhel8-power9le
module load py-pip/23.1.2-gcc-10.2.0-h6vgzq4
# py-setuptools@=69.2.0%gcc@=10.2.0 build_system=generic arch=linux-rhel8-power9le
module load py-setuptools/69.2.0-gcc-10.2.0-rohuias
# py-wheel@=0.41.2%gcc@=10.2.0 build_system=generic arch=linux-rhel8-power9le
module load py-wheel/0.41.2-gcc-10.2.0-u3zzpw7
# py-cython@=0.29.36%gcc@=10.2.0 build_system=python_pip patches=c4369ad arch=linux-rhel8-power9le
module load py-cython/0.29.36-gcc-10.2.0-ccntxzk
# py-mpi4py@=3.1.6%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-mpi4py/3.1.6-gcc-10.2.0-ca4xqhn
# py-flit-core@=3.9.0%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-flit-core/3.9.0-gcc-10.2.0-n5aiq2z
# git@=2.31.1%gcc@=10.2.0+man+nls+perl+subtree~svn~tcltk build_system=autotools arch=linux-rhel8-power9le
module load git/2.31.1-gcc-10.2.0-ghacwlt
# py-packaging@=23.1%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-packaging/23.1-gcc-10.2.0-swqk3sy
# py-tomli@=2.0.1%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-tomli/2.0.1-gcc-10.2.0-rz6scxp
# py-typing-extensions@=4.8.0%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-typing-extensions/4.8.0-gcc-10.2.0-vxrgpiv
# py-setuptools-scm@=8.0.4%gcc@=10.2.0+toml build_system=python_pip arch=linux-rhel8-power9le
module load py-setuptools-scm/8.0.4-gcc-10.2.0-maruw4n
# py-flit-scm@=1.7.0%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-flit-scm/1.7.0-gcc-10.2.0-jczguhi
# py-exceptiongroup@=1.1.1%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-exceptiongroup/1.1.1-gcc-10.2.0-tfpftyl
# py-editables@=0.5%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-editables/0.5-gcc-10.2.0-ygrqbnu
# py-pathspec@=0.11.1%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-pathspec/0.11.1-gcc-10.2.0-pb6lexf
# py-pluggy@=1.5.0%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-pluggy/1.5.0-gcc-10.2.0-tww7bn7
# py-calver@=2022.6.26%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-calver/2022.6.26-gcc-10.2.0-j4n6pgj
# py-trove-classifiers@=2023.8.7%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-trove-classifiers/2023.8.7-gcc-10.2.0-4a2tyf4
# py-hatchling@=1.21.0%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-hatchling/1.21.0-gcc-10.2.0-7w3bstn
# py-hatch-vcs@=0.3.0%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-hatch-vcs/0.3.0-gcc-10.2.0-habcfns
# py-iniconfig@=2.0.0%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-iniconfig/2.0.0-gcc-10.2.0-7cmoqd3
# py-pytest@=8.2.1%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-pytest/8.2.1-gcc-10.2.0-p4ovyyp
# exago@=develop%gcc@=10.2.0+cuda+hiop~ipo+ipopt+logging+mpi+python+raja~rocm build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-rhel8-power9le
## module load exago/develop-gcc-10.2.0-ieiuv7j
