module use -a /qfs/projects/exasgd/src/deception-ci/install/ci-modules/linux-centos7-zen2
# cmake@=3.26.3%gcc@=9.1.0~doc+ncurses+ownlibs build_system=generic build_type=Release arch=linux-centos7-zen2
module load cmake/3.26.3-gcc-9.1.0-wryjrgo
# blt@=0.4.1%gcc@=9.1.0 build_system=generic arch=linux-centos7-zen2
module load blt/0.4.1-gcc-9.1.0-2kcdf2y
# cub@=2.1.0%gcc@=9.1.0 build_system=generic arch=linux-centos7-zen2
module load cub/2.1.0-gcc-9.1.0-wx4parb
# cuda@=11.4%gcc@=9.1.0~allow-unsupported-compilers~dev build_system=generic arch=linux-centos7-zen2
module load cuda/11.4-gcc-9.1.0-tmrg4gr
# gmake@=4.4.1%gcc@=9.1.0~guile build_system=generic arch=linux-centos7-zen2
module load gmake/4.4.1-gcc-9.1.0-2hyua4w
# camp@=0.2.3%gcc@=9.1.0+cuda~ipo+openmp~rocm~tests build_system=cmake build_type=Release cuda_arch=60,70,75,80 generator=make arch=linux-centos7-zen2
module load camp/0.2.3-gcc-9.1.0-hyzpmsn
# perl@=5.26.0%gcc@=9.1.0+cpanm+opcode+open+shared+threads build_system=generic patches=0eac10e,8cf4302 arch=linux-centos7-zen2
module load perl/5.26.0-gcc-9.1.0-aca44by
# openblas@=0.3.25%gcc@=9.1.0~bignuma~consistent_fpcsr+fortran~ilp64+locking+pic+shared build_system=makefile symbol_suffix=none threads=none arch=linux-centos7-zen2
module load openblas/0.3.25-gcc-9.1.0-edjd2fv
# coinhsl@=2019.05.21%gcc@=9.1.0+blas build_system=autotools arch=linux-centos7-zen2
module load coinhsl/2019.05.21-gcc-9.1.0-xttde3f
# ginkgo@=1.5.0.glu_experimental%gcc@=9.1.0+cuda~develtools~full_optimizations~hwloc~ipo~mpi+openmp~rocm+shared~sycl build_system=cmake build_type=Release cuda_arch=60,70,75,80 generator=make arch=linux-centos7-zen2
module load ginkgo/1.5.0.glu_experimental-gcc-9.1.0-byoaioj
# magma@=2.6.2%gcc@=9.1.0+cuda+fortran~ipo~rocm+shared build_system=cmake build_type=Release cuda_arch=60,70,75,80 generator=make arch=linux-centos7-zen2
module load magma/2.6.2-gcc-9.1.0-eybt5h4
# metis@=5.1.0%gcc@=9.1.0~gdb~int64~ipo~real64+shared build_system=cmake build_type=Release generator=make patches=4991da9,93a7903,b1225da arch=linux-centos7-zen2
module load metis/5.1.0-gcc-9.1.0-yjhcoa3
# openmpi@=4.1.0mlx5.0%gcc@=9.1.0~atomics~cuda~cxx~cxx_exceptions~gpfs~internal-hwloc~internal-pmix~java~legacylaunchers~lustre~memchecker~openshmem~orterunprefix+romio+rsh~singularity+static+vt+wrapper-rpath build_system=autotools fabrics=none patches=60ce20b schedulers=none arch=linux-centos7-zen2
module load openmpi/4.1.0mlx5.0-gcc-9.1.0-vj5ufod
# raja@=0.14.0%gcc@=9.1.0+cuda~examples~exercises~ipo+openmp~plugins~rocm+shared~tests build_system=cmake build_type=Release cuda_arch=60,70,75,80 generator=make arch=linux-centos7-zen2
module load raja/0.14.0-gcc-9.1.0-tq2hxbw
# libiconv@=1.17%gcc@=9.1.0 build_system=autotools libs=shared,static arch=linux-centos7-zen2
module load libiconv/1.17-gcc-9.1.0-ryrjbul
# diffutils@=3.9%gcc@=9.1.0 build_system=autotools arch=linux-centos7-zen2
module load diffutils/3.9-gcc-9.1.0-4r6gzdy
# libsigsegv@=2.14%gcc@=9.1.0 build_system=autotools arch=linux-centos7-zen2
module load libsigsegv/2.14-gcc-9.1.0-cy447uk
# m4@=1.4.19%gcc@=9.1.0+sigsegv build_system=autotools patches=9dc5fbd,bfdffa7 arch=linux-centos7-zen2
module load m4/1.4.19-gcc-9.1.0-yf2slnz
# autoconf@=2.69%gcc@=9.1.0 build_system=autotools patches=35c4492,7793209,a49dd5b arch=linux-centos7-zen2
module load autoconf/2.69-gcc-9.1.0-fj76bwy
# automake@=1.16.5%gcc@=9.1.0 build_system=autotools arch=linux-centos7-zen2
module load automake/1.16.5-gcc-9.1.0-6ciswaj
# libtool@=2.4.7%gcc@=9.1.0 build_system=autotools arch=linux-centos7-zen2
module load libtool/2.4.7-gcc-9.1.0-tsf7hsq
# gmp@=6.2.1%gcc@=9.1.0+cxx build_system=autotools libs=shared,static patches=69ad2e2 arch=linux-centos7-zen2
module load gmp/6.2.1-gcc-9.1.0-tsjd6bm
# autoconf-archive@=2023.02.20%gcc@=9.1.0 build_system=autotools arch=linux-centos7-zen2
module load autoconf-archive/2023.02.20-gcc-9.1.0-ebkhgoh
# bzip2@=1.0.8%gcc@=9.1.0~debug~pic+shared build_system=generic arch=linux-centos7-zen2
module load bzip2/1.0.8-gcc-9.1.0-mfp2wwg
# pkgconf@=1.9.5%gcc@=9.1.0 build_system=autotools arch=linux-centos7-zen2
module load pkgconf/1.9.5-gcc-9.1.0-bg6z6dj
# xz@=5.4.1%gcc@=9.1.0~pic build_system=autotools libs=shared,static arch=linux-centos7-zen2
module load xz/5.4.1-gcc-9.1.0-cflgxw5
# zlib@=1.3%gcc@=9.1.0+optimize+pic+shared build_system=makefile arch=linux-centos7-zen2
module load zlib/1.3-gcc-9.1.0-3ckc3rk
# libxml2@=2.10.3%gcc@=9.1.0+pic~python+shared build_system=autotools arch=linux-centos7-zen2
module load libxml2/2.10.3-gcc-9.1.0-kel6jjt
# ncurses@=6.4%gcc@=9.1.0~symlinks+termlib abi=none build_system=autotools arch=linux-centos7-zen2
module load ncurses/6.4-gcc-9.1.0-aurkwew
# pigz@=2.7%gcc@=9.1.0 build_system=makefile arch=linux-centos7-zen2
module load pigz/2.7-gcc-9.1.0-pvzit7o
# zstd@=1.5.5%gcc@=9.1.0+programs build_system=makefile compression=none libs=shared,static arch=linux-centos7-zen2
module load zstd/1.5.5-gcc-9.1.0-lcexipe
# tar@=1.34%gcc@=9.1.0 build_system=autotools zip=pigz arch=linux-centos7-zen2
module load tar/1.34-gcc-9.1.0-vpzflun
# gettext@=0.22.4%gcc@=9.1.0+bzip2+curses+git~libunistring+libxml2+pic+shared+tar+xz build_system=autotools arch=linux-centos7-zen2
module load gettext/0.22.4-gcc-9.1.0-6bwgojr
# texinfo@=7.0.3%gcc@=9.1.0 build_system=autotools arch=linux-centos7-zen2
module load texinfo/7.0.3-gcc-9.1.0-kwmj254
# mpfr@=4.2.0%gcc@=9.1.0 build_system=autotools libs=shared,static arch=linux-centos7-zen2
module load mpfr/4.2.0-gcc-9.1.0-2d5o4yj
# suite-sparse@=5.13.0%gcc@=9.1.0~cuda~graphblas~openmp+pic build_system=generic arch=linux-centos7-zen2
module load suite-sparse/5.13.0-gcc-9.1.0-aytsg72
# umpire@=6.0.0%gcc@=9.1.0+c+cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp~rocm~shared build_system=cmake build_type=Release cuda_arch=60,70,75,80 generator=make tests=none arch=linux-centos7-zen2
module load umpire/6.0.0-gcc-9.1.0-zh3kabb
# hiop@=develop%gcc@=9.1.0+cuda~cusolver_lu~deepchecking+ginkgo~ipo~jsrun+kron+mpi+raja~rocm~shared+sparse build_system=cmake build_type=MinSizeRel cuda_arch=60,70,75,80 dev_path=/people/svcexasgd/gitlab/26233/spack_deception/hiop_dev generator=make arch=linux-centos7-zen2
module load hiop/develop-gcc-9.1.0-kllcazs
# ipopt@=3.12.10%gcc@=9.1.0+coinhsl~debug~metis~mumps build_system=autotools arch=linux-centos7-zen2
module load ipopt/3.12.10-gcc-9.1.0-li7gpqb
# python@=3.9.12%gcc@=9.1.0+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=0d98e93,4c24573,ebdca64,f2fd060 arch=linux-centos7-zen2
module load python/3.9.12-gcc-9.1.0-a5w6cuq
# petsc@=3.20.1%gcc@=9.1.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~sycl~tetgen~trilinos~valgrind build_system=generic clanguage=C memalign=none arch=linux-centos7-zen2
module load petsc/3.20.1-gcc-9.1.0-llz33ok
# py-pip@=23.1.2%gcc@=9.1.0 build_system=generic arch=linux-centos7-zen2
module load py-pip/23.1.2-gcc-9.1.0-pvy73ms
# py-setuptools@=68.0.0%gcc@=9.1.0 build_system=generic arch=linux-centos7-zen2
module load py-setuptools/68.0.0-gcc-9.1.0-fef5ebz
# py-wheel@=0.41.2%gcc@=9.1.0 build_system=generic arch=linux-centos7-zen2
module load py-wheel/0.41.2-gcc-9.1.0-ycwjnxk
# py-cython@=0.29.36%gcc@=9.1.0 build_system=python_pip patches=c4369ad arch=linux-centos7-zen2
module load py-cython/0.29.36-gcc-9.1.0-lhxhitm
# py-mpi4py@=3.1.4%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-mpi4py/3.1.4-gcc-9.1.0-y5dj2ru
# py-flit-core@=3.9.0%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-flit-core/3.9.0-gcc-9.1.0-tgih3ig
# git@=2.37.3%gcc@=9.1.0+man+nls+perl+subtree~svn~tcltk build_system=autotools arch=linux-centos7-zen2
module load git/2.37.3-gcc-9.1.0-qxn54jq
# py-packaging@=23.1%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-packaging/23.1-gcc-9.1.0-3uczz47
# py-tomli@=2.0.1%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-tomli/2.0.1-gcc-9.1.0-6w6ssdm
# py-typing-extensions@=4.8.0%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-typing-extensions/4.8.0-gcc-9.1.0-iv3gy6a
# py-setuptools-scm@=7.1.0%gcc@=9.1.0+toml build_system=python_pip arch=linux-centos7-zen2
module load py-setuptools-scm/7.1.0-gcc-9.1.0-huz56hs
# py-flit-scm@=1.7.0%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-flit-scm/1.7.0-gcc-9.1.0-ihypnl3
# py-exceptiongroup@=1.1.1%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-exceptiongroup/1.1.1-gcc-9.1.0-2uuhbhc
# py-editables@=0.3%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-editables/0.3-gcc-9.1.0-6c6asxe
# py-pathspec@=0.11.1%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-pathspec/0.11.1-gcc-9.1.0-tja5anm
# py-pluggy@=1.0.0%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-pluggy/1.0.0-gcc-9.1.0-gmrzy4m
# py-calver@=2022.6.26%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-calver/2022.6.26-gcc-9.1.0-6pb422e
# py-trove-classifiers@=2023.8.7%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-trove-classifiers/2023.8.7-gcc-9.1.0-wvlyozv
# py-hatchling@=1.18.0%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-hatchling/1.18.0-gcc-9.1.0-qpw4or6
# py-hatch-vcs@=0.3.0%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-hatch-vcs/0.3.0-gcc-9.1.0-hwwyqic
# py-iniconfig@=2.0.0%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-iniconfig/2.0.0-gcc-9.1.0-raa2jbg
# py-pytest@=7.3.2%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-pytest/7.3.2-gcc-9.1.0-2nh22mh
# exago@=develop%gcc@=9.1.0+cuda+hiop~ipo+ipopt+logging+mpi+python+raja~rocm build_system=cmake build_type=Debug cuda_arch=60,70,75,80 dev_path=/people/svcexasgd/gitlab/26233/spack_deception generator=make arch=linux-centos7-zen2
## module load exago/develop-gcc-9.1.0-27tyuch
