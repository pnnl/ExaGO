module use -a /qfs/projects/earthshot/src/deception-ci/install/ci-modules/linux-centos7-zen2
# cmake@=3.26.3%gcc@=9.1.0~doc+ncurses+ownlibs build_system=generic build_type=Release arch=linux-centos7-zen2
module load cmake/3.26.3-gcc-9.1.0-wryjrgo
# gcc-runtime@=9.1.0%gcc@=9.1.0 build_system=generic arch=linux-centos7-zen2
module load gcc-runtime/9.1.0-gcc-9.1.0-abiv3jd
# blt@=0.4.1%gcc@=9.1.0 build_system=generic arch=linux-centos7-zen2
module load blt/0.4.1-gcc-9.1.0-r65qngy
# cub@=2.1.0%gcc@=9.1.0 build_system=generic arch=linux-centos7-zen2
module load cub/2.1.0-gcc-9.1.0-toinftv
# cuda@=11.4%gcc@=9.1.0~allow-unsupported-compilers~dev build_system=generic arch=linux-centos7-zen2
module load cuda/11.4-gcc-9.1.0-tmrg4gr
# gmake@=4.4.1%gcc@=9.1.0~guile build_system=generic arch=linux-centos7-zen2
module load gmake/4.4.1-gcc-9.1.0-zro7edd
# camp@=0.2.3%gcc@=9.1.0+cuda~ipo+openmp~rocm~tests build_system=cmake build_type=Release cuda_arch=60,70,75,80 generator=make patches=cb9e25b arch=linux-centos7-zen2
module load camp/0.2.3-gcc-9.1.0-wrcymwz
# ginkgo@=1.5.0.glu_experimental%gcc@=9.1.0+cuda~develtools~full_optimizations~hwloc~ipo~mpi+openmp~rocm+shared~sycl build_system=cmake build_type=Release cuda_arch=60,70,75,80 generator=make arch=linux-centos7-zen2
module load ginkgo/1.5.0.glu_experimental-gcc-9.1.0-uxadao2
# perl@=5.26.0%gcc@=9.1.0+cpanm+opcode+open+shared+threads build_system=generic patches=0eac10e,8cf4302 arch=linux-centos7-zen2
module load perl/5.26.0-gcc-9.1.0-cjrkygi
# openblas@=0.3.26%gcc@=9.1.0~bignuma~consistent_fpcsr+fortran~ilp64+locking+pic+shared build_system=makefile symbol_suffix=none threads=none arch=linux-centos7-zen2
module load openblas/0.3.26-gcc-9.1.0-ksf53ii
# coinhsl@=2019.05.21%gcc@=9.1.0+blas build_system=autotools arch=linux-centos7-zen2
module load coinhsl/2019.05.21-gcc-9.1.0-gpha2c7
# magma@=2.6.2%gcc@=9.1.0+cuda+fortran~ipo~rocm+shared build_system=cmake build_type=Release cuda_arch=60,70,75,80 generator=make arch=linux-centos7-zen2
module load magma/2.6.2-gcc-9.1.0-5sbhjkx
# metis@=5.1.0%gcc@=9.1.0~gdb~int64~ipo~real64+shared build_system=cmake build_type=Release generator=make patches=4991da9,93a7903,b1225da arch=linux-centos7-zen2
module load metis/5.1.0-gcc-9.1.0-f7o2p3m
# openmpi@=4.1.0mlx5.0%gcc@=9.1.0~atomics~cuda~cxx~cxx_exceptions~gpfs~internal-hwloc~internal-libevent~internal-pmix~java~legacylaunchers~lustre~memchecker~openshmem~orterunprefix~romio+rsh~singularity~static+vt+wrapper-rpath build_system=autotools fabrics=none patches=60ce20b schedulers=none arch=linux-centos7-zen2
module load openmpi/4.1.0mlx5.0-gcc-9.1.0-wxzawxc
# raja@=0.14.0%gcc@=9.1.0+cuda~examples~exercises~ipo+openmp~plugins~rocm+shared~tests build_system=cmake build_type=Release cuda_arch=60,70,75,80 generator=make arch=linux-centos7-zen2
module load raja/0.14.0-gcc-9.1.0-poskbuz
# libiconv@=1.17%gcc@=9.1.0 build_system=autotools libs=shared,static arch=linux-centos7-zen2
module load libiconv/1.17-gcc-9.1.0-2st3zf2
# diffutils@=3.10%gcc@=9.1.0 build_system=autotools arch=linux-centos7-zen2
module load diffutils/3.10-gcc-9.1.0-ackf5zf
# libsigsegv@=2.14%gcc@=9.1.0 build_system=autotools arch=linux-centos7-zen2
module load libsigsegv/2.14-gcc-9.1.0-6jgesu5
# m4@=1.4.19%gcc@=9.1.0+sigsegv build_system=autotools patches=9dc5fbd,bfdffa7 arch=linux-centos7-zen2
module load m4/1.4.19-gcc-9.1.0-upgj45g
# autoconf@=2.72%gcc@=9.1.0 build_system=autotools arch=linux-centos7-zen2
module load autoconf/2.72-gcc-9.1.0-f7ut5ij
# automake@=1.16.5%gcc@=9.1.0 build_system=autotools arch=linux-centos7-zen2
module load automake/1.16.5-gcc-9.1.0-rt4hcgd
# findutils@=4.9.0%gcc@=9.1.0 build_system=autotools patches=440b954 arch=linux-centos7-zen2
module load findutils/4.9.0-gcc-9.1.0-265efmp
# libtool@=2.4.7%gcc@=9.1.0 build_system=autotools arch=linux-centos7-zen2
module load libtool/2.4.7-gcc-9.1.0-7gpcxtf
# gmp@=6.2.1%gcc@=9.1.0+cxx build_system=autotools libs=shared,static patches=69ad2e2 arch=linux-centos7-zen2
module load gmp/6.2.1-gcc-9.1.0-q5u7hor
# autoconf-archive@=2023.02.20%gcc@=9.1.0 build_system=autotools arch=linux-centos7-zen2
module load autoconf-archive/2023.02.20-gcc-9.1.0-v26a2vn
# bzip2@=1.0.8%gcc@=9.1.0~debug~pic+shared build_system=generic arch=linux-centos7-zen2
module load bzip2/1.0.8-gcc-9.1.0-ckcsq7c
# pkgconf@=1.9.5%gcc@=9.1.0 build_system=autotools arch=linux-centos7-zen2
module load pkgconf/1.9.5-gcc-9.1.0-zmbhfgn
# xz@=5.4.6%gcc@=9.1.0~pic build_system=autotools libs=shared,static arch=linux-centos7-zen2
module load xz/5.4.6-gcc-9.1.0-nqmk7cs
# zlib@=1.3%gcc@=9.1.0+optimize+pic+shared build_system=makefile arch=linux-centos7-zen2
module load zlib/1.3-gcc-9.1.0-glaosn3
# libxml2@=2.10.3%gcc@=9.1.0+pic~python+shared build_system=autotools arch=linux-centos7-zen2
module load libxml2/2.10.3-gcc-9.1.0-2tijqp2
# ncurses@=6.4%gcc@=9.1.0~symlinks+termlib abi=none build_system=autotools arch=linux-centos7-zen2
module load ncurses/6.4-gcc-9.1.0-ue2vvf3
# pigz@=2.8%gcc@=9.1.0 build_system=makefile arch=linux-centos7-zen2
module load pigz/2.8-gcc-9.1.0-mvwpico
# zstd@=1.5.5%gcc@=9.1.0+programs build_system=makefile compression=none libs=shared,static arch=linux-centos7-zen2
module load zstd/1.5.5-gcc-9.1.0-twxts5j
# tar@=1.34%gcc@=9.1.0 build_system=autotools zip=pigz arch=linux-centos7-zen2
module load tar/1.34-gcc-9.1.0-yivbsi7
# gettext@=0.22.4%gcc@=9.1.0+bzip2+curses+git~libunistring+libxml2+pic+shared+tar+xz build_system=autotools arch=linux-centos7-zen2
module load gettext/0.22.4-gcc-9.1.0-ivckr74
# texinfo@=7.0.3%gcc@=9.1.0 build_system=autotools arch=linux-centos7-zen2
module load texinfo/7.0.3-gcc-9.1.0-phuxsxp
# mpfr@=4.2.1%gcc@=9.1.0 build_system=autotools libs=shared,static arch=linux-centos7-zen2
module load mpfr/4.2.1-gcc-9.1.0-gbnsni6
# suite-sparse@=5.13.0%gcc@=9.1.0~cuda~graphblas~openmp+pic build_system=generic arch=linux-centos7-zen2
module load suite-sparse/5.13.0-gcc-9.1.0-jwgbvrn
# umpire@=6.0.0%gcc@=9.1.0+c+cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp~rocm~shared build_system=cmake build_type=Release cuda_arch=60,70,75,80 generator=make tests=none arch=linux-centos7-zen2
module load umpire/6.0.0-gcc-9.1.0-o2367nu
# hiop@=develop%gcc@=9.1.0+cuda~cusolver_lu~deepchecking+ginkgo~ipo~jsrun+kron+mpi+raja~rocm~shared+sparse build_system=cmake build_type=Release cuda_arch=60,70,75,80 dev_path=/people/svcearthshot/gitlab/31746/spack_deception/hiop_dev generator=make arch=linux-centos7-zen2
module load hiop/develop-gcc-9.1.0-xa72s2d
# ipopt@=3.12.10%gcc@=9.1.0+coinhsl~debug~metis~mumps build_system=autotools arch=linux-centos7-zen2
module load ipopt/3.12.10-gcc-9.1.0-zdw2vxv
# python@=3.9.12%gcc@=9.1.0+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=0d98e93,4c24573,ebdca64,f2fd060 arch=linux-centos7-zen2
module load python/3.9.12-gcc-9.1.0-xlhqfig
# petsc@=3.20.4%gcc@=9.1.0~X+batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~sycl~tetgen~trilinos~valgrind~zoltan build_system=generic clanguage=C memalign=none arch=linux-centos7-zen2
module load petsc/3.20.4-gcc-9.1.0-qvnl7v7
# py-pip@=23.1.2%gcc@=9.1.0 build_system=generic arch=linux-centos7-zen2
module load py-pip/23.1.2-gcc-9.1.0-nrjrbog
# py-setuptools@=69.2.0%gcc@=9.1.0 build_system=generic arch=linux-centos7-zen2
module load py-setuptools/69.2.0-gcc-9.1.0-wruumo7
# py-wheel@=0.41.2%gcc@=9.1.0 build_system=generic arch=linux-centos7-zen2
module load py-wheel/0.41.2-gcc-9.1.0-ouezqxe
# py-cython@=0.29.36%gcc@=9.1.0 build_system=python_pip patches=c4369ad arch=linux-centos7-zen2
module load py-cython/0.29.36-gcc-9.1.0-v4ukl22
# py-mpi4py@=3.1.5%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-mpi4py/3.1.5-gcc-9.1.0-kuiww6g
# py-flit-core@=3.9.0%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-flit-core/3.9.0-gcc-9.1.0-2dpwfr7
# git@=2.37.3%gcc@=9.1.0+man+nls+perl+subtree~svn~tcltk build_system=autotools arch=linux-centos7-zen2
module load git/2.37.3-gcc-9.1.0-qxn54jq
# py-packaging@=23.1%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-packaging/23.1-gcc-9.1.0-mnl5iwv
# py-tomli@=2.0.1%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-tomli/2.0.1-gcc-9.1.0-24mbqlo
# py-typing-extensions@=4.8.0%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-typing-extensions/4.8.0-gcc-9.1.0-2lhchlw
# py-setuptools-scm@=7.1.0%gcc@=9.1.0+toml build_system=python_pip arch=linux-centos7-zen2
module load py-setuptools-scm/7.1.0-gcc-9.1.0-jw2iop6
# py-flit-scm@=1.7.0%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-flit-scm/1.7.0-gcc-9.1.0-kbgm7pu
# py-exceptiongroup@=1.1.1%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-exceptiongroup/1.1.1-gcc-9.1.0-qs5eakf
# py-editables@=0.3%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-editables/0.3-gcc-9.1.0-arv3gc3
# py-pathspec@=0.11.1%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-pathspec/0.11.1-gcc-9.1.0-3vaax3z
# py-pluggy@=1.4.0%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-pluggy/1.4.0-gcc-9.1.0-b6gb62e
# py-calver@=2022.6.26%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-calver/2022.6.26-gcc-9.1.0-q26c62c
# py-trove-classifiers@=2023.8.7%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-trove-classifiers/2023.8.7-gcc-9.1.0-wumpenj
# py-hatchling@=1.21.0%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-hatchling/1.21.0-gcc-9.1.0-6rpq2y3
# py-hatch-vcs@=0.3.0%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-hatch-vcs/0.3.0-gcc-9.1.0-tfm37qn
# py-iniconfig@=2.0.0%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-iniconfig/2.0.0-gcc-9.1.0-mtbo5hz
# py-pytest@=8.0.0%gcc@=9.1.0 build_system=python_pip arch=linux-centos7-zen2
module load py-pytest/8.0.0-gcc-9.1.0-c72o7pz
# exago@=develop%gcc@=9.1.0+cuda+hiop~ipo+ipopt+logging+mpi+python+raja~rocm build_system=cmake build_type=Release cuda_arch=60,70,75,80 dev_path=/people/svcearthshot/gitlab/31746/spack_deception generator=make arch=linux-centos7-zen2
## module load exago/develop-gcc-9.1.0-3h2zyl5
