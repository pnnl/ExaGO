module use -a /gpfs/wolf/proj-shared/csc359/exago/spack-ci/install/modules/linux-rhel8-power9le
# cmake@=3.22.2%gcc@=11.2.0~doc+ncurses+ownlibs build_system=generic build_type=Release arch=linux-rhel8-power9le
module load cmake/3.22.2-gcc-11.2.0-c6karkt
# blt@=0.4.1%gcc@=11.2.0 build_system=generic arch=linux-rhel8-power9le
module load blt/0.4.1-gcc-11.2.0-m6zk6cg
# cub@=2.1.0%gcc@=11.2.0 build_system=generic arch=linux-rhel8-power9le
module load cub/2.1.0-gcc-11.2.0-nqylh3l
# gnuconfig@=2022-09-17%gcc@=11.2.0 build_system=generic arch=linux-rhel8-power9le
module load gnuconfig/2022-09-17-gcc-11.2.0-2giv246
# libiconv@=1.17%gcc@=11.2.0 build_system=autotools libs=shared,static arch=linux-rhel8-power9le
module load libiconv/1.17-gcc-11.2.0-banm3b5
# pkgconf@=1.9.5%gcc@=11.2.0 build_system=autotools arch=linux-rhel8-power9le
module load pkgconf/1.9.5-gcc-11.2.0-lvy36r4
# xz@=5.4.1%gcc@=11.2.0~pic build_system=autotools libs=shared,static arch=linux-rhel8-power9le
module load xz/5.4.1-gcc-11.2.0-r2loadq
# zlib-ng@=2.1.3%gcc@=11.2.0+compat+opt build_system=autotools patches=299b958,ae9077a,b692621 arch=linux-rhel8-power9le
module load zlib-ng/2.1.3-gcc-11.2.0-qut5h6m
# libxml2@=2.10.3%gcc@=11.2.0+pic~python+shared build_system=autotools arch=linux-rhel8-power9le
module load libxml2/2.10.3-gcc-11.2.0-pznasx7
# cuda@=11.8.0%gcc@=11.2.0~allow-unsupported-compilers~dev build_system=generic arch=linux-rhel8-power9le
module load cuda/11.8.0-gcc-11.2.0-pjldssb
# gmake@=4.4.1%gcc@=11.2.0~guile build_system=autotools arch=linux-rhel8-power9le
module load gmake/4.4.1-gcc-11.2.0-mo7higu
# camp@=0.2.3%gcc@=11.2.0+cuda~ipo+openmp~rocm~tests build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-rhel8-power9le
module load camp/0.2.3-gcc-11.2.0-w6fb4j6
# perl@=5.30.1%gcc@=11.2.0+cpanm+opcode+open+shared+threads build_system=generic arch=linux-rhel8-power9le
module load perl/5.30.1-gcc-11.2.0-h4oqz4a
# openblas@=0.3.24%gcc@=11.2.0~bignuma~consistent_fpcsr+fortran~ilp64+locking+pic+shared build_system=makefile symbol_suffix=none threads=none arch=linux-rhel8-power9le
module load openblas/0.3.24-gcc-11.2.0-mx4owbc
# coinhsl@=2019.05.21%gcc@=11.2.0+blas build_system=autotools arch=linux-rhel8-power9le
module load coinhsl/2019.05.21-gcc-11.2.0-47yd5po
# ginkgo@=1.5.0.glu_experimental%gcc@=11.2.0+cuda~develtools~full_optimizations~hwloc~ipo~mpi~oneapi+openmp~rocm+shared build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-rhel8-power9le
module load ginkgo/1.5.0.glu_experimental-gcc-11.2.0-e3i6qrb
# magma@=2.6.2%gcc@=11.2.0+cuda+fortran~ipo~rocm+shared build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-rhel8-power9le
module load magma/2.6.2-gcc-11.2.0-tiugqi3
# metis@=5.1.0%gcc@=11.2.0~gdb~int64~ipo~real64+shared build_system=cmake build_type=Release generator=make patches=4991da9,93a7903,b1225da arch=linux-rhel8-power9le
module load metis/5.1.0-gcc-11.2.0-phx5jh2
# raja@=0.14.0%gcc@=11.2.0+cuda~examples~exercises~ipo+openmp~rocm+shared~tests build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-rhel8-power9le
module load raja/0.14.0-gcc-11.2.0-6kcew5j
# spectrum-mpi@=10.4.0.3-20210112%gcc@=11.2.0 build_system=bundle arch=linux-rhel8-power9le
module load spectrum-mpi/10.4.0.3-20210112-gcc-11.2.0-jflmvka
# gmp@=6.2.1%gcc@=11.2.0+cxx build_system=autotools libs=shared,static patches=69ad2e2 arch=linux-rhel8-power9le
module load gmp/6.2.1-gcc-11.2.0-acpul5s
# diffutils@=3.9%gcc@=11.2.0 build_system=autotools arch=linux-rhel8-power9le
module load diffutils/3.9-gcc-11.2.0-dayncqp
# libsigsegv@=2.14%gcc@=11.2.0 build_system=autotools arch=linux-rhel8-power9le
module load libsigsegv/2.14-gcc-11.2.0-ofen7pj
# m4@=1.4.19%gcc@=11.2.0+sigsegv build_system=autotools patches=9dc5fbd,bfdffa7 arch=linux-rhel8-power9le
module load m4/1.4.19-gcc-11.2.0-zi7es42
# autoconf@=2.69%gcc@=11.2.0 build_system=autotools patches=35c4492,7793209,a49dd5b arch=linux-rhel8-power9le
module load autoconf/2.69-gcc-11.2.0-unnjts5
# autoconf-archive@=2023.02.20%gcc@=11.2.0 build_system=autotools arch=linux-rhel8-power9le
module load autoconf-archive/2023.02.20-gcc-11.2.0-3eft6su
# automake@=1.16.5%gcc@=11.2.0 build_system=autotools arch=linux-rhel8-power9le
module load automake/1.16.5-gcc-11.2.0-zup3cph
# libtool@=2.4.7%gcc@=11.2.0 build_system=autotools arch=linux-rhel8-power9le
module load libtool/2.4.7-gcc-11.2.0-g2c3rnj
# texinfo@=6.5%gcc@=11.2.0 build_system=autotools patches=12f6edb,1732115 arch=linux-rhel8-power9le
module load texinfo/6.5-gcc-11.2.0-r37p4yk
# mpfr@=4.2.0%gcc@=11.2.0 build_system=autotools libs=shared,static arch=linux-rhel8-power9le
module load mpfr/4.2.0-gcc-11.2.0-kom6hev
# suite-sparse@=5.13.0%gcc@=11.2.0~cuda~graphblas~openmp+pic build_system=generic arch=linux-rhel8-power9le
module load suite-sparse/5.13.0-gcc-11.2.0-nxauhvy
# umpire@=6.0.0%gcc@=11.2.0+c+cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp~rocm~shared build_system=cmake build_type=Release cuda_arch=70 generator=make tests=none arch=linux-rhel8-power9le
module load umpire/6.0.0-gcc-11.2.0-bg3pjfg
# hiop@=1.0.0%gcc@=11.2.0+cuda+cusolver_lu~deepchecking+ginkgo~ipo~jsrun+kron+mpi+raja~rocm~shared+sparse build_system=cmake build_type=MinSizeRel cuda_arch=70 generator=make arch=linux-rhel8-power9le
module load hiop/1.0.0-gcc-11.2.0-fc36rfn
# ipopt@=3.12.10%gcc@=11.2.0+coinhsl~debug~metis~mumps build_system=autotools arch=linux-rhel8-power9le
module load ipopt/3.12.10-gcc-11.2.0-7fkxunc
# bzip2@=1.0.8%gcc@=11.2.0~debug~pic+shared build_system=generic arch=linux-rhel8-power9le
module load bzip2/1.0.8-gcc-11.2.0-mkrdkja
# libmd@=1.0.4%gcc@=11.2.0 build_system=autotools arch=linux-rhel8-power9le
module load libmd/1.0.4-gcc-11.2.0-7fmyd5e
# libbsd@=0.11.7%gcc@=11.2.0 build_system=autotools arch=linux-rhel8-power9le
module load libbsd/0.11.7-gcc-11.2.0-l7yvt7d
# expat@=2.5.0%gcc@=11.2.0+libbsd build_system=autotools arch=linux-rhel8-power9le
module load expat/2.5.0-gcc-11.2.0-v4j5xwo
# ncurses@=6.4%gcc@=11.2.0~symlinks+termlib abi=none build_system=autotools arch=linux-rhel8-power9le
module load ncurses/6.4-gcc-11.2.0-zer2yt3
# readline@=8.2%gcc@=11.2.0 build_system=autotools patches=bbf97f1 arch=linux-rhel8-power9le
module load readline/8.2-gcc-11.2.0-r2dcinw
# gdbm@=1.23%gcc@=11.2.0 build_system=autotools arch=linux-rhel8-power9le
module load gdbm/1.23-gcc-11.2.0-qwankyx
# tar@=1.34%gcc@=11.2.0 build_system=autotools zip=pigz arch=linux-rhel8-power9le
module load tar/1.34-gcc-11.2.0-2jrcbem
# gettext@=0.21.1%gcc@=11.2.0+bzip2+curses+git~libunistring+libxml2+tar+xz build_system=autotools arch=linux-rhel8-power9le
module load gettext/0.21.1-gcc-11.2.0-hq5nrhb
# libffi@=3.4.4%gcc@=11.2.0 build_system=autotools arch=linux-rhel8-power9le
module load libffi/3.4.4-gcc-11.2.0-n4x452n
# libxcrypt@=4.4.35%gcc@=11.2.0~obsolete_api build_system=autotools patches=4885da3 arch=linux-rhel8-power9le
module load libxcrypt/4.4.35-gcc-11.2.0-x7eom22
# ca-certificates-mozilla@=2023-05-30%gcc@=11.2.0 build_system=generic arch=linux-rhel8-power9le
module load ca-certificates-mozilla/2023-05-30-gcc-11.2.0-iekj6zv
# openssl@=3.1.3%gcc@=11.2.0~docs+shared build_system=generic certs=mozilla arch=linux-rhel8-power9le
## module load openssl/3.1.3-gcc-11.2.0-gwvthxv
# sqlite@=3.42.0%gcc@=11.2.0+column_metadata+dynamic_extensions+fts~functions+rtree build_system=autotools arch=linux-rhel8-power9le
module load sqlite/3.42.0-gcc-11.2.0-ycvjqs4
# util-linux-uuid@=2.38.1%gcc@=11.2.0 build_system=autotools arch=linux-rhel8-power9le
module load util-linux-uuid/2.38.1-gcc-11.2.0-5tjc54d
# python@=3.11.4%gcc@=11.2.0+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=13fa8bf,b0615b2,ebdca64,f2fd060 arch=linux-rhel8-power9le
module load python/3.11.4-gcc-11.2.0-ttme25e
# petsc@=3.19.6%gcc@=11.2.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~tetgen~trilinos~valgrind build_system=generic clanguage=C memalign=none arch=linux-rhel8-power9le
module load petsc/3.19.6-gcc-11.2.0-eriyan4
# py-pip@=23.1.2%gcc@=11.2.0 build_system=generic arch=linux-rhel8-power9le
module load py-pip/23.1.2-gcc-11.2.0-ic7sml4
# py-setuptools@=68.0.0%gcc@=11.2.0 build_system=generic arch=linux-rhel8-power9le
module load py-setuptools/68.0.0-gcc-11.2.0-bvs3vhu
# py-wheel@=0.37.1%gcc@=11.2.0 build_system=generic arch=linux-rhel8-power9le
module load py-wheel/0.37.1-gcc-11.2.0-k4awc4n
# py-cython@=0.29.36%gcc@=11.2.0 build_system=python_pip patches=c4369ad arch=linux-rhel8-power9le
module load py-cython/0.29.36-gcc-11.2.0-ukkzjdw
# py-mpi4py@=3.1.4%gcc@=11.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-mpi4py/3.1.4-gcc-11.2.0-vjak4q5
# py-editables@=0.3%gcc@=11.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-editables/0.3-gcc-11.2.0-64km524
# py-flit-core@=3.9.0%gcc@=11.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-flit-core/3.9.0-gcc-11.2.0-i75hw7u
# py-packaging@=23.1%gcc@=11.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-packaging/23.1-gcc-11.2.0-gwwtxns
# py-pathspec@=0.11.1%gcc@=11.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-pathspec/0.11.1-gcc-11.2.0-j3cvp7f
# git@=2.35.1%gcc@=11.2.0+man+nls+perl+subtree~svn~tcltk build_system=autotools arch=linux-rhel8-power9le
module load git/2.35.1-gcc-11.2.0-vpvo7mj
# py-tomli@=2.0.1%gcc@=11.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-tomli/2.0.1-gcc-11.2.0-j3ewsfv
# py-typing-extensions@=4.8.0%gcc@=11.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-typing-extensions/4.8.0-gcc-11.2.0-2ptsvzd
# py-setuptools-scm@=7.1.0%gcc@=11.2.0+toml build_system=python_pip arch=linux-rhel8-power9le
module load py-setuptools-scm/7.1.0-gcc-11.2.0-naetuvu
# py-pluggy@=1.0.0%gcc@=11.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-pluggy/1.0.0-gcc-11.2.0-z4rnrss
# py-calver@=2022.6.26%gcc@=11.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-calver/2022.6.26-gcc-11.2.0-c7wrei6
# py-trove-classifiers@=2023.8.7%gcc@=11.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-trove-classifiers/2023.8.7-gcc-11.2.0-gig62u3
# py-hatchling@=1.18.0%gcc@=11.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-hatchling/1.18.0-gcc-11.2.0-cukj4b5
# py-hatch-vcs@=0.3.0%gcc@=11.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-hatch-vcs/0.3.0-gcc-11.2.0-oy5cwb2
# py-iniconfig@=2.0.0%gcc@=11.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-iniconfig/2.0.0-gcc-11.2.0-zlv4gcn
# py-pytest@=7.3.2%gcc@=11.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-pytest/7.3.2-gcc-11.2.0-ctedzyz
# exago@=develop%gcc@=11.2.0+cuda+hiop~ipo+ipopt~logging+mpi+python+raja~rocm build_system=cmake build_type=Debug cuda_arch=70 dev_path=/gpfs/wolf/proj-shared/csc359/ci/467269 generator=make arch=linux-rhel8-power9le
## module load exago/develop-gcc-11.2.0-ycftvjb
