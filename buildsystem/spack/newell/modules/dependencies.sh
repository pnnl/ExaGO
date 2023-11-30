module use -a /qfs/projects/exasgd/src/ci-newll/spack-install/ci-modules/linux-centos8-power9le
# gmake@=4.4.1%gcc@=8.5.0~guile build_system=generic arch=linux-centos8-power9le
module load gmake/4.4.1-gcc-8.5.0-uysj4bc
# gnuconfig@=2022-09-17%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load gnuconfig/2022-09-17-gcc-8.5.0-brd4tz7
# pkgconf@=1.9.5%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load pkgconf/1.9.5-gcc-8.5.0-nrcoqfd
# nghttp2@=1.57.0%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load nghttp2/1.57.0-gcc-8.5.0-pbbbhf7
# ca-certificates-mozilla@=2023-05-30%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load ca-certificates-mozilla/2023-05-30-gcc-8.5.0-ueq22oa
# berkeley-db@=18.1.40%gcc@=8.5.0+cxx~docs+stl build_system=autotools patches=26090f4,b231fcc arch=linux-centos8-power9le
module load berkeley-db/18.1.40-gcc-8.5.0-44mqiak
# libiconv@=1.17%gcc@=8.5.0 build_system=autotools libs=shared,static arch=linux-centos8-power9le
module load libiconv/1.17-gcc-8.5.0-wtzwh5x
# diffutils@=3.9%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load diffutils/3.9-gcc-8.5.0-vl62cx7
# bzip2@=1.0.8%gcc@=8.5.0~debug~pic+shared build_system=generic arch=linux-centos8-power9le
module load bzip2/1.0.8-gcc-8.5.0-rqfpz3l
# ncurses@=6.4%gcc@=8.5.0~symlinks+termlib abi=none build_system=autotools arch=linux-centos8-power9le
module load ncurses/6.4-gcc-8.5.0-t7cgvnx
# readline@=8.2%gcc@=8.5.0 build_system=autotools patches=bbf97f1 arch=linux-centos8-power9le
module load readline/8.2-gcc-8.5.0-ekbdswm
# gdbm@=1.23%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load gdbm/1.23-gcc-8.5.0-bnlfqr7
# zlib-ng@=2.1.4%gcc@=8.5.0+compat+opt build_system=autotools arch=linux-centos8-power9le
module load zlib-ng/2.1.4-gcc-8.5.0-w7n6mdy
# perl@=5.38.0%gcc@=8.5.0+cpanm+opcode+open+shared+threads build_system=generic patches=714e4d1 arch=linux-centos8-power9le
module load perl/5.38.0-gcc-8.5.0-2r4twvm
# openssl@=3.1.3%gcc@=8.5.0~docs+shared build_system=generic certs=mozilla arch=linux-centos8-power9le
## module load openssl/3.1.3-gcc-8.5.0-slvshxg
# curl@=8.4.0%gcc@=8.5.0~gssapi~ldap~libidn2~librtmp~libssh~libssh2+nghttp2 build_system=autotools libs=shared,static tls=openssl arch=linux-centos8-power9le
module load curl/8.4.0-gcc-8.5.0-fj3pe3l
# cmake@=3.27.7%gcc@=8.5.0~doc+ncurses+ownlibs build_system=generic build_type=Release arch=linux-centos8-power9le
module load cmake/3.27.7-gcc-8.5.0-3xxpipf
# blt@=0.4.1%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load blt/0.4.1-gcc-8.5.0-yadj2vt
# cub@=2.1.0%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load cub/2.1.0-gcc-8.5.0-lx5hjs2
# cuda@=11.4%gcc@=8.5.0~allow-unsupported-compilers~dev build_system=generic arch=linux-centos8-power9le
module load cuda/11.4-gcc-8.5.0-kh4c7o4
# camp@=0.2.3%gcc@=8.5.0+cuda~ipo+openmp~rocm~tests build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-centos8-power9le
module load camp/0.2.3-gcc-8.5.0-z42wdnn
# openblas@=0.3.25%gcc@=8.5.0~bignuma~consistent_fpcsr+fortran~ilp64+locking+pic+shared build_system=makefile symbol_suffix=none threads=none arch=linux-centos8-power9le
module load openblas/0.3.25-gcc-8.5.0-wokmlm3
# coinhsl@=2019.05.21%gcc@=8.5.0+blas build_system=autotools arch=linux-centos8-power9le
module load coinhsl/2019.05.21-gcc-8.5.0-z5bfsep
# ginkgo@=1.5.0.glu_experimental%gcc@=8.5.0+cuda~develtools~full_optimizations~hwloc~ipo~mpi+openmp~rocm+shared~sycl build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-centos8-power9le
module load ginkgo/1.5.0.glu_experimental-gcc-8.5.0-to767fj
# magma@=2.6.2%gcc@=8.5.0+cuda+fortran~ipo~rocm+shared build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-centos8-power9le
module load magma/2.6.2-gcc-8.5.0-o7f67se
# metis@=5.1.0%gcc@=8.5.0~gdb~int64~ipo~real64+shared build_system=cmake build_type=Release generator=make patches=4991da9,93a7903,b1225da arch=linux-centos8-power9le
module load metis/5.1.0-gcc-8.5.0-qmqpgih
# openmpi@=4.1.4%gcc@=8.5.0~atomics~cuda~cxx~cxx_exceptions~gpfs~internal-hwloc~internal-pmix~java~legacylaunchers~lustre~memchecker~openshmem~orterunprefix+romio+rsh~singularity+static+vt+wrapper-rpath build_system=autotools fabrics=none schedulers=none arch=linux-centos8-power9le
module load openmpi/4.1.4-gcc-8.5.0-b4va3xy
# raja@=0.14.0%gcc@=8.5.0+cuda~examples~exercises~ipo+openmp~plugins~rocm+shared~tests build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-centos8-power9le
module load raja/0.14.0-gcc-8.5.0-akyygh5
# libsigsegv@=2.14%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load libsigsegv/2.14-gcc-8.5.0-2igjfex
# m4@=1.4.19%gcc@=8.5.0+sigsegv build_system=autotools patches=9dc5fbd,bfdffa7 arch=linux-centos8-power9le
module load m4/1.4.19-gcc-8.5.0-zmwmrjj
# autoconf@=2.69%gcc@=8.5.0 build_system=autotools patches=35c4492,7793209,a49dd5b arch=linux-centos8-power9le
module load autoconf/2.69-gcc-8.5.0-qss6fq3
# automake@=1.16.5%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load automake/1.16.5-gcc-8.5.0-5ak4ma2
# libtool@=2.4.7%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load libtool/2.4.7-gcc-8.5.0-el3zfel
# gmp@=6.2.1%gcc@=8.5.0+cxx build_system=autotools libs=shared,static patches=69ad2e2 arch=linux-centos8-power9le
module load gmp/6.2.1-gcc-8.5.0-shj27b4
# autoconf-archive@=2023.02.20%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load autoconf-archive/2023.02.20-gcc-8.5.0-h3vkl44
# xz@=5.4.1%gcc@=8.5.0~pic build_system=autotools libs=shared,static arch=linux-centos8-power9le
module load xz/5.4.1-gcc-8.5.0-4k2rc3w
# libxml2@=2.10.3%gcc@=8.5.0+pic~python+shared build_system=autotools arch=linux-centos8-power9le
module load libxml2/2.10.3-gcc-8.5.0-wighgm3
# pigz@=2.7%gcc@=8.5.0 build_system=makefile arch=linux-centos8-power9le
module load pigz/2.7-gcc-8.5.0-272vqg3
# zstd@=1.5.5%gcc@=8.5.0+programs build_system=makefile compression=none libs=shared,static arch=linux-centos8-power9le
module load zstd/1.5.5-gcc-8.5.0-p52n5sn
# tar@=1.34%gcc@=8.5.0 build_system=autotools zip=pigz arch=linux-centos8-power9le
module load tar/1.34-gcc-8.5.0-dt37vzs
# gettext@=0.22.4%gcc@=8.5.0+bzip2+curses+git~libunistring+libxml2+pic+shared+tar+xz build_system=autotools arch=linux-centos8-power9le
module load gettext/0.22.4-gcc-8.5.0-abrp5rm
# texinfo@=7.0.3%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load texinfo/7.0.3-gcc-8.5.0-wimw3l7
# mpfr@=4.2.0%gcc@=8.5.0 build_system=autotools libs=shared,static arch=linux-centos8-power9le
module load mpfr/4.2.0-gcc-8.5.0-zts65hz
# suite-sparse@=5.13.0%gcc@=8.5.0~cuda~graphblas~openmp+pic build_system=generic arch=linux-centos8-power9le
module load suite-sparse/5.13.0-gcc-8.5.0-yzt7mdy
# umpire@=6.0.0%gcc@=8.5.0+c+cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp~rocm~shared build_system=cmake build_type=Release cuda_arch=70 generator=make tests=none arch=linux-centos8-power9le
module load umpire/6.0.0-gcc-8.5.0-3hm2whr
# hiop@=develop%gcc@=8.5.0+cuda~cusolver_lu~deepchecking+ginkgo~ipo~jsrun+kron+mpi+raja~rocm~shared+sparse build_system=cmake build_type=Debug cuda_arch=70 dev_path=/people/svcexasgd/gitlab/26190/spack_newell/hiop_dev generator=make arch=linux-centos8-power9le
module load hiop/develop-gcc-8.5.0-7ck7er5
# ipopt@=3.12.10%gcc@=8.5.0+coinhsl~debug~metis~mumps build_system=autotools arch=linux-centos8-power9le
module load ipopt/3.12.10-gcc-8.5.0-rucines
# python@=3.8.5%gcc@=8.5.0+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=0d98e93,4c24573,ebdca64,f2fd060 arch=linux-centos8-power9le
module load python/3.8.5-gcc-8.5.0-dwf3kgs
# petsc@=3.20.1%gcc@=8.5.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~sycl~tetgen~trilinos~valgrind build_system=generic clanguage=C memalign=none arch=linux-centos8-power9le
module load petsc/3.20.1-gcc-8.5.0-6wcyj55
# py-pip@=23.1.2%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load py-pip/23.1.2-gcc-8.5.0-45skzep
# py-setuptools@=68.0.0%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load py-setuptools/68.0.0-gcc-8.5.0-be735hh
# py-wheel@=0.41.2%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load py-wheel/0.41.2-gcc-8.5.0-o5hqddm
# py-cython@=0.29.36%gcc@=8.5.0 build_system=python_pip patches=c4369ad arch=linux-centos8-power9le
module load py-cython/0.29.36-gcc-8.5.0-4lthysc
# py-mpi4py@=3.1.4%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-mpi4py/3.1.4-gcc-8.5.0-aaxbl2w
# py-flit-core@=3.9.0%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-flit-core/3.9.0-gcc-8.5.0-7qlf3hc
# libmd@=1.0.4%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load libmd/1.0.4-gcc-8.5.0-s5jz3si
# libbsd@=0.11.7%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load libbsd/0.11.7-gcc-8.5.0-nqmfv4p
# expat@=2.5.0%gcc@=8.5.0+libbsd build_system=autotools arch=linux-centos8-power9le
module load expat/2.5.0-gcc-8.5.0-dl7oihc
# libunistring@=1.1%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load libunistring/1.1-gcc-8.5.0-y3gjdbv
# libidn2@=2.3.4%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load libidn2/2.3.4-gcc-8.5.0-seplt7v
# bison@=3.8.2%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load bison/3.8.2-gcc-8.5.0-p53xq7y
# findutils@=4.9.0%gcc@=8.5.0 build_system=autotools patches=440b954 arch=linux-centos8-power9le
module load findutils/4.9.0-gcc-8.5.0-j6yj47y
# krb5@=1.20.1%gcc@=8.5.0+shared build_system=autotools arch=linux-centos8-power9le
module load krb5/1.20.1-gcc-8.5.0-lqo7ggz
# libedit@=3.1-20210216%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load libedit/3.1-20210216-gcc-8.5.0-xz3heyj
# libxcrypt@=4.4.35%gcc@=8.5.0~obsolete_api build_system=autotools patches=4885da3 arch=linux-centos8-power9le
module load libxcrypt/4.4.35-gcc-8.5.0-nnuqrig
# openssh@=9.5p1%gcc@=8.5.0+gssapi build_system=autotools arch=linux-centos8-power9le
module load openssh/9.5p1-gcc-8.5.0-rq4fwuu
# pcre2@=10.42%gcc@=8.5.0~jit+multibyte build_system=autotools arch=linux-centos8-power9le
module load pcre2/10.42-gcc-8.5.0-lvf7rst
# git@=2.42.0%gcc@=8.5.0+man+nls+perl+subtree~svn~tcltk build_system=autotools arch=linux-centos8-power9le
module load git/2.42.0-gcc-8.5.0-b5ddzfl
# py-packaging@=23.1%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-packaging/23.1-gcc-8.5.0-fupho4o
# py-tomli@=2.0.1%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-tomli/2.0.1-gcc-8.5.0-a4g62fz
# py-typing-extensions@=4.8.0%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-typing-extensions/4.8.0-gcc-8.5.0-usi6lq5
# py-setuptools-scm@=7.1.0%gcc@=8.5.0+toml build_system=python_pip arch=linux-centos8-power9le
module load py-setuptools-scm/7.1.0-gcc-8.5.0-u7dcrps
# py-flit-scm@=1.7.0%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-flit-scm/1.7.0-gcc-8.5.0-ov6nqzs
# py-exceptiongroup@=1.1.1%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-exceptiongroup/1.1.1-gcc-8.5.0-uk34au2
# py-editables@=0.3%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-editables/0.3-gcc-8.5.0-dcmsfjz
# py-pathspec@=0.11.1%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-pathspec/0.11.1-gcc-8.5.0-7xbvbhg
# py-pluggy@=1.0.0%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-pluggy/1.0.0-gcc-8.5.0-wlbjopm
# py-calver@=2022.6.26%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-calver/2022.6.26-gcc-8.5.0-hgxbsu2
# py-trove-classifiers@=2023.8.7%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-trove-classifiers/2023.8.7-gcc-8.5.0-5wwsb65
# py-hatchling@=1.18.0%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-hatchling/1.18.0-gcc-8.5.0-agn6can
# py-hatch-vcs@=0.3.0%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-hatch-vcs/0.3.0-gcc-8.5.0-aiflljq
# py-iniconfig@=2.0.0%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-iniconfig/2.0.0-gcc-8.5.0-f45q4wi
# py-pytest@=7.3.2%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-pytest/7.3.2-gcc-8.5.0-fnzfnsk
# exago@=develop%gcc@=8.5.0+cuda+hiop~ipo+ipopt+logging+mpi+python+raja~rocm build_system=cmake build_type=Debug cuda_arch=70 dev_path=/people/svcexasgd/gitlab/26190/spack_newell generator=make arch=linux-centos8-power9le
## module load exago/develop-gcc-8.5.0-ypyriqf
