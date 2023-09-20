module use -a /qfs/projects/exasgd/src/ci-newll/spack-install/ci-modules/linux-centos8-power9le
# gnuconfig@=2022-09-17%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load gnuconfig/2022-09-17-gcc-8.5.0-brd4tz7
# pkgconf@=1.9.5%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load pkgconf/1.9.5-gcc-8.5.0-7iwzavg
# nghttp2@=1.52.0%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load nghttp2/1.52.0-gcc-8.5.0-djfawnd
# ca-certificates-mozilla@=2023-05-30%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load ca-certificates-mozilla/2023-05-30-gcc-8.5.0-ueq22oa
# berkeley-db@=18.1.40%gcc@=8.5.0+cxx~docs+stl build_system=autotools patches=26090f4,b231fcc arch=linux-centos8-power9le
module load berkeley-db/18.1.40-gcc-8.5.0-z3z2uin
# libiconv@=1.17%gcc@=8.5.0 build_system=autotools libs=shared,static arch=linux-centos8-power9le
module load libiconv/1.17-gcc-8.5.0-nzdc6r3
# diffutils@=3.9%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load diffutils/3.9-gcc-8.5.0-nlw7zfo
# bzip2@=1.0.8%gcc@=8.5.0~debug~pic+shared build_system=generic arch=linux-centos8-power9le
module load bzip2/1.0.8-gcc-8.5.0-fbtjslf
# ncurses@=6.4%gcc@=8.5.0~symlinks+termlib abi=none build_system=autotools arch=linux-centos8-power9le
module load ncurses/6.4-gcc-8.5.0-bfcysa2
# readline@=8.2%gcc@=8.5.0 build_system=autotools patches=bbf97f1 arch=linux-centos8-power9le
module load readline/8.2-gcc-8.5.0-yku2ukd
# gdbm@=1.23%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load gdbm/1.23-gcc-8.5.0-tfwl4o4
# zlib-ng@=2.1.3%gcc@=8.5.0+compat+opt build_system=autotools patches=299b958,ae9077a,b692621 arch=linux-centos8-power9le
module load zlib-ng/2.1.3-gcc-8.5.0-kuviknf
# perl@=5.38.0%gcc@=8.5.0+cpanm+opcode+open+shared+threads build_system=generic patches=714e4d1 arch=linux-centos8-power9le
module load perl/5.38.0-gcc-8.5.0-hlhx7fa
# openssl@=3.1.2%gcc@=8.5.0~docs+shared build_system=generic certs=mozilla arch=linux-centos8-power9le
## module load openssl/3.1.2-gcc-8.5.0-zfrzvxo
# curl@=8.1.2%gcc@=8.5.0~gssapi~ldap~libidn2~librtmp~libssh~libssh2+nghttp2 build_system=autotools libs=shared,static tls=openssl arch=linux-centos8-power9le
module load curl/8.1.2-gcc-8.5.0-kd3kffs
# cmake@=3.27.4%gcc@=8.5.0~doc+ncurses+ownlibs build_system=generic build_type=Release arch=linux-centos8-power9le
module load cmake/3.27.4-gcc-8.5.0-zccvhqs
# blt@=0.4.1%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load blt/0.4.1-gcc-8.5.0-odrdc4r
# cub@=2.1.0%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load cub/2.1.0-gcc-8.5.0-lx5hjs2
# cuda@=11.4%gcc@=8.5.0~allow-unsupported-compilers~dev build_system=generic arch=linux-centos8-power9le
module load cuda/11.4-gcc-8.5.0-vfytmhm
# gmake@=4.4.1%gcc@=8.5.0~guile build_system=autotools arch=linux-centos8-power9le
module load gmake/4.4.1-gcc-8.5.0-vdsjgdd
# camp@=0.2.3%gcc@=8.5.0+cuda~ipo+openmp~rocm~tests build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-centos8-power9le
module load camp/0.2.3-gcc-8.5.0-27w2i42
# openblas@=0.3.23%gcc@=8.5.0~bignuma~consistent_fpcsr+fortran~ilp64+locking+pic+shared build_system=makefile symbol_suffix=none threads=none arch=linux-centos8-power9le
module load openblas/0.3.23-gcc-8.5.0-jyppjuy
# coinhsl@=2019.05.21%gcc@=8.5.0+blas build_system=autotools arch=linux-centos8-power9le
module load coinhsl/2019.05.21-gcc-8.5.0-5wx4j3p
# ginkgo@=1.5.0.glu_experimental%gcc@=8.5.0+cuda~develtools~full_optimizations~hwloc~ipo~mpi~oneapi+openmp~rocm+shared build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-centos8-power9le
module load ginkgo/1.5.0.glu_experimental-gcc-8.5.0-23k5d3p
# magma@=2.6.2%gcc@=8.5.0+cuda+fortran~ipo~rocm+shared build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-centos8-power9le
module load magma/2.6.2-gcc-8.5.0-vlyt3nf
# metis@=5.1.0%gcc@=8.5.0~gdb~int64~ipo~real64+shared build_system=cmake build_type=Release generator=make patches=4991da9,93a7903,b1225da arch=linux-centos8-power9le
module load metis/5.1.0-gcc-8.5.0-kilafih
# openmpi@=4.1.4%gcc@=8.5.0~atomics~cuda~cxx~cxx_exceptions~gpfs~internal-hwloc~internal-pmix~java~legacylaunchers~lustre~memchecker~openshmem~orterunprefix+romio+rsh~singularity+static+vt+wrapper-rpath build_system=autotools fabrics=none schedulers=none arch=linux-centos8-power9le
module load openmpi/4.1.4-gcc-8.5.0-b4va3xy
# raja@=0.14.0%gcc@=8.5.0+cuda~examples~exercises~ipo+openmp~rocm+shared~tests build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-centos8-power9le
module load raja/0.14.0-gcc-8.5.0-qz5ppvh
# libsigsegv@=2.14%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load libsigsegv/2.14-gcc-8.5.0-fw5y3ai
# m4@=1.4.19%gcc@=8.5.0+sigsegv build_system=autotools patches=9dc5fbd,bfdffa7 arch=linux-centos8-power9le
module load m4/1.4.19-gcc-8.5.0-4r5k5rv
# autoconf@=2.69%gcc@=8.5.0 build_system=autotools patches=35c4492,7793209,a49dd5b arch=linux-centos8-power9le
module load autoconf/2.69-gcc-8.5.0-wpknfu7
# automake@=1.16.5%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load automake/1.16.5-gcc-8.5.0-3ml7b4b
# libtool@=2.4.7%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load libtool/2.4.7-gcc-8.5.0-ouhpeo3
# gmp@=6.2.1%gcc@=8.5.0+cxx build_system=autotools libs=shared,static patches=69ad2e2 arch=linux-centos8-power9le
module load gmp/6.2.1-gcc-8.5.0-bswbnus
# autoconf-archive@=2023.02.20%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load autoconf-archive/2023.02.20-gcc-8.5.0-tk6s3qc
# xz@=5.4.1%gcc@=8.5.0~pic build_system=autotools libs=shared,static arch=linux-centos8-power9le
module load xz/5.4.1-gcc-8.5.0-4wstnwi
# libxml2@=2.10.3%gcc@=8.5.0~python build_system=autotools arch=linux-centos8-power9le
module load libxml2/2.10.3-gcc-8.5.0-6efn6lp
# pigz@=2.7%gcc@=8.5.0 build_system=makefile arch=linux-centos8-power9le
module load pigz/2.7-gcc-8.5.0-wya247t
# zstd@=1.5.5%gcc@=8.5.0+programs build_system=makefile compression=none libs=shared,static arch=linux-centos8-power9le
module load zstd/1.5.5-gcc-8.5.0-apvtmrg
# tar@=1.34%gcc@=8.5.0 build_system=autotools zip=pigz arch=linux-centos8-power9le
module load tar/1.34-gcc-8.5.0-rzjrqwk
# gettext@=0.21.1%gcc@=8.5.0+bzip2+curses+git~libunistring+libxml2+tar+xz build_system=autotools arch=linux-centos8-power9le
module load gettext/0.21.1-gcc-8.5.0-42wk4og
# texinfo@=7.0.3%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load texinfo/7.0.3-gcc-8.5.0-xp6rohx
# mpfr@=4.2.0%gcc@=8.5.0 build_system=autotools libs=shared,static arch=linux-centos8-power9le
module load mpfr/4.2.0-gcc-8.5.0-tyh35ce
# suite-sparse@=5.13.0%gcc@=8.5.0~cuda~graphblas~openmp+pic build_system=generic arch=linux-centos8-power9le
module load suite-sparse/5.13.0-gcc-8.5.0-uft6tn5
# umpire@=6.0.0%gcc@=8.5.0+c+cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp~rocm~shared build_system=cmake build_type=Release cuda_arch=70 generator=make tests=none arch=linux-centos8-power9le
module load umpire/6.0.0-gcc-8.5.0-lh4bj5g
# hiop@=develop%gcc@=8.5.0+cuda~cusolver_lu~deepchecking+ginkgo~ipo~jsrun+kron+mpi+raja~rocm~shared+sparse build_system=cmake build_type=Debug cuda_arch=70 generator=make arch=linux-centos8-power9le
module load hiop/develop-gcc-8.5.0-26hrqbh
# ipopt@=3.12.10%gcc@=8.5.0+coinhsl~debug~metis~mumps build_system=autotools arch=linux-centos8-power9le
module load ipopt/3.12.10-gcc-8.5.0-4hqffpd
# python@=3.8.5%gcc@=8.5.0+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=0d98e93,4c24573,f2fd060 arch=linux-centos8-power9le
module load python/3.8.5-gcc-8.5.0-ogeppsp
# petsc@=3.19.4%gcc@=8.5.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~tetgen~trilinos~valgrind build_system=generic clanguage=C memalign=none arch=linux-centos8-power9le
module load petsc/3.19.4-gcc-8.5.0-bevhrwm
# py-pip@=23.1.2%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load py-pip/23.1.2-gcc-8.5.0-4kzc7dp
# py-setuptools@=68.0.0%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load py-setuptools/68.0.0-gcc-8.5.0-eo2egu5
# py-wheel@=0.37.1%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load py-wheel/0.37.1-gcc-8.5.0-hs4dzbb
# py-cython@=0.29.36%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-cython/0.29.36-gcc-8.5.0-ugkqxqj
# py-mpi4py@=3.1.4%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-mpi4py/3.1.4-gcc-8.5.0-hymjh25
# py-flit-core@=3.9.0%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-flit-core/3.9.0-gcc-8.5.0-5vymygb
# libmd@=1.0.4%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load libmd/1.0.4-gcc-8.5.0-ychwxkl
# libbsd@=0.11.7%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load libbsd/0.11.7-gcc-8.5.0-67zh7gq
# expat@=2.5.0%gcc@=8.5.0+libbsd build_system=autotools arch=linux-centos8-power9le
module load expat/2.5.0-gcc-8.5.0-ijq6ufa
# libunistring@=1.1%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load libunistring/1.1-gcc-8.5.0-bjw2hvg
# libidn2@=2.3.4%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load libidn2/2.3.4-gcc-8.5.0-nwmhm35
# bison@=3.8.2%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load bison/3.8.2-gcc-8.5.0-2ift2tq
# findutils@=4.9.0%gcc@=8.5.0 build_system=autotools patches=440b954 arch=linux-centos8-power9le
module load findutils/4.9.0-gcc-8.5.0-plaf76d
# krb5@=1.20.1%gcc@=8.5.0+shared build_system=autotools arch=linux-centos8-power9le
module load krb5/1.20.1-gcc-8.5.0-gkiv2fi
# libedit@=3.1-20210216%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load libedit/3.1-20210216-gcc-8.5.0-6pxpqfh
# libxcrypt@=4.4.35%gcc@=8.5.0~obsolete_api build_system=autotools patches=4885da3 arch=linux-centos8-power9le
module load libxcrypt/4.4.35-gcc-8.5.0-nqk6zhv
# openssh@=9.3p1%gcc@=8.5.0+gssapi build_system=autotools arch=linux-centos8-power9le
module load openssh/9.3p1-gcc-8.5.0-2r62a3r
# pcre2@=10.42%gcc@=8.5.0~jit+multibyte build_system=autotools arch=linux-centos8-power9le
module load pcre2/10.42-gcc-8.5.0-egehak2
# git@=2.41.0%gcc@=8.5.0+man+nls+perl+subtree~svn~tcltk build_system=autotools arch=linux-centos8-power9le
module load git/2.41.0-gcc-8.5.0-ft37kpd
# py-packaging@=23.1%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-packaging/23.1-gcc-8.5.0-xcviq66
# py-tomli@=2.0.1%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-tomli/2.0.1-gcc-8.5.0-yw34kvp
# py-typing-extensions@=4.6.3%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-typing-extensions/4.6.3-gcc-8.5.0-6zasnpy
# py-setuptools-scm@=7.1.0%gcc@=8.5.0+toml build_system=python_pip arch=linux-centos8-power9le
module load py-setuptools-scm/7.1.0-gcc-8.5.0-x3pe74m
# py-flit-scm@=1.7.0%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-flit-scm/1.7.0-gcc-8.5.0-4bxhula
# py-exceptiongroup@=1.1.1%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-exceptiongroup/1.1.1-gcc-8.5.0-7kvna4h
# py-editables@=0.3%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-editables/0.3-gcc-8.5.0-7l4cwk6
# py-pathspec@=0.11.1%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-pathspec/0.11.1-gcc-8.5.0-u2routk
# py-pluggy@=1.0.0%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-pluggy/1.0.0-gcc-8.5.0-ssp5oep
# py-calver@=2022.6.26%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-calver/2022.6.26-gcc-8.5.0-p4pnje4
# py-trove-classifiers@=2023.3.9%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-trove-classifiers/2023.3.9-gcc-8.5.0-m3hs64b
# py-hatchling@=1.17.0%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-hatchling/1.17.0-gcc-8.5.0-5brbik5
# py-hatch-vcs@=0.3.0%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-hatch-vcs/0.3.0-gcc-8.5.0-lyq3fgv
# py-iniconfig@=2.0.0%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-iniconfig/2.0.0-gcc-8.5.0-6exw3my
# py-pytest@=7.3.2%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-pytest/7.3.2-gcc-8.5.0-hafdero
# exago@=develop%gcc@=8.5.0+cuda+hiop~ipo+ipopt+logging+mpi+python+raja~rocm build_system=cmake build_type=MinSizeRel cuda_arch=70 dev_path=/people/svcexasgd/gitlab/22741/spack_newell generator=make arch=linux-centos8-power9le
## module load exago/develop-gcc-8.5.0-dfiemk3
