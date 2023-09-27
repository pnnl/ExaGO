module use -a /gpfs/alpine/proj-shared/csc359/cameron/spack-install/summit-modules/linux-rhel8-power9le
# cmake@=3.21.3%gcc@=10.2.0~doc+ncurses+ownlibs~qt build_type=Release arch=linux-rhel8-power9le
module load cmake/3.21.3-gcc-10.2.0-ijmjum5
# blt@=0.4.1%gcc@=10.2.0 arch=linux-rhel8-power9le
module load blt/0.4.1-gcc-10.2.0-ymjl6ab
# cub@=1.16.0%gcc@=10.2.0 arch=linux-rhel8-power9le
module load cub/1.16.0-gcc-10.2.0-a5uzxeq
# cuda@=11.4.2%gcc@=10.2.0~allow-unsupported-compilers~dev arch=linux-rhel8-power9le
module load cuda/11.4.2-gcc-10.2.0-4bvhk3u
# camp@=0.2.3%gcc@=10.2.0+cuda~ipo+openmp~rocm~tests build_type=RelWithDebInfo cuda_arch=70 arch=linux-rhel8-power9le
module load camp/0.2.3-gcc-10.2.0-fnpymoo
# gnuconfig@=2021-08-14%gcc@=10.2.0 arch=linux-rhel8-power9le
module load gnuconfig/2021-08-14-gcc-10.2.0-4gsxwpk
# gmake@=4.4.1%gcc@=10.2.0~guile build_system=autotools arch=linux-rhel8-power9le
module load gmake/4.4.1-gcc-10.2.0-m2nvote
# openblas@=0.17.0%gcc@=10.2.0~bignuma~consistent_fpcsr+fortran~ilp64+locking+pic+shared symbol_suffix=none threads=none arch=linux-rhel8-power9le
module load openblas/0.17.0-gcc-10.2.0-qcs6vhe
# magma@=2.6.2%gcc@=10.2.0+cuda+fortran~ipo~rocm+shared build_type=Release cuda_arch=70 arch=linux-rhel8-power9le
module load magma/2.6.2-gcc-10.2.0-wfpnmgj
# raja@=0.14.0%gcc@=10.2.0+cuda+examples+exercises~ipo+openmp~rocm+shared~tests build_type=Release cuda_arch=70 arch=linux-rhel8-power9le
module load raja/0.14.0-gcc-10.2.0-r7kozj7
# spectrum-mpi@=10.4.0.3-20210112%gcc@=10.2.0 arch=linux-rhel8-power9le
module load spectrum-mpi/10.4.0.3-20210112-gcc-10.2.0-rhfl5sr
# umpire@=6.0.0%gcc@=10.2.0~c+cuda~device_alloc~deviceconst+examples~fortran~ipo~numa~openmp~rocm~shared build_type=Release cuda_arch=70 tests=none arch=linux-rhel8-power9le
module load umpire/6.0.0-gcc-10.2.0-eaxxvet
# hiop@=1.0.0%gcc@=10.2.0+cuda~cusolver_lu~deepchecking~ginkgo~ipo~jsrun~kron+mpi+raja~rocm~shared~sparse build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-rhel8-power9le
module load hiop/1.0.0-gcc-10.2.0-pomojvf
# coinhsl@=2015.06.23%gcc@=10.2.0+blas arch=linux-rhel8-power9le
module load coinhsl/2015.06.23-gcc-10.2.0-vctp77v
# metis@=5.1.0%gcc@=10.2.0~gdb~int64~real64+shared build_type=Release patches=4991da9,b1225da arch=linux-rhel8-power9le
module load metis/5.1.0-gcc-10.2.0-lpkktw5
# pkgconf@=1.8.0%gcc@=10.2.0 arch=linux-rhel8-power9le
module load pkgconf/1.8.0-gcc-10.2.0-ygtljd4
# ipopt@=3.12.10%gcc@=10.2.0+coinhsl~debug+metis~mumps arch=linux-rhel8-power9le
module load ipopt/3.12.10-gcc-10.2.0-zl5xybc
# libiconv@=1.16%gcc@=10.2.0 libs=shared,static arch=linux-rhel8-power9le
module load libiconv/1.16-gcc-10.2.0-wyjuzaw
# diffutils@=3.8%gcc@=10.2.0 arch=linux-rhel8-power9le
module load diffutils/3.8-gcc-10.2.0-gkkihj7
# parmetis@=4.0.3%gcc@=10.2.0~gdb~int64~ipo+shared build_type=RelWithDebInfo patches=4f89253,50ed208,704b84f arch=linux-rhel8-power9le
module load parmetis/4.0.3-gcc-10.2.0-by3w24c
# python@=3.8.10%gcc@=10.2.0+bz2+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tix~tkinter~ucs4+uuid+zlib patches=0d98e93,4c24573,f2fd060 arch=linux-rhel8-power9le
module load python/3.8.10-gcc-10.2.0-yeap5du
# petsc@=3.18.0%gcc@=10.2.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind+metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~tetgen~trilinos~valgrind clanguage=C patches=2daeca7 arch=linux-rhel8-power9le
module load petsc/3.18.0-gcc-10.2.0-6zx23r5
# python@=3.8.10%gcc@=10.2.0+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=0d98e93,4c24573,f2fd060 arch=linux-rhel8-power9le
module load python/3.8.10-gcc-10.2.0-5pl4xoh
# py-pip@=23.1.2%gcc@=10.2.0 build_system=generic arch=linux-rhel8-power9le
module load py-pip/23.1.2-gcc-10.2.0-tk6lrbb
# py-setuptools@=68.0.0%gcc@=10.2.0 build_system=generic arch=linux-rhel8-power9le
module load py-setuptools/68.0.0-gcc-10.2.0-52xm36g
# py-wheel@=0.37.1%gcc@=10.2.0 build_system=generic arch=linux-rhel8-power9le
module load py-wheel/0.37.1-gcc-10.2.0-gzdghyv
# py-cython@=0.29.36%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-cython/0.29.36-gcc-10.2.0-futh2n7
# py-mpi4py@=3.1.4%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-mpi4py/3.1.4-gcc-10.2.0-giulifn
# py-flit-core@=3.9.0%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-flit-core/3.9.0-gcc-10.2.0-rkvgviw
# libsigsegv@=2.13%gcc@=10.2.0 arch=linux-rhel8-power9le
module load libsigsegv/2.13-gcc-10.2.0-7j4qaaw
# m4@=1.4.19%gcc@=10.2.0+sigsegv patches=9dc5fbd,bfdffa7 arch=linux-rhel8-power9le
module load m4/1.4.19-gcc-10.2.0-meue3ml
# perl@=5.30.1%gcc@=10.2.0+cpanm+opcode+open+shared+threads build_system=generic arch=linux-rhel8-power9le
module load perl/5.30.1-gcc-10.2.0-jh4r2fu
# autoconf@=2.69%gcc@=10.2.0 build_system=autotools patches=35c4492,7793209,a49dd5b arch=linux-rhel8-power9le
module load autoconf/2.69-gcc-10.2.0-45t46f3
# automake@=1.16.5%gcc@=10.2.0 build_system=autotools arch=linux-rhel8-power9le
module load automake/1.16.5-gcc-10.2.0-qsilt77
# nghttp2@=1.52.0%gcc@=10.2.0 build_system=autotools arch=linux-rhel8-power9le
module load nghttp2/1.52.0-gcc-10.2.0-ryjqxcn
# ca-certificates-mozilla@=2023-05-30%gcc@=10.2.0 build_system=generic arch=linux-rhel8-power9le
module load ca-certificates-mozilla/2023-05-30-gcc-10.2.0-v7vz3c5
# zlib-ng@=2.1.3%gcc@=10.2.0+compat+opt build_system=autotools patches=299b958,ae9077a,b692621 arch=linux-rhel8-power9le
module load zlib-ng/2.1.3-gcc-10.2.0-56obt5c
# openssl@=3.1.2%gcc@=10.2.0~docs+shared build_system=generic certs=mozilla arch=linux-rhel8-power9le
## module load openssl/3.1.2-gcc-10.2.0-tlszcll
# curl@=8.1.2%gcc@=10.2.0~gssapi~ldap~libidn2~librtmp~libssh~libssh2+nghttp2 build_system=autotools libs=shared,static tls=openssl arch=linux-rhel8-power9le
module load curl/8.1.2-gcc-10.2.0-qtqogpi
# libmd@=1.0.4%gcc@=10.2.0 build_system=autotools arch=linux-rhel8-power9le
module load libmd/1.0.4-gcc-10.2.0-vc6uybq
# libbsd@=0.11.7%gcc@=10.2.0 build_system=autotools arch=linux-rhel8-power9le
module load libbsd/0.11.7-gcc-10.2.0-sfbxsnu
# expat@=2.5.0%gcc@=10.2.0+libbsd build_system=autotools arch=linux-rhel8-power9le
module load expat/2.5.0-gcc-10.2.0-qz7ylbw
# bzip2@=1.0.8%gcc@=10.2.0~debug~pic+shared build_system=generic arch=linux-rhel8-power9le
module load bzip2/1.0.8-gcc-10.2.0-2hkazde
# xz@=5.4.1%gcc@=10.2.0~pic build_system=autotools libs=shared,static arch=linux-rhel8-power9le
module load xz/5.4.1-gcc-10.2.0-fi6tzmu
# libxml2@=2.10.3%gcc@=10.2.0~python build_system=autotools arch=linux-rhel8-power9le
module load libxml2/2.10.3-gcc-10.2.0-unxhqfx
# ncurses@=6.4%gcc@=10.2.0~symlinks+termlib abi=none build_system=autotools arch=linux-rhel8-power9le
module load ncurses/6.4-gcc-10.2.0-kw5rm5w
# pigz@=2.7%gcc@=10.2.0 build_system=makefile arch=linux-rhel8-power9le
module load pigz/2.7-gcc-10.2.0-3wopxwk
# zstd@=1.5.5%gcc@=10.2.0+programs build_system=makefile compression=none libs=shared,static arch=linux-rhel8-power9le
module load zstd/1.5.5-gcc-10.2.0-uxatbd2
# tar@=1.34%gcc@=10.2.0 build_system=autotools zip=pigz arch=linux-rhel8-power9le
module load tar/1.34-gcc-10.2.0-sfabuja
# gettext@=0.21.1%gcc@=10.2.0+bzip2+curses+git~libunistring+libxml2+tar+xz build_system=autotools arch=linux-rhel8-power9le
module load gettext/0.21.1-gcc-10.2.0-t5i6mgx
# libunistring@=1.1%gcc@=10.2.0 build_system=autotools arch=linux-rhel8-power9le
module load libunistring/1.1-gcc-10.2.0-dbahaeb
# libidn2@=2.3.4%gcc@=10.2.0 build_system=autotools arch=linux-rhel8-power9le
module load libidn2/2.3.4-gcc-10.2.0-vydwqsc
# libtool@=2.4.7%gcc@=10.2.0 arch=linux-rhel8-power9le
module load libtool/2.4.7-gcc-10.2.0-n2qpycy
# bison@=3.8.2%gcc@=10.2.0 build_system=autotools arch=linux-rhel8-power9le
module load bison/3.8.2-gcc-10.2.0-ijrzj6d
# findutils@=4.9.0%gcc@=10.2.0 build_system=autotools patches=440b954 arch=linux-rhel8-power9le
module load findutils/4.9.0-gcc-10.2.0-lr4vpa7
# krb5@=1.20.1%gcc@=10.2.0+shared build_system=autotools arch=linux-rhel8-power9le
module load krb5/1.20.1-gcc-10.2.0-2fd22ss
# libedit@=3.1-20210216%gcc@=10.2.0 build_system=autotools arch=linux-rhel8-power9le
module load libedit/3.1-20210216-gcc-10.2.0-zvx7vgn
# libxcrypt@=4.4.35%gcc@=10.2.0~obsolete_api build_system=autotools patches=4885da3 arch=linux-rhel8-power9le
module load libxcrypt/4.4.35-gcc-10.2.0-cgn72av
# openssh@=9.3p1%gcc@=10.2.0+gssapi build_system=autotools arch=linux-rhel8-power9le
module load openssh/9.3p1-gcc-10.2.0-ah5wkbh
# pcre2@=10.42%gcc@=10.2.0~jit+multibyte build_system=autotools arch=linux-rhel8-power9le
module load pcre2/10.42-gcc-10.2.0-sss2r4u
# git@=2.41.0%gcc@=10.2.0+man+nls+perl+subtree~svn~tcltk build_system=autotools arch=linux-rhel8-power9le
module load git/2.41.0-gcc-10.2.0-h4x4wvc
# py-packaging@=23.1%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-packaging/23.1-gcc-10.2.0-duwxwhx
# py-tomli@=2.0.1%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-tomli/2.0.1-gcc-10.2.0-5rf2kwr
# py-typing-extensions@=4.6.3%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-typing-extensions/4.6.3-gcc-10.2.0-fn237ld
# py-setuptools-scm@=7.1.0%gcc@=10.2.0+toml build_system=python_pip arch=linux-rhel8-power9le
module load py-setuptools-scm/7.1.0-gcc-10.2.0-jiqp4sa
# py-flit-scm@=1.7.0%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-flit-scm/1.7.0-gcc-10.2.0-m4j6dwa
# py-exceptiongroup@=1.1.1%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-exceptiongroup/1.1.1-gcc-10.2.0-5xyi4i2
# py-editables@=0.3%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-editables/0.3-gcc-10.2.0-poertub
# py-pathspec@=0.11.1%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-pathspec/0.11.1-gcc-10.2.0-zhfln35
# py-pluggy@=1.0.0%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-pluggy/1.0.0-gcc-10.2.0-c5x5qj2
# py-calver@=2022.6.26%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-calver/2022.6.26-gcc-10.2.0-kmyz4d5
# py-trove-classifiers@=2023.3.9%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-trove-classifiers/2023.3.9-gcc-10.2.0-tpphxxs
# py-hatchling@=1.17.0%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-hatchling/1.17.0-gcc-10.2.0-e5iupfm
# py-hatch-vcs@=0.3.0%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-hatch-vcs/0.3.0-gcc-10.2.0-n4c3tt6
# py-iniconfig@=2.0.0%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-iniconfig/2.0.0-gcc-10.2.0-n5eln73
# py-pytest@=7.3.2%gcc@=10.2.0 build_system=python_pip arch=linux-rhel8-power9le
module load py-pytest/7.3.2-gcc-10.2.0-l2ofgsi
# exago@=develop%gcc@=10.2.0+cuda+hiop~ipo+ipopt+logging+mpi+python+raja~rocm build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-rhel8-power9le
## module load exago/develop-gcc-10.2.0-lpqzdxo
