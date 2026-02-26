module use -a /qfs/projects/earthshot/src/newell-ci/spack-install/ci-modules/linux-centos8-power9le
# gcc-runtime@=8.5.0%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load gcc-runtime/8.5.0-gcc-8.5.0-pmgzgcn
# gmake@=4.4.1%gcc@=8.5.0~guile build_system=generic arch=linux-centos8-power9le
module load gmake/4.4.1-gcc-8.5.0-3binpl5
# gnuconfig@=2022-09-17%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load gnuconfig/2022-09-17-gcc-8.5.0-amqajvg
# pkgconf@=1.9.5%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load pkgconf/1.9.5-gcc-8.5.0-f3il7dq
# nghttp2@=1.57.0%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load nghttp2/1.57.0-gcc-8.5.0-awrz4bk
# ca-certificates-mozilla@=2023-05-30%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load ca-certificates-mozilla/2023-05-30-gcc-8.5.0-lw3r7th
# berkeley-db@=18.1.40%gcc@=8.5.0+cxx~docs+stl build_system=autotools patches=26090f4,b231fcc arch=linux-centos8-power9le
module load berkeley-db/18.1.40-gcc-8.5.0-u5uqtkn
# libiconv@=1.17%gcc@=8.5.0 build_system=autotools libs=shared,static arch=linux-centos8-power9le
module load libiconv/1.17-gcc-8.5.0-tteytxo
# diffutils@=3.9%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load diffutils/3.9-gcc-8.5.0-3ggffsl
# bzip2@=1.0.8%gcc@=8.5.0~debug~pic+shared build_system=generic arch=linux-centos8-power9le
module load bzip2/1.0.8-gcc-8.5.0-voqeke2
# ncurses@=6.4%gcc@=8.5.0~symlinks+termlib abi=none build_system=autotools arch=linux-centos8-power9le
module load ncurses/6.4-gcc-8.5.0-kczmwq4
# readline@=8.2%gcc@=8.5.0 build_system=autotools patches=bbf97f1 arch=linux-centos8-power9le
module load readline/8.2-gcc-8.5.0-ynosqes
# gdbm@=1.23%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load gdbm/1.23-gcc-8.5.0-ypnlczn
# zlib-ng@=2.1.5%gcc@=8.5.0+compat+opt build_system=autotools arch=linux-centos8-power9le
module load zlib-ng/2.1.5-gcc-8.5.0-dsplf72
# perl@=5.38.0%gcc@=8.5.0+cpanm+opcode+open+shared+threads build_system=generic patches=714e4d1 arch=linux-centos8-power9le
module load perl/5.38.0-gcc-8.5.0-h2oo6ir
# openssl@=3.1.3%gcc@=8.5.0~docs+shared build_system=generic certs=mozilla arch=linux-centos8-power9le
## module load openssl/3.1.3-gcc-8.5.0-pswhwzo
# curl@=8.4.0%gcc@=8.5.0~gssapi~ldap~libidn2~librtmp~libssh~libssh2+nghttp2 build_system=autotools libs=shared,static tls=openssl arch=linux-centos8-power9le
module load curl/8.4.0-gcc-8.5.0-mwzlu55
# cmake@=3.27.9%gcc@=8.5.0~doc+ncurses+ownlibs build_system=generic build_type=Release arch=linux-centos8-power9le
module load cmake/3.27.9-gcc-8.5.0-3vm6qun
# blt@=0.4.1%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load blt/0.4.1-gcc-8.5.0-qbqjtr2
# cub@=2.1.0%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load cub/2.1.0-gcc-8.5.0-n4h3qsz
# cuda@=11.4%gcc@=8.5.0~allow-unsupported-compilers~dev build_system=generic arch=linux-centos8-power9le
module load cuda/11.4-gcc-8.5.0-fcxmvvd
# camp@=0.2.3%gcc@=8.5.0+cuda~ipo+openmp~rocm~tests build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-centos8-power9le
module load camp/0.2.3-gcc-8.5.0-me7qrtf
# ginkgo@=1.5.0.glu_experimental%gcc@=8.5.0+cuda~develtools~full_optimizations~hwloc~ipo~mpi+openmp~rocm+shared~sycl build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-centos8-power9le
module load ginkgo/1.5.0.glu_experimental-gcc-8.5.0-g5n6tgr
# openblas@=0.3.26%gcc@=8.5.0~bignuma~consistent_fpcsr+fortran~ilp64+locking+pic+shared build_system=makefile symbol_suffix=none threads=none arch=linux-centos8-power9le
module load openblas/0.3.26-gcc-8.5.0-atmmzrk
# coinhsl@=2019.05.21%gcc@=8.5.0+blas build_system=autotools arch=linux-centos8-power9le
module load coinhsl/2019.05.21-gcc-8.5.0-p7hwz42
# magma@=2.6.2%gcc@=8.5.0+cuda+fortran~ipo~rocm+shared build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-centos8-power9le
module load magma/2.6.2-gcc-8.5.0-wcsaklt
# metis@=5.1.0%gcc@=8.5.0~gdb~int64~ipo~real64+shared build_system=cmake build_type=Release generator=make patches=4991da9,93a7903,b1225da arch=linux-centos8-power9le
module load metis/5.1.0-gcc-8.5.0-7u2ygpy
# openmpi@=4.1.4%gcc@=8.5.0~atomics~cuda~cxx~cxx_exceptions~gpfs~internal-hwloc~internal-libevent~internal-pmix~java~legacylaunchers~lustre~memchecker~openshmem~orterunprefix~romio+rsh~singularity~static+vt+wrapper-rpath build_system=autotools fabrics=none schedulers=none arch=linux-centos8-power9le
module load openmpi/4.1.4-gcc-8.5.0-ntdlotr
# raja@=0.14.0%gcc@=8.5.0+cuda~examples~exercises~ipo+openmp~plugins~rocm+shared~tests build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-centos8-power9le
module load raja/0.14.0-gcc-8.5.0-gw6xepf
# libsigsegv@=2.14%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load libsigsegv/2.14-gcc-8.5.0-zynce5o
# m4@=1.4.19%gcc@=8.5.0+sigsegv build_system=autotools patches=9dc5fbd,bfdffa7 arch=linux-centos8-power9le
module load m4/1.4.19-gcc-8.5.0-m4hoz3i
# autoconf@=2.72%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load autoconf/2.72-gcc-8.5.0-datymaf
# automake@=1.16.5%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load automake/1.16.5-gcc-8.5.0-xtl6kc7
# findutils@=4.9.0%gcc@=8.5.0 build_system=autotools patches=440b954 arch=linux-centos8-power9le
module load findutils/4.9.0-gcc-8.5.0-3vzp5dt
# libtool@=2.4.7%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load libtool/2.4.7-gcc-8.5.0-fbte22i
# gmp@=6.2.1%gcc@=8.5.0+cxx build_system=autotools libs=shared,static patches=69ad2e2 arch=linux-centos8-power9le
module load gmp/6.2.1-gcc-8.5.0-z2jr2nf
# autoconf-archive@=2023.02.20%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load autoconf-archive/2023.02.20-gcc-8.5.0-lba25r6
# xz@=5.4.1%gcc@=8.5.0~pic build_system=autotools libs=shared,static arch=linux-centos8-power9le
module load xz/5.4.1-gcc-8.5.0-655mase
# libxml2@=2.10.3%gcc@=8.5.0+pic~python+shared build_system=autotools arch=linux-centos8-power9le
module load libxml2/2.10.3-gcc-8.5.0-kt7uf6w
# pigz@=2.8%gcc@=8.5.0 build_system=makefile arch=linux-centos8-power9le
module load pigz/2.8-gcc-8.5.0-2nhvqh2
# zstd@=1.5.5%gcc@=8.5.0+programs build_system=makefile compression=none libs=shared,static arch=linux-centos8-power9le
module load zstd/1.5.5-gcc-8.5.0-l5rig24
# tar@=1.34%gcc@=8.5.0 build_system=autotools zip=pigz arch=linux-centos8-power9le
module load tar/1.34-gcc-8.5.0-ovvdvs3
# gettext@=0.22.4%gcc@=8.5.0+bzip2+curses+git~libunistring+libxml2+pic+shared+tar+xz build_system=autotools arch=linux-centos8-power9le
module load gettext/0.22.4-gcc-8.5.0-wacf7kj
# texinfo@=7.0.3%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load texinfo/7.0.3-gcc-8.5.0-3tatox3
# mpfr@=4.2.1%gcc@=8.5.0 build_system=autotools libs=shared,static arch=linux-centos8-power9le
module load mpfr/4.2.1-gcc-8.5.0-qp26atl
# suite-sparse@=5.13.0%gcc@=8.5.0~cuda~graphblas~openmp+pic build_system=generic arch=linux-centos8-power9le
module load suite-sparse/5.13.0-gcc-8.5.0-e2alwtx
# umpire@=6.0.0%gcc@=8.5.0+c+cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp~rocm~shared build_system=cmake build_type=Release cuda_arch=70 generator=make tests=none arch=linux-centos8-power9le
module load umpire/6.0.0-gcc-8.5.0-k2dmyxy
# hiop@=develop%gcc@=8.5.0+cuda~cusolver_lu~deepchecking+ginkgo~ipo~jsrun+kron+mpi+raja~rocm~shared+sparse build_system=cmake build_type=Release cuda_arch=70 dev_path=/people/svcearthshot/gitlab/30290/spack_newell/hiop_dev generator=make arch=linux-centos8-power9le
module load hiop/develop-gcc-8.5.0-rhaypd6
# ipopt@=3.12.10%gcc@=8.5.0+coinhsl~debug~metis~mumps build_system=autotools arch=linux-centos8-power9le
module load ipopt/3.12.10-gcc-8.5.0-7sqpcwh
# python@=3.8.5%gcc@=8.5.0+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=0d98e93,4c24573,ebdca64,f2fd060 arch=linux-centos8-power9le
module load python/3.8.5-gcc-8.5.0-3ohhmxo
# petsc@=3.20.3%gcc@=8.5.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~sycl~tetgen~trilinos~valgrind~zoltan build_system=generic clanguage=C memalign=none arch=linux-centos8-power9le
module load petsc/3.20.3-gcc-8.5.0-qxxu574
# py-pip@=23.1.2%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load py-pip/23.1.2-gcc-8.5.0-altra54
# py-setuptools@=68.0.0%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load py-setuptools/68.0.0-gcc-8.5.0-rndmqnv
# py-wheel@=0.41.2%gcc@=8.5.0 build_system=generic arch=linux-centos8-power9le
module load py-wheel/0.41.2-gcc-8.5.0-x26btdo
# py-cython@=0.29.36%gcc@=8.5.0 build_system=python_pip patches=c4369ad arch=linux-centos8-power9le
module load py-cython/0.29.36-gcc-8.5.0-nusktwp
# py-mpi4py@=3.1.5%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-mpi4py/3.1.5-gcc-8.5.0-rf7gwc6
# py-flit-core@=3.9.0%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-flit-core/3.9.0-gcc-8.5.0-c6bnetl
# libmd@=1.0.4%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load libmd/1.0.4-gcc-8.5.0-ghqr67q
# libbsd@=0.11.7%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load libbsd/0.11.7-gcc-8.5.0-ji4eh3o
# expat@=2.5.0%gcc@=8.5.0+libbsd build_system=autotools arch=linux-centos8-power9le
module load expat/2.5.0-gcc-8.5.0-udfpyas
# libunistring@=1.1%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load libunistring/1.1-gcc-8.5.0-epmig4y
# libidn2@=2.3.4%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load libidn2/2.3.4-gcc-8.5.0-iu6yp7z
# bison@=3.8.2%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load bison/3.8.2-gcc-8.5.0-d4pxux7
# krb5@=1.20.1%gcc@=8.5.0+shared build_system=autotools arch=linux-centos8-power9le
module load krb5/1.20.1-gcc-8.5.0-4toeerd
# libedit@=3.1-20210216%gcc@=8.5.0 build_system=autotools arch=linux-centos8-power9le
module load libedit/3.1-20210216-gcc-8.5.0-dbftmgn
# libxcrypt@=4.4.35%gcc@=8.5.0~obsolete_api build_system=autotools patches=4885da3 arch=linux-centos8-power9le
module load libxcrypt/4.4.35-gcc-8.5.0-34deseu
# openssh@=9.5p1%gcc@=8.5.0+gssapi build_system=autotools arch=linux-centos8-power9le
module load openssh/9.5p1-gcc-8.5.0-iypfxiy
# pcre2@=10.42%gcc@=8.5.0~jit+multibyte build_system=autotools arch=linux-centos8-power9le
module load pcre2/10.42-gcc-8.5.0-igublbw
# git@=2.42.0%gcc@=8.5.0+man+nls+perl+subtree~svn~tcltk build_system=autotools arch=linux-centos8-power9le
module load git/2.42.0-gcc-8.5.0-wnjkvsr
# py-packaging@=23.1%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-packaging/23.1-gcc-8.5.0-faorpri
# py-tomli@=2.0.1%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-tomli/2.0.1-gcc-8.5.0-zenze6c
# py-typing-extensions@=4.8.0%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-typing-extensions/4.8.0-gcc-8.5.0-d2bqqhv
# py-setuptools-scm@=7.1.0%gcc@=8.5.0+toml build_system=python_pip arch=linux-centos8-power9le
module load py-setuptools-scm/7.1.0-gcc-8.5.0-3tgsvnx
# py-flit-scm@=1.7.0%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-flit-scm/1.7.0-gcc-8.5.0-irmu4g4
# py-exceptiongroup@=1.1.1%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-exceptiongroup/1.1.1-gcc-8.5.0-iffozcl
# py-editables@=0.3%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-editables/0.3-gcc-8.5.0-i6vrh67
# py-pathspec@=0.11.1%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-pathspec/0.11.1-gcc-8.5.0-fenlcx7
# py-pluggy@=1.4.0%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-pluggy/1.4.0-gcc-8.5.0-b6rdncv
# py-calver@=2022.6.26%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-calver/2022.6.26-gcc-8.5.0-d6vrsrb
# py-trove-classifiers@=2023.8.7%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-trove-classifiers/2023.8.7-gcc-8.5.0-mwtxbac
# py-hatchling@=1.21.0%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-hatchling/1.21.0-gcc-8.5.0-dmwtut6
# py-hatch-vcs@=0.3.0%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-hatch-vcs/0.3.0-gcc-8.5.0-c5ejpb5
# py-iniconfig@=2.0.0%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-iniconfig/2.0.0-gcc-8.5.0-7tok6sm
# py-pytest@=8.0.0%gcc@=8.5.0 build_system=python_pip arch=linux-centos8-power9le
module load py-pytest/8.0.0-gcc-8.5.0-pttcqxm
# exago@=develop%gcc@=8.5.0+cuda+hiop~ipo+ipopt+logging+mpi+python+raja~rocm build_system=cmake build_type=Release cuda_arch=70 dev_path=/people/svcearthshot/gitlab/30290/spack_newell generator=make arch=linux-centos8-power9le
## module load exago/develop-gcc-8.5.0-4d5hzqt
