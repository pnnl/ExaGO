module use -a /qfs/projects/exasgd/src/ci-deception/ci-modules/linux-centos7-zen2
# pkgconf@1.8.0%gcc@10.2.0 arch=linux-centos7-zen2
module load pkgconf-1.8.0-gcc-10.2.0-fuflwbl
# ncurses@6.3%gcc@10.2.0~symlinks+termlib abi=none arch=linux-centos7-zen2
module load ncurses-6.3-gcc-10.2.0-4wlnxto
# ca-certificates-mozilla@2022-07-19%gcc@10.2.0 arch=linux-centos7-zen2
module load ca-certificates-mozilla-2022-07-19-gcc-10.2.0-h2opehw
# berkeley-db@18.1.40%gcc@10.2.0+cxx~docs+stl patches=b231fcc arch=linux-centos7-zen2
module load berkeley-db-18.1.40-gcc-10.2.0-hltd4j3
# libiconv@1.16%gcc@10.2.0 libs=shared,static arch=linux-centos7-zen2
module load libiconv-1.16-gcc-10.2.0-gbg7l5p
# diffutils@3.8%gcc@10.2.0 arch=linux-centos7-zen2
module load diffutils-3.8-gcc-10.2.0-mjfwces
# bzip2@1.0.8%gcc@10.2.0~debug~pic+shared arch=linux-centos7-zen2
module load bzip2-1.0.8-gcc-10.2.0-bxh46iv
# readline@8.1.2%gcc@10.2.0 arch=linux-centos7-zen2
module load readline-8.1.2-gcc-10.2.0-vtya5ay
# gdbm@1.19%gcc@10.2.0 arch=linux-centos7-zen2
module load gdbm-1.19-gcc-10.2.0-efj5agg
# zlib@1.2.12%gcc@10.2.0+optimize+pic+shared patches=0d38234 arch=linux-centos7-zen2
module load zlib-1.2.12-gcc-10.2.0-gnkqokp
# perl@5.34.1%gcc@10.2.0+cpanm+shared+threads arch=linux-centos7-zen2
module load perl-5.34.1-gcc-10.2.0-xp4fpdr
# openssl@1.1.1q%gcc@10.2.0~docs~shared certs=mozilla patches=3fdcf2d arch=linux-centos7-zen2
## module load openssl-1.1.1q-gcc-10.2.0-xhxspos
# cmake@3.23.3%gcc@10.2.0~doc+ncurses+ownlibs~qt build_type=Release arch=linux-centos7-zen2
module load cmake-3.23.3-gcc-10.2.0-ggyj7bs
# blt@0.4.1%gcc@10.2.0 arch=linux-centos7-zen2
module load blt-0.4.1-gcc-10.2.0-oabae2w
# cub@1.16.0%gcc@10.2.0 arch=linux-centos7-zen2
module load cub-1.16.0-gcc-10.2.0-ovgrtom
# cuda@11.4%gcc@10.2.0~allow-unsupported-compilers~dev arch=linux-centos7-zen2
module load cuda-11.4-gcc-10.2.0-ewurpsv
# camp@0.2.3%gcc@10.2.0+cuda~ipo+openmp~rocm~tests build_type=Release cuda_arch=60 arch=linux-centos7-zen2
module load camp-0.2.3-gcc-10.2.0-w5sq7yo
# openblas@0.3.20%gcc@10.2.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared patches=9f12903 symbol_suffix=none threads=none arch=linux-centos7-zen2
module load openblas-0.3.20-gcc-10.2.0-x6v3mwm
# coinhsl@2019.05.21%gcc@10.2.0+blas arch=linux-centos7-zen2
module load coinhsl-2019.05.21-gcc-10.2.0-gkzkws6
# ginkgo@1.5.0.glu_experimental%gcc@10.2.0+cuda~develtools+full_optimizations~hwloc~ipo~oneapi+openmp~rocm+shared build_type=Release cuda_arch=60 arch=linux-centos7-zen2
module load ginkgo-1.5.0.glu_experimental-gcc-10.2.0-wcabkw2
# magma@2.6.2%gcc@10.2.0+cuda+fortran~ipo~rocm+shared build_type=Release cuda_arch=60 arch=linux-centos7-zen2
module load magma-2.6.2-gcc-10.2.0-sotxtwz
# metis@5.1.0%gcc@10.2.0~gdb~int64~real64+shared build_type=Release patches=4991da9,b1225da arch=linux-centos7-zen2
module load metis-5.1.0-gcc-10.2.0-k4z4v6l
# openmpi@4.1.0mlx5.0%gcc@10.2.0~atomics~cuda~cxx~cxx_exceptions~gpfs~internal-hwloc~java~legacylaunchers~lustre~memchecker+romio+rsh~singularity+static+vt+wrapper-rpath fabrics=none patches=60ce20b schedulers=none arch=linux-centos7-zen2
module load openmpi-4.1.0mlx5.0-gcc-10.2.0-ytj7jxb
# raja@0.14.0%gcc@10.2.0+cuda~examples~exercises~ipo+openmp~rocm+shared~tests build_type=Release cuda_arch=60 arch=linux-centos7-zen2
module load raja-0.14.0-gcc-10.2.0-h4zcwkw
# libsigsegv@2.13%gcc@10.2.0 arch=linux-centos7-zen2
module load libsigsegv-2.13-gcc-10.2.0-aj5goyi
# m4@1.4.19%gcc@10.2.0+sigsegv patches=9dc5fbd,bfdffa7 arch=linux-centos7-zen2
module load m4-1.4.19-gcc-10.2.0-k5kkyx6
# autoconf@2.69%gcc@10.2.0 patches=35c4492,7793209,a49dd5b arch=linux-centos7-zen2
module load autoconf-2.69-gcc-10.2.0-jnh4mbw
# automake@1.16.5%gcc@10.2.0 arch=linux-centos7-zen2
module load automake-1.16.5-gcc-10.2.0-pgpzgqq
# libtool@2.4.7%gcc@10.2.0 arch=linux-centos7-zen2
module load libtool-2.4.7-gcc-10.2.0-mzc2mvw
# gmp@6.2.1%gcc@10.2.0 libs=shared,static arch=linux-centos7-zen2
module load gmp-6.2.1-gcc-10.2.0-tpo7i4x
# autoconf-archive@2022.02.11%gcc@10.2.0 patches=139214f arch=linux-centos7-zen2
module load autoconf-archive-2022.02.11-gcc-10.2.0-tirhdzr
# texinfo@6.5%gcc@10.2.0 patches=12f6edb,1732115 arch=linux-centos7-zen2
module load texinfo-6.5-gcc-10.2.0-mcrbwnj
# mpfr@4.1.0%gcc@10.2.0 libs=shared,static arch=linux-centos7-zen2
module load mpfr-4.1.0-gcc-10.2.0-3yutkz3
# suite-sparse@5.10.1%gcc@10.2.0~cuda~graphblas~openmp+pic~tbb arch=linux-centos7-zen2
module load suite-sparse-5.10.1-gcc-10.2.0-add65sb
# umpire@6.0.0%gcc@10.2.0+c+cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp~rocm~shared build_type=Release cuda_arch=60 tests=none arch=linux-centos7-zen2
module load umpire-6.0.0-gcc-10.2.0-ae6qhla
# hiop@develop%gcc@10.2.0+cuda~cusolver~deepchecking+full_optimizations+ginkgo~ipo~jsrun+kron+mpi+raja~rocm~shared+sparse build_type=Release cuda_arch=60 arch=linux-centos7-zen2
module load hiop-develop-gcc-10.2.0-iuotcdz
# ipopt@3.12.10%gcc@10.2.0+coinhsl~debug~metis~mumps arch=linux-centos7-zen2
module load ipopt-3.12.10-gcc-10.2.0-c3chb6f
# python@3.9.12%gcc@10.2.0+bz2+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tix~tkinter~ucs4+uuid+zlib patches=0d98e93,4c24573,f2fd060 arch=linux-centos7-zen2
module load python-3.9.12-gcc-10.2.0-d4tck3b
# petsc@3.16.6%gcc@10.2.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~tetgen~trilinos~valgrind clanguage=C arch=linux-centos7-zen2
module load petsc-3.16.6-gcc-10.2.0-5gsowcj
# py-pip@22.1.2%gcc@10.2.0 arch=linux-centos7-zen2
module load py-pip-22.1.2-gcc-10.2.0-2kqxny4
# py-setuptools@65.0.0%gcc@10.2.0 arch=linux-centos7-zen2
module load py-setuptools-65.0.0-gcc-10.2.0-ohs4k5t
# py-wheel@0.37.1%gcc@10.2.0 arch=linux-centos7-zen2
module load py-wheel-0.37.1-gcc-10.2.0-q6kgedz
# py-mpi4py@3.1.2%gcc@10.2.0 arch=linux-centos7-zen2
module load py-mpi4py-3.1.2-gcc-10.2.0-odfcvnk
# py-attrs@21.4.0%gcc@10.2.0 arch=linux-centos7-zen2
module load py-attrs-21.4.0-gcc-10.2.0-w3og53i
# py-iniconfig@1.1.1%gcc@10.2.0 arch=linux-centos7-zen2
module load py-iniconfig-1.1.1-gcc-10.2.0-bz5lvsm
# py-pyparsing@3.0.6%gcc@10.2.0 arch=linux-centos7-zen2
module load py-pyparsing-3.0.6-gcc-10.2.0-hti3bid
# py-packaging@21.3%gcc@10.2.0 arch=linux-centos7-zen2
module load py-packaging-21.3-gcc-10.2.0-fvfwsrs
# py-tomli@1.2.2%gcc@10.2.0 arch=linux-centos7-zen2
module load py-tomli-1.2.2-gcc-10.2.0-2rcj7pu
# py-flit-core@3.6.0%gcc@10.2.0 arch=linux-centos7-zen2
module load py-flit-core-3.6.0-gcc-10.2.0-plq5vpq
# py-typing-extensions@4.3.0%gcc@10.2.0 arch=linux-centos7-zen2
module load py-typing-extensions-4.3.0-gcc-10.2.0-2wwhyww
# py-setuptools-scm@7.0.3%gcc@10.2.0+toml arch=linux-centos7-zen2
module load py-setuptools-scm-7.0.3-gcc-10.2.0-2xvjkec
# py-pluggy@1.0.0%gcc@10.2.0 arch=linux-centos7-zen2
module load py-pluggy-1.0.0-gcc-10.2.0-m2l6skg
# py-py@1.11.0%gcc@10.2.0 arch=linux-centos7-zen2
module load py-py-1.11.0-gcc-10.2.0-meh5son
# py-toml@0.10.2%gcc@10.2.0 arch=linux-centos7-zen2
module load py-toml-0.10.2-gcc-10.2.0-gqgcbhp
# py-pytest@6.2.5%gcc@10.2.0 arch=linux-centos7-zen2
module load py-pytest-6.2.5-gcc-10.2.0-qnxtpdk
# exago@develop%gcc@10.2.0+cuda+full_optimizations+hiop~ipo+ipopt+mpi+python+raja~rocm build_type=Release cuda_arch=60 dev_path=/people/svcexasgd/gitlab/17410/spack_deception arch=linux-centos7-zen2
## module load exago-develop-gcc-10.2.0-cuqll76
