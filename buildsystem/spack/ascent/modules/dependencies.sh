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
# hiop@=develop%gcc@=11.2.0+cuda+cusolver_lu~deepchecking+ginkgo~ipo~jsrun+kron+mpi+raja~rocm~shared+sparse build_system=cmake build_type=MinSizeRel cuda_arch=70 dev_path=/gpfs/wolf/proj-shared/csc359/ci/466959/hiop generator=make arch=linux-rhel8-power9le
module load hiop/develop-gcc-11.2.0-rcgwvnh
# ipopt@=3.12.10%gcc@=11.2.0+coinhsl~debug~metis~mumps build_system=autotools arch=linux-rhel8-power9le
module load ipopt/3.12.10-gcc-11.2.0-7fkxunc
# python@=3.9.7%gcc@=11.2.0+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=0d98e93,4c24573,ebdca64,f2fd060 arch=linux-rhel8-power9le
module load python/3.9.7-gcc-11.2.0-gruthp3
# petsc@=3.19.6%gcc@=11.2.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~tetgen~trilinos~valgrind build_system=generic clanguage=C memalign=none arch=linux-rhel8-power9le
module load petsc/3.19.6-gcc-11.2.0-qj3xkii
# exago@=develop%gcc@=11.2.0+cuda+hiop~ipo+ipopt~logging+mpi~python+raja~rocm build_system=cmake build_type=Debug cuda_arch=70 dev_path=/gpfs/wolf/proj-shared/csc359/ci/466959 generator=make arch=linux-rhel8-power9le
## module load exago/develop-gcc-11.2.0-uweodnf
