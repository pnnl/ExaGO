module use -a /gpfs/wolf/proj-shared/csc359/cameron/spack-install/ascent-modules/linux-rhel8-power9le
# cmake@=3.22.2%gcc@=10.2.0~doc+ncurses+ownlibs build_system=generic build_type=Release arch=linux-rhel8-power9le
module load cmake/3.22.2-gcc-10.2.0-k2nj5og
# blt@=0.4.1%gcc@=10.2.0 build_system=generic arch=linux-rhel8-power9le
module load blt/0.4.1-gcc-10.2.0-ffe5pzm
# cub@=2.1.0%gcc@=10.2.0 build_system=generic arch=linux-rhel8-power9le
module load cub/2.1.0-gcc-10.2.0-wqobtj5
# cuda@=11.4.2%gcc@=10.2.0~allow-unsupported-compilers~dev build_system=generic arch=linux-rhel8-power9le
module load cuda/11.4.2-gcc-10.2.0-mkeyz4o
# gnuconfig@=2022-09-17%gcc@=10.2.0 build_system=generic arch=linux-rhel8-power9le
module load gnuconfig/2022-09-17-gcc-10.2.0-errn6vw
# gmake@=4.4.1%gcc@=10.2.0~guile build_system=autotools arch=linux-rhel8-power9le
module load gmake/4.4.1-gcc-10.2.0-h6hscfq
# camp@=0.2.3%gcc@=10.2.0+cuda~ipo+openmp~rocm~tests build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-rhel8-power9le
module load camp/0.2.3-gcc-10.2.0-o4mm4xo
# perl@=5.30.1%gcc@=10.2.0+cpanm+opcode+open+shared+threads build_system=generic arch=linux-rhel8-power9le
module load perl/5.30.1-gcc-10.2.0-j665lil
# openblas@=0.3.23%gcc@=10.2.0~bignuma~consistent_fpcsr+fortran~ilp64+locking+pic+shared build_system=makefile symbol_suffix=none threads=none arch=linux-rhel8-power9le
module load openblas/0.3.23-gcc-10.2.0-55qdzyk
# magma@=2.7.2%gcc@=10.2.0+cuda+fortran~ipo~rocm+shared build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-rhel8-power9le
module load magma/2.7.2-gcc-10.2.0-5y555ji
# raja@=0.14.0%gcc@=10.2.0+cuda+examples+exercises~ipo+openmp~rocm+shared~tests build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-rhel8-power9le
module load raja/0.14.0-gcc-10.2.0-hqvuwn5
# spectrum-mpi@=10.4.0.3-20210112%gcc@=10.2.0 build_system=bundle arch=linux-rhel8-power9le
module load spectrum-mpi/10.4.0.3-20210112-gcc-10.2.0-tzrotbp
# umpire@=6.0.0%gcc@=10.2.0~c+cuda~device_alloc~deviceconst+examples~fortran~ipo~numa~openmp~rocm~shared build_system=cmake build_type=Release cuda_arch=70 generator=make tests=none arch=linux-rhel8-power9le
module load umpire/6.0.0-gcc-10.2.0-lhsel3h
# hiop@=1.0.0%gcc@=10.2.0+cuda~cusolver_lu~deepchecking~ginkgo~ipo~jsrun~kron+mpi+raja~rocm~shared~sparse build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-rhel8-power9le
module load hiop/1.0.0-gcc-10.2.0-g52fx4j
# coinhsl@=2015.06.23%gcc@=10.2.0+blas build_system=autotools arch=linux-rhel8-power9le
module load coinhsl/2015.06.23-gcc-10.2.0-aavytsk
# metis@=5.1.0%gcc@=10.2.0~gdb~int64~ipo~real64+shared build_system=cmake build_type=Release generator=make patches=4991da9,93a7903,b1225da arch=linux-rhel8-power9le
module load metis/5.1.0-gcc-10.2.0-wpcxkuq
# pkgconf@=1.9.5%gcc@=10.2.0 build_system=autotools arch=linux-rhel8-power9le
module load pkgconf/1.9.5-gcc-10.2.0-kisdb24
# ipopt@=3.12.10%gcc@=10.2.0+coinhsl~debug+metis~mumps build_system=autotools arch=linux-rhel8-power9le
module load ipopt/3.12.10-gcc-10.2.0-cbreqjj
# libiconv@=1.17%gcc@=10.2.0 build_system=autotools libs=shared,static arch=linux-rhel8-power9le
module load libiconv/1.17-gcc-10.2.0-ocxrxp4
# diffutils@=3.9%gcc@=10.2.0 build_system=autotools arch=linux-rhel8-power9le
module load diffutils/3.9-gcc-10.2.0-7e3rocg
# parmetis@=4.0.3%gcc@=10.2.0~gdb~int64~ipo+shared build_system=cmake build_type=Release generator=make patches=4f89253,50ed208,704b84f arch=linux-rhel8-power9le
module load parmetis/4.0.3-gcc-10.2.0-fkeq46a
# python@=3.9.7%gcc@=10.2.0+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=0d98e93,4c24573,f2fd060 arch=linux-rhel8-power9le
module load python/3.9.7-gcc-10.2.0-ej2ff3b
# petsc@=3.18.1%gcc@=10.2.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind+metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~tetgen~trilinos~valgrind build_system=generic clanguage=C memalign=none arch=linux-rhel8-power9le
module load petsc/3.18.1-gcc-10.2.0-vf4cchc
# exago@=develop%gcc@=10.2.0+cuda+hiop~ipo+ipopt+logging+mpi~python+raja~rocm build_system=cmake build_type=Release cuda_arch=70 generator=make arch=linux-rhel8-power9le
## module load exago/develop-gcc-10.2.0-rc6eoxa
