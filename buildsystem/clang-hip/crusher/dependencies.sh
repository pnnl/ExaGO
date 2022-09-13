# blt@0.4.1%clang@14.0.0 arch=cray-sles15-zen3
module load blt-0.4.1-clang-14.0.0-eez6hwz
# camp@0.2.3%clang@14.0.0~cuda~ipo+rocm~tests amdgpu_target=gfx90a build_type=RelWithDebInfo arch=cray-sles15-zen3
module load camp-0.2.3-clang-14.0.0-gmhjhg7
# cmake@3.20.6%clang@14.0.0~doc+ncurses+ownlibs~qt build_type=Release arch=cray-sles15-zen3
module load cmake-3.20.6-clang-14.0.0-r3spc3q
# coinhsl@2019.05.21%clang@14.0.0+blas arch=cray-sles15-zen3
module load coinhsl-2019.05.21-clang-14.0.0-zk44j3w
# cray-mpich@8.1.16%gcc@11.2.0+wrappers arch=cray-sles15-zen3
module load cray-mpich-8.1.16-gcc-11.2.0-yfbzrpg
# diffutils@3.8%clang@14.0.0 arch=cray-sles15-zen3
module load diffutils-3.8-clang-14.0.0-xels4fd
# exago@git.crusher-dev%clang@14.0.0~cuda+hiop~ipo+ipopt+mpi~python+raja+rocm amdgpu_target=gfx90a build_type=RelWithDebInfo arch=cray-sles15-zen3
## module load exago-git.crusher-dev-clang-14.0.0-on6o47n
# ginkgo@glu_experimental%clang@14.0.0~cuda~develtools+full_optimizations~hwloc~ipo~oneapi+openmp+rocm+shared amdgpu_target=gfx90a build_type=Release arch=cray-sles15-zen3
module load ginkgo-glu_experimental-clang-14.0.0-i3pjtma
# hiop@develop%clang@14.0.0~cuda~cusolver~deepchecking+ginkgo~ipo~jsrun+kron+mpi+raja+rocm~shared+sparse amdgpu_target=gfx90a build_type=RelWithDebInfo arch=cray-sles15-zen3
module load hiop-develop-clang-14.0.0-dz3mvok
# hip@5.2.0%clang@14.0.0~ipo build_type=Release patches=959d1fe arch=cray-sles15-zen3
module load hip-5.2.0-clang-14.0.0-cvwypef
# hipblas@5.2.0%clang@14.0.0~ipo build_type=Release arch=cray-sles15-zen3
module load hipblas-5.2.0-clang-14.0.0-evzuomg
# hipsparse@5.2.0%clang@14.0.0~ipo build_type=Release arch=cray-sles15-zen3
module load hipsparse-5.2.0-clang-14.0.0-eth6np3
# hsa-rocr-dev@5.2.0%clang@14.0.0+image~ipo+shared build_type=Release patches=71e6851 arch=cray-sles15-zen3
module load hsa-rocr-dev-5.2.0-clang-14.0.0-mllpr5y
# ipopt@3.12.10%clang@14.0.0+coinhsl~debug~metis~mumps arch=cray-sles15-zen3
module load ipopt-3.12.10-clang-14.0.0-hiwa3zh
# libiconv@1.16%gcc@11.2.0 libs=shared,static arch=cray-sles15-zen3
module load libiconv-1.16-gcc-11.2.0-po4cjxv
# llvm-amdgpu@5.2.0%clang@14.0.0~ipo~link_llvm_dylib~llvm_dylib~openmp+rocm-device-libs build_type=Release patches=a08bbe1 arch=cray-sles15-zen3
module load llvm-amdgpu-5.2.0-clang-14.0.0-bzexqxu
# magma@2.6.2%clang@14.0.0~cuda+fortran~ipo+rocm+shared amdgpu_target=gfx90a build_type=RelWithDebInfo arch=cray-sles15-zen3
module load magma-2.6.2-clang-14.0.0-j5phv6v
# metis@5.1.0%clang@14.0.0~gdb~int64~real64+shared build_type=Release patches=4991da9 arch=cray-sles15-zen3
module load metis-5.1.0-clang-14.0.0-ua6eslw
# openblas@0.3.20%clang@14.0.0~bignuma~consistent_fpcsr~ilp64+locking+pic+shared symbol_suffix=none threads=openmp arch=cray-sles15-zen3
module load openblas-0.3.20-clang-14.0.0-4vd4qti
# parmetis@4.0.3%clang@14.0.0~gdb~int64~ipo+shared build_type=RelWithDebInfo patches=4f89253,50ed208,704b84f arch=cray-sles15-zen3
module load parmetis-4.0.3-clang-14.0.0-r4nxjdr
# perl@5.34.0%clang@14.0.0+cpanm+shared+threads arch=cray-sles15-zen3
module load perl-5.34.0-clang-14.0.0-mdatxiw
# petsc@3.16.5%clang@14.0.0~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind+metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~tetgen~trilinos~valgrind clanguage=C arch=cray-sles15-zen3
module load petsc-3.16.5-clang-14.0.0-5bqymj4
# pkgconf@1.8.0%clang@14.0.0 arch=cray-sles15-zen3
module load pkgconf-1.8.0-clang-14.0.0-67lrujf
# python@3.9.12%clang@14.0.0+bz2+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tix~tkinter~ucs4+uuid+zlib patches=0d98e93,4c24573,ebdca64,f2fd060 arch=cray-sles15-zen3
module load python-3.9.12-clang-14.0.0-wi6ggso
# raja@0.14.0%clang@14.0.0~cuda~examples~exercises~ipo~openmp+rocm+shared~tests amdgpu_target=gfx90a build_type=RelWithDebInfo arch=cray-sles15-zen3
module load raja-0.14.0-clang-14.0.0-lpvqr6p
# rocrand@5.2.0%clang@14.0.0~ipo amdgpu_target=auto build_type=Release arch=cray-sles15-zen3
module load rocrand-5.2.0-clang-14.0.0-lvitlo3
# rocthrust@5.2.0%clang@14.0.0~ipo build_type=Release arch=cray-sles15-zen3
module load rocthrust-5.2.0-clang-14.0.0-x5cspi3
# suite-sparse@4.5.6%clang@14.0.0~cuda~graphblas~openmp+pic~tbb arch=cray-sles15-zen3
module load suite-sparse-4.5.6-clang-14.0.0-kn37dfp
# umpire@6.0.0%clang@14.0.0+c~cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp+rocm+shared amdgpu_target=gfx90a build_type=RelWithDebInfo tests=none arch=cray-sles15-zen3
module load umpire-6.0.0-clang-14.0.0-3cgj2f3
