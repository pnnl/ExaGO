module use -a /lustre/orion/csc359/proj-shared/nkouk/spack-install/modules/linux-sles15-zen3
# pkgconf@=1.9.5%clang@=16.0.0-rocm5.6.0-mixed build_system=autotools arch=linux-sles15-zen3
module load pkgconf/1.9.5-clang-16.0.0-rocm5.6.0-mixed-vp3ll6l
# nghttp2@=1.52.0%clang@=16.0.0-rocm5.6.0-mixed build_system=autotools arch=linux-sles15-zen3
module load nghttp2/1.52.0-clang-16.0.0-rocm5.6.0-mixed-prks6rd
# ca-certificates-mozilla@=2023-05-30%clang@=16.0.0-rocm5.6.0-mixed build_system=generic arch=linux-sles15-zen3
module load ca-certificates-mozilla/2023-05-30-clang-16.0.0-rocm5.6.0-mixed-swebvgy
# perl@=5.34.0%clang@=16.0.0-rocm5.6.0-mixed+cpanm+opcode+open+shared+threads build_system=generic arch=linux-sles15-zen3
module load perl/5.34.0-clang-16.0.0-rocm5.6.0-mixed-hyokfvp
# zlib-ng@=2.1.3%clang@=16.0.0-rocm5.6.0-mixed+compat+opt build_system=autotools patches=299b958,ae9077a,b692621 arch=linux-sles15-zen3
module load zlib-ng/2.1.3-clang-16.0.0-rocm5.6.0-mixed-kjg4hdy
# openssl@=3.1.3%clang@=16.0.0-rocm5.6.0-mixed~docs+shared build_system=generic certs=mozilla arch=linux-sles15-zen3
## module load openssl/3.1.3-clang-16.0.0-rocm5.6.0-mixed-q3rmy4x
# curl@=8.1.2%clang@=16.0.0-rocm5.6.0-mixed~gssapi~ldap~libidn2~librtmp~libssh~libssh2+nghttp2 build_system=autotools libs=shared,static tls=openssl arch=linux-sles15-zen3
module load curl/8.1.2-clang-16.0.0-rocm5.6.0-mixed-hpeeijp
# ncurses@=6.4%clang@=16.0.0-rocm5.6.0-mixed~symlinks+termlib abi=none build_system=autotools arch=linux-sles15-zen3
module load ncurses/6.4-clang-16.0.0-rocm5.6.0-mixed-j2aesfx
# cmake@=3.20.6%clang@=16.0.0-rocm5.6.0-mixed~doc+ncurses+ownlibs build_system=generic build_type=Release arch=linux-sles15-zen3
module load cmake/3.20.6-clang-16.0.0-rocm5.6.0-mixed-xcb2k24
# blt@=0.4.1%clang@=16.0.0-rocm5.6.0-mixed build_system=generic arch=linux-sles15-zen3
module load blt/0.4.1-clang-16.0.0-rocm5.6.0-mixed-wcabdav
# gmake@=4.4.1%clang@=16.0.0-rocm5.6.0-mixed~guile build_system=autotools arch=linux-sles15-zen3
module load gmake/4.4.1-clang-16.0.0-rocm5.6.0-mixed-jtzohhr
# hip@=5.6.0%clang@=16.0.0-rocm5.6.0-mixed~cuda+rocm build_system=cmake build_type=Release generator=make patches=aee7249,c2ee21c,e73e91b arch=linux-sles15-zen3
module load hip/5.6.0-clang-16.0.0-rocm5.6.0-mixed-lvcatvl
# hsa-rocr-dev@=5.6.0%clang@=16.0.0-rocm5.6.0-mixed+image+shared build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load hsa-rocr-dev/5.6.0-clang-16.0.0-rocm5.6.0-mixed-26mgwl7
# llvm-amdgpu@=5.6.0%clang@=16.0.0-rocm5.6.0-mixed~link_llvm_dylib~llvm_dylib~openmp+rocm-device-libs build_system=cmake build_type=Release generator=ninja patches=a08bbe1,c4750bb,d35aec9 arch=linux-sles15-zen3
module load llvm-amdgpu/5.6.0-clang-16.0.0-rocm5.6.0-mixed-tydmudc
# camp@=0.2.3%clang@=16.0.0-rocm5.6.0-mixed~cuda~ipo~openmp+rocm~tests amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load camp/0.2.3-clang-16.0.0-rocm5.6.0-mixed-sfn2hme
# cray-mpich@=8.1.25%clang@=16.0.0-rocm5.6.0-mixed+wrappers build_system=generic arch=linux-sles15-zen3
module load cray-mpich/8.1.25-clang-16.0.0-rocm5.6.0-mixed-6teohk2
# openblas@=0.3.20%gcc@=12.2.0-mixed~bignuma~consistent_fpcsr~ilp64+locking+pic+shared build_system=makefile patches=9f12903 symbol_suffix=none threads=none arch=linux-sles15-zen3
module load openblas/0.3.20-gcc-12.2.0-mixed-7hydqmq
# coinhsl@=2019.05.21%gcc@=12.2.0-mixed+blas build_system=autotools arch=linux-sles15-zen3
module load coinhsl/2019.05.21-gcc-12.2.0-mixed-pdovwcj
# hipblas@=5.6.0%clang@=16.0.0-rocm5.6.0-mixed~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load hipblas/5.6.0-clang-16.0.0-rocm5.6.0-mixed-oeelt53
# hipsparse@=5.6.0%clang@=16.0.0-rocm5.6.0-mixed~cuda+rocm amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load hipsparse/5.6.0-clang-16.0.0-rocm5.6.0-mixed-upordfk
# magma@=2.6.2%clang@=16.0.0-rocm5.6.0-mixed~cuda+fortran~ipo+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load magma/2.6.2-clang-16.0.0-rocm5.6.0-mixed-x2tryz7
# metis@=5.1.0%clang@=16.0.0-rocm5.6.0-mixed~gdb~int64~ipo~real64+shared build_system=cmake build_type=Release generator=make patches=4991da9,93a7903 arch=linux-sles15-zen3
module load metis/5.1.0-clang-16.0.0-rocm5.6.0-mixed-kzayikt
# rocprim@=5.6.0%clang@=16.0.0-rocm5.6.0-mixed amdgpu_target=auto build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load rocprim/5.6.0-clang-16.0.0-rocm5.6.0-mixed-ulpi7by
# raja@=0.14.0%clang@=16.0.0-rocm5.6.0-mixed~cuda~examples~exercises~ipo~openmp+rocm+shared~tests amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load raja/0.14.0-clang-16.0.0-rocm5.6.0-mixed-3t3z5vf
# suite-sparse@=4.5.6%clang@=16.0.0-rocm5.6.0-mixed~cuda~graphblas~openmp+pic~tbb build_system=generic arch=linux-sles15-zen3
module load suite-sparse/4.5.6-clang-16.0.0-rocm5.6.0-mixed-ercyneq
# umpire@=6.0.0%clang@=16.0.0-rocm5.6.0-mixed+c~cuda~device_alloc~deviceconst~examples~fortran~ipo~numa~openmp+rocm+shared amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make tests=none arch=linux-sles15-zen3
module load umpire/6.0.0-clang-16.0.0-rocm5.6.0-mixed-hre5hhv
# hiop@=0.7.2%clang@=16.0.0-rocm5.6.0-mixed~cuda~deepchecking~ginkgo~ipo~jsrun+kron+mpi+raja+rocm~shared+sparse amdgpu_target=gfx90a build_system=cmake build_type=Release generator=make arch=linux-sles15-zen3
module load hiop/0.7.2-clang-16.0.0-rocm5.6.0-mixed-i6l3zkn
# ipopt@=3.12.10%clang@=16.0.0-rocm5.6.0-mixed+coinhsl~debug~metis~mumps build_system=autotools arch=linux-sles15-zen3
module load ipopt/3.12.10-clang-16.0.0-rocm5.6.0-mixed-5qwrkim
# libiconv@=1.17%clang@=14.0.0-rocm5.2.0-mixed build_system=autotools libs=shared,static arch=linux-sles15-zen3
module load libiconv/1.17-clang-14.0.0-rocm5.2.0-mixed-ulmz25p
# diffutils@=3.9%clang@=14.0.0-rocm5.2.0-mixed build_system=autotools arch=linux-sles15-zen3
module load diffutils/3.9-clang-14.0.0-rocm5.2.0-mixed-scfhpum
# python@=3.9.12%clang@=14.0.0-rocm5.2.0-mixed+bz2+crypt+ctypes+dbm~debug+libxml2+lzma~nis~optimizations+pic+pyexpat+pythoncmd+readline+shared+sqlite3+ssl~tkinter+uuid+zlib build_system=generic patches=0d98e93,4c24573,ebdca64,f2fd060 arch=linux-sles15-zen3
module load python/3.9.12-clang-14.0.0-rocm5.2.0-mixed-xgn77vb
# petsc@=3.19.6%clang@=14.0.0-rocm5.2.0-mixed~X~batch~cgns~complex~cuda~debug+double~exodusii~fftw+fortran~giflib~hdf5~hpddm~hwloc~hypre~int64~jpeg~knl~kokkos~libpng~libyaml~memkind~metis~mkl-pardiso~mmg~moab~mpfr+mpi~mumps~openmp~p4est~parmmg~ptscotch~random123~rocm~saws~scalapack+shared~strumpack~suite-sparse~superlu-dist~tetgen~trilinos~valgrind build_system=generic clanguage=C memalign=none arch=linux-sles15-zen3
module load petsc/3.19.6-clang-14.0.0-rocm5.2.0-mixed-7dso2zm
# exago@=develop%clang@=16.0.0-rocm5.6.0-mixed~cuda+hiop~ipo+ipopt~logging+mpi~python+raja+rocm amdgpu_target=gfx90a build_system=cmake build_type=Release dev_path=/lustre/orion/scratch/nkouk/csc359/exago-frontier-amd-gfortran-github generator=make arch=linux-sles15-zen3
## module load exago/develop-clang-16.0.0-rocm5.6.0-mixed-7zsfkec
