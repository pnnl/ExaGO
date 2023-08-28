module use -a /qfs/projects/exasgd/src/ci-incline/install/modules/linux-centos7-zen
# exago@=develop%clang@=15.0.0-rocm5.3.0 cxxflags="--gcc-toolchain=/share/apps/gcc/8.4.0/" ~cuda+hiop~ipo+ipopt+logging+mpi~python+raja+rocm amdgpu_target=gfx908 build_system=cmake build_type=Release dev_path=/people/svcexasgd/gitlab/22333/spack_incline generator=make arch=linux-centos7-zen
module load exago/develop-clang-15.0.0-rocm5.3.0-z2sfasr
