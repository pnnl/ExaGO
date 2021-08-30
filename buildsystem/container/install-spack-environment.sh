git clone --branch ${SPACK_VERSION} ${SPACK_REPO} ${SPACK_INSTALL_LOCATION}
. "${SPACK_INSTALL_LOCATION}/share/spack/setup-env.sh"
spack env activate .

# Find a few system packages to improve build times
for pkg in cmake gmake python git; do
  spack external find $pkg
done
spack compiler find

# Show the full concretization
spack concretize -f
spack install --source --keep-stage -j `nproc`

spack test exago
