# GitHub Actions Documentation

## `file_naming.yaml`

Runs on push, this action runs our perl script in `buildsystem/tools/file_naming_conventions.pl` to enforce certain name restrictions defined in our developer guidelines (P003). We use one action to checkout the code, then another to install perl before running the script.

## `ornl_ascent_mirror.yaml`

This pushes to a GitLab at ORNL and runs CI/CD on Ascent. This also can re-build modules for testing newer versions of ExaGO and it's dependencies without needing to monitor builds by hand.

## `pnnl_mirror.yaml`

Similar to `ornl_ascent_mirror.yaml`, this mirrors to PNNL GitLab, but also supports Incline, Decpeption and Newell.

## `pre_commit.yaml`

This enforces and runs pre-commit, and automatically commits fixes to any tests that it can. Noteably applies clang formatting, and cmake formatting most often, and requires developers to either install locally, or rebase to incorporate changes.

## `spack_cpu_build.yaml`

Logically does the following:
- Build a base image with Linux deps for things like mpi / gcc that are runtime dependencies for exago
- Build binaries in ubuntu GitHub actions runner for a matrix of ExaGO configs, force pushing each time to refresh binaries
- Push said binaries with the custom base image to the ghcr packages for exago

This also leverages the spack public mirror binaries as well as ExaGO's GHCR binaries, to drop builds down from taking 1.5hrs down to <10 minutes!! Now the slowest part is the concretization...

To pull the binaries and run, you can consult some more verbose docs:
- From GitHub https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry
- From Spack https://spack.readthedocs.io/en/latest/binary_caches.html#oci-docker-v2-registries-as-build-cache
