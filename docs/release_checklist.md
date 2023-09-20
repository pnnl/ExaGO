## Release Checklist

This document outlines the required checks performed before and after a release of ExaGO.

### Testing

Ensure all continuous integration pipelines are passing. 
As per P025, all CI pipelines should pass before an MR is merged.
As per P029, ensure Spack environments are re-built in CI semi frequently, especially before a release.
This is especially important when merging `develop` -> `master`.
Ensure all pipelines are passing before a release.

### Documentation

#### Changelog

As per ExaGO developer policy P024, this project adheres to the [Keep-a-Changelog](https://keepachangelog.com/en/1.0.0/) guidelines.
If the changelog has not been updated, review merged MRs and update the changelog accordingly.
It is especially important to document changes in versions of our dependencies in the changelog.

#### Manual

Ensure the user manual has the appropriate version and has been updated with any API changes.

### Tag

This projects aims to adhere to [version 2.0.0 of the Semantic Versioning guidelines](https://semver.org/spec/v2.0.0.html).
The most important details of the guideslines are included here:

> 1. MAJOR version when you make incompatible API changes,
> 1. MINOR version when you add functionality in a backwards compatible manner, and
> 1. PATCH version when you make backwards compatible bug fixes.
> 
> Additional labels for pre-release and build metadata are available as extensions to the MAJOR.MINOR.PATCH format.

Tags should only ever be made on the master branch.
Before tagging, `CHANGELOG.md` must be updated such that the items under the `develop` or `unrealeased` header are moved to a new header indicative of the new version being tagged.

### CMake version

Before merging `develop`->`master` ensure that the version number in the project definition in main [CMakeLists.txt](./CMakeLists.txt) is upated in accordance with the planned upcoming release.

### Spack

Update the Spack package for ExaGO.
Modify the ExaGO spack package and build the newest version on at least one CI or target platform, and ensure that all went smoothly.
Spack PRs sometimes take a long time to get reviewed, so it's important to get the PR opened quickly.

Once a release is made, and the Spack PR is merged, make sure to update the ExaGO spack submodule and rebuild CI per P029.
