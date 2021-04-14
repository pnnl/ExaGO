
# Developer Guidelines

## Contributing

### Summary

- Where a clear precedent exists in the project, follow it.
- Otherwise, where a good public C++ style guide is relevant and clear, follow it. [Google's](https://google.github.io/styleguide/cppguide.html) is pretty good and comes with lots of justifications for its rules. [Flang's](https://github.com/llvm/llvm-project/blob/main/flang/docs/C%2B%2Bstyle.md) is also very thorough and is heavily referenced here.

### Development Workflow

To get started contributing to ExaGO, check out a feature or fix branch off of
the `develop` branch:

```shell
$ cd exago
$ git checkout develop
$ git checkout -b my-cool-feature-dev
```

Feature branches must have the suffix `-dev` and fix branches `-fix`.

To submit a merge request, please make sure your branch is up to date with
current development like so:

```shell
$ git pull --rebase origin develop
```

This will ensure you are rebased on the most recent development branch.
If you see [open merge requests on ExaGO's repository](https://gitlab.pnnl.gov/exasgd/frameworks/exago/-/merge_requests) that touch the same lines of code that you are working on, please coordinate with the developers of that merge request so your work doesn't conflict.
Ideally, only one person is working on a single part of the codebase - this
prevents situations which require difficult rebases, which are extremely
error-prone.

If your commit only touches comments and/or documentation, please add `[skip ci]` to your commit message.
This reduces unneeded burden on our CI system.

### Commit Messages

Commit messages should use the imperative mood.
[This widely referenced article on commit messages](https://chris.beams.io/posts/git-commit/#imperative)
goes into greater detail and describes other best practices that are also
recommended.
*Merge request titles should also follow the same conventions.*

### Merge Request Conventions

- ***IMPORTANT: Patches should be rather small***. Extremely large merge requests are helpful as proofs-of-concept, but when merging into key development branches, large merge requests should be broken into smaller ones. Smaller merge requests are easier to review and test, and when something goes wrong it's much easier to track down a bug. Rarely is a patch *too small* to constitute it's own MR.
- Each MR should usually have a corresponding issue.

### Releases

Each month, we tag our master branch after incrementing the patch number.
We also push a release if we merge a feature significant enough to have its own
release.

### Files

- Header files should be idempotent.

### Coding Conventions

- Don't use dynamic solutions to solve problems that can be solved at build time; don't solve build time problems by writing programs that produce source code when macros and templates suffice; don't write macros when templates suffice.
- Indent with two spaces.
- Don't indent `public:`, `protected:`, and `private:` accessibility labels.
- Never use more than 80 characters per source line.
- Don't use tabs.
- Don't indent the bodies of namespaces, even when nested.
- Function result types go on the same line as the function and argument names.

### CMake

- Use target-oriented cmake.
- A `.cmake-format` file exists in the top-level directory - if you change any cmake code, please run `cmake-format -i <file name>` before pushing your changes. The build script may also apply these changes for you by using the job specified by `./buildsystem/build.sh --job=cmake-lint-apply`.

## References

1. [LLVM Flang's style guide](https://github.com/llvm/llvm-project/blob/main/flang/docs/C%2B%2Bstyle.md)
1. [Google's style guide](https://google.github.io/styleguide/cppguide.html)
1. [ExaGO's public repository](https://gitlab.pnnl.gov/exasgd/frameworks/exago/-/merge_requests)
1. [Chris Beams's blog post on writing commit messages](https://chris.beams.io/posts/git-commit)
1. [Pablo Ariasblog post on cmake best practices](https://pabloariasal.github.io/2018/02/19/its-time-to-do-cmake-right/)
