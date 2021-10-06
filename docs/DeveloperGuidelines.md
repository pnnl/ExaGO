
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

### Coding Conventions

- Use `PascalCase` for functions/classes, and `snake_case` for variables.
- Use underscores (`_`) and all lowercase characters for filenames.
- Prefer `enum class` over a bare `enum`.
- Use `constexpr`, `static`, and `const` wherever appropriate. It is good practice to make all variables and functions `static`, `const`, and/or `constexpr` by default, and only remove these qualifiers when you have good reason to do so.
- Use `auto` to declare a variable when the type is obvious. For example, when dynamic casting from a pointer, there is no need to specify the type twice:
```cpp
auto *my_var = dynamic_cast<MyType*>(my_other_var);
```
- Functions should usually be less than 50 lines of code. If they are longer, consider breaking the function down into smaller functions.
- All C++ APIs should be wrapped in the `ExaGO::` namespace, followed by the module name. For example, the `ExaGO::OPFLOW::` namespaces should wrap any public APIs in the `src/opflow` directory.
- Use anonymous namespaces for functions and variables which should have internal linkage.
- Don't use dynamic solutions to solve problems that can be solved at build time; don't solve build time problems by writing programs that produce source code when macros and templates suffice; don't write macros when templates suffice.
- Indent with two spaces.
- Don't indent `public:`, `protected:`, and `private:` accessibility labels.
- Never use more than 80 characters per source line, unless you're editing documentation.
- Don't use tabs.
- Don't indent the bodies of namespaces, even when nested.
- Function result types go on the same line as the function and argument names.

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
release. This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html) and the [Keep-a-Changelog](https://keepachangelog.com/en/1.0.0/) guidelines.

### Files

- Header files should be idempotent and use header guards.

### CMake

- Use target-oriented cmake.
- A `.cmake-format` file exists in the top-level directory - if you change any cmake code, please run `cmake-format -i <file name>` before pushing your changes. The build script may also apply these changes for you by using the job specified by `./buildsystem/build.sh --job=cmake-lint-apply`. See the build script's help message (`./buildsystem/build.sh --help`) for more information on this.

## References

1. [LLVM Flang's style guide](https://github.com/llvm/llvm-project/blob/main/flang/docs/C%2B%2Bstyle.md)
1. [Google's style guide](https://google.github.io/styleguide/cppguide.html)
1. [ExaGO's public repository](https://gitlab.pnnl.gov/exasgd/frameworks/exago/-/merge_requests)
1. [Chris Beams's blog post on writing commit messages](https://chris.beams.io/posts/git-commit)
1. [Pablo Ariasblog post on cmake best practices](https://pabloariasal.github.io/2018/02/19/its-time-to-do-cmake-right/)
1. [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
1. [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
1. [xSDK Project](https://xsdk.info/policies/)
