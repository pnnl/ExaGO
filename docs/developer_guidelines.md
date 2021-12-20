
# Developer Guidelines

## Contributing

### Summary

- Where an ExaGO policy exists, follow it.
- Otherwise, where a clear precedent exists in the project, follow it.
- Otherwise, where a good public C++ style guide is relevant and clear, follow it. [Google's](https://google.github.io/styleguide/cppguide.html) is very good and comes with justification for its rules. [Flang's](https://github.com/llvm/llvm-project/blob/main/flang/docs/C%2B%2Bstyle.md) is also very thorough and is heavily referenced here.
- Policies specific to ExaGO are numbered so they may be referenced elsewhere. If someone or a tool refers to ExaGO developer policy P157, that refers to policy number 157 in this document.

### Development Workflow

#### P022: Branch Off Of `develop` Unless You Have a Good Reason Not To

To get started contributing to ExaGO, check out a feature or fix branch off of
the `develop` branch:

```shell
$ cd exago
$ git checkout develop
$ git checkout -b my-cool-feature-dev
```

#### P023: Keep MRs Up To Date With Target Branch

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

#### P001: Don't Run Continuous Integration Pipelines If You Only Changed Documentation

If your commit only touches comments and/or documentation, please add `[skip ci]` to your commit message.
This reduces unneeded burden on our CI system.

### Coding Conventions

#### P002: Variable, Function, and Class Names

Use `PascalCase` for functions/classes, and `snake_case` for variables.

#### P004: Enums

Prefer `enum class` over a bare `enum`.

#### P005: Const, Static, and Constexpr Qualifiers

Use `constexpr`, `static`, and `const` wherever appropriate. It is good practice to make all variables and functions `static`, `const`, and/or `constexpr` by default, and only remove these qualifiers when you have good reason to do so.

#### P006: Auto

Use `auto` to declare a variable when the type is obvious. For example, when dynamic casting from a pointer, there is no need to specify the type twice:
```cpp
auto *my_var = dynamic_cast<MyType*>(my_other_var);
```

#### P007: Function Length

Functions should usually be less than 50 lines of code. If they are longer, consider breaking the function down into smaller functions.

#### P008: Anonymous Namespaces

Use anonymous namespaces for functions and variables which should have internal linkage.

#### P009: Perform as Much Work as Possible Before Runtime

Don't use dynamic solutions to solve problems that can be solved at build time; don't solve build time problems by writing programs that produce source code when macros and templates suffice; don't write macros when templates suffice.

#### P010: Line Length

Never use more than 80 characters per source line, unless you're editing documentation.

Use `clang-format` (or the utility script `buildsystem/tools/clang-format.pl` with the `-i` flag) to format your code.

Example usage of the utility script:
```console
$ cd exago/build
$ # do some development
$ make && ctest

$ # format ExaGO source in-place
$ ../buildsystem/tools/all.pl -i

$ # check to make sure everything is formatted (no flags means check if already formatted):
$ ../buildsystem/tools/all.pl
```

This will format all the C++ source code in the ExaGO source tree.

### Commit Messages

#### P011: Commit messages should use the imperative mood.

[This widely referenced article on commit messages](https://chris.beams.io/posts/git-commit/#imperative)
goes into greater detail and describes other best practices that are also
recommended.
*Merge request titles should also follow the same conventions.*

### Merge Request Conventions

#### P025: CI and Testing

If your merge request changes only documentation, you may add `[skip ci]` to your commit message to skip CI.
Otherwise, your MR must pass our continuous integration pipelines before being merged.

#### P012: Patches Should Be Small

***IMPORTANT: Patches should be rather small***.
Extremely large merge requests are helpful as proofs-of-concept, but when merging into key development branches, large merge requests should be broken into smaller ones. Smaller merge requests are easier to review and test, and when something goes wrong it's much easier to track down a bug. Rarely is a patch *too small* to constitute it's own MR.
A notable exception to this rule is when merging key development branches back into develop and/or master. It's hard to avoid large MRs in this case.

#### P013: MRs Resolve Issues

Each MR should usually have a corresponding issue.

### Releases

#### P014: Follow the Release Checklist

[The release checklist is linked here. Refer to this document for release documentation.](/docs/release_checklist.md)

#### P024: Changelog

This project adheres to the [Keep-a-Changelog](https://keepachangelog.com/en/1.0.0/) guidelines.
When making significant changes, summarize the changes in `CHANGELOG.md` under the `develop` or `unreleased` header.

### Files

#### P017: Header Idempotency

Header files should be idempotent.

#### P018: Header Guards

Header files should use header guards.

#### P003: File Naming Conventions

Use underscores (`_`) and all lowercase characters for C, C++, and Python filenames and paths.
Use pascal case for CMake code.
Use the `.h` extension for all C and C++ headers.
Use the `.cpp` extension for all C++ sources.

#### P021: Driver Naming Conventions

All binaries, applications, drivers, tests, and libraries should use underscores (`_`) and all lowercase characters.

### CMake

#### P019: Target-Oriented CMake

Use target-oriented cmake wherever possible.

#### P020: CMake Formatting

A `.cmake-format.py` file exists in the top-level directory - if you change any cmake code, please run `cmake-format -i <file name>` before pushing your changes.
You may also run `buildsystem/tools/cmake_format.pl -i`.

## References

1. [LLVM Flang's style guide](https://github.com/llvm/llvm-project/blob/main/flang/docs/C%2B%2Bstyle.md)
1. [Google's style guide](https://google.github.io/styleguide/cppguide.html)
1. [ExaGO's public repository](https://gitlab.pnnl.gov/exasgd/frameworks/exago/-/merge_requests)
1. [Chris Beams's blog post on writing commit messages](https://chris.beams.io/posts/git-commit)
1. [Pablo Ariasblog post on cmake best practices](https://pabloariasal.github.io/2018/02/19/its-time-to-do-cmake-right/)
1. [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
1. [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
1. [xSDK Project](https://xsdk.info/policies/)
