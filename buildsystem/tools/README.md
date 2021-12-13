## ExaGO Code Quality Tools

This directory contains tools used to manage the code quality of the ExaGO codebase.
Many of these tools enforce the developer guidelines in `docs/developer_guidelines.md`.

An individual tool is a submodule of the ExaGO perl module found in the `lib/` directory.
The tool drivers found in this directory run one or more tool from the `lib/ExaGO/` directory.

### Adding a Tool

To add a tool to the perl module, create a new file `lib/ExaGO/MyTool.pm` with the package name `package ExaGO::MyTool;`.

At the very least, your tool should have this header at the top of the module:
```perl
package ExaGO::MyTool;
use strict;
use warnings;
use Exporter;
our @ISA    = qw( Exporter );
our @EXPORT = qw(tool);
```
This will export the `tool` subroutine from your perl module.

The next step is to create the `tool` subroutine:
```perl
sub tool {
  say "Running MyTool!";
}
```

You may also create a driver for your tool in the top-level directory which runs your tool after taking some arguments.

If you make edits to a code quality tool in this directory, you should lint your code with the `perltidy` tool like so:

```console
find \
  buildsystem/tools -name '*.pl' -o -name '*.pm' \
  -exec perltidy -i=2 -b {} \;
```

You will have to install the `Perl::Tidy` module from CPAN to run `perltidy`.
