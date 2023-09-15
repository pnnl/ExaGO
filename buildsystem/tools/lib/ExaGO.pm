# Use `perltidy -i=2 -b` on all perl files to format this file after editing.
# Must install the perl module Perl::Tidy to format.
#
# Asher Mancinelli <asher.mancinelli@pnnl.gov>
#
# Apply cmake format to all source files in ExaGO, or check for formatting.

package ExaGO;

use strict;
use warnings;
use v5.16;

use ExaGO::Env;
use ExaGO::FileNamingConventions;

use Exporter;
our @ISA    = qw( Exporter );
our @EXPORT = qw($root printoneline $host $el);

1;
