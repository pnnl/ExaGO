source /etc/profile.d/modules.sh
module purge
module load python
set -x
pip install cmake-format --user
export PATH="$PATH:~/.local/bin"
export srcdir=${$srcdir:-$PWD}
srcdir=${srcdir:-$PWD}
