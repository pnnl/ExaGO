source /etc/profile.d/modules.sh
module purge
module load python
export PATH="$PATH:~/.local/bin"
export srcdir=${srcdir:-$PWD}
python -m venv python-env
source python-env/bin/activate
pip install cmake-format
