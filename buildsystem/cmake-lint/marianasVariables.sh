source /etc/profile.d/modules.sh
module purge
module load python/miniconda3.8
source /share/apps/python/miniconda3.8/etc/profile.d/conda.sh
conda activate /qfs/projects/exasgd/marianas/cmake-format-python-env
export srcdir=${srcdir:-$PWD}