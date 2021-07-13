# Python Wrapper

## Environment
Currently exago requires python version greater than [3.6](https://gitlab.pnnl.gov/exasgd/frameworks/exago/-/blob/python-wrapper-rebased-dev/CMakeLists.txt#L73). Once Python is configured (`module load python`), you will need to the ExaGO python libraries to your `PYTHONPATH` according to your target python version:
```bash
# For Python 3.6
export PYTHONPATH=/<exago_install>/lib/python3.6/site-packages:$PYTHONPATH
# For Python 3.7
export PYTHONPATH=/<exago_install>/lib/python3.7/site-packages:$PYTHONPATH
```
For ExaGO installations in non-standard locations, the `LD_LIBRARY_PATH` (`DYLD_LIBRARY_PATH` for MacOS) will also need to point to the ExaGO install location:
```bash
# For MacOS
# export DYLD_LIBRARY_PATH=/<exago_install>/lib:$DYLD_LIBRARY_PATH
export LD_LIBRARY_PATH=/<exago_install>/lib:$LD_LIBRARY_PATH
```
A sample opflow script `test.py` is provided, with `config.prefix()` providing the path to the installation directory of exago. For more example usage and for the tests that cover this code, see `tests/interfaces/python`:
