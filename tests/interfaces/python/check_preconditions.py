import sys
import os


def check_import():
    '''Check that ExaGO python bindings can be imported with the current
    PYTHONPATH'''
    try:
        import exago
    except ImportError as e:
        print('Could not import ExaGO python bindings!\n'
              'Have you run `make install` yet?\n'
              'The python tests require the python bindings to have already '
              'been installed.\n')
        sys.exit(1)


def check_preconditions():
    '''Ensure all preconditions for testing python bindings have been met.'''
    check_import()
