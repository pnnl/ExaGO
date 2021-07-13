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
              'The python tests require the python bindings to have already been installed.\n')
        sys.exit(1)


def check_exago_config():
    '''Check that ExaGOConfig can be imported with the current PYTHONPATH'''
    try:
        from exago import config
    except ImportError as e:
        print('Could not import exago_config!')
        sys.exit(1)


def check_datafiles():
    '''Check for datafiles needed by test suite. These should be installed with
    `make install` in the build phase.'''

    from exago import config
    for path in (
            config.prefix(),
            os.path.join(config.prefix(), 'lib'),
            os.path.join(config.prefix(), 'share', 'exago', 'datafiles'),
            os.path.join(config.prefix(), 'share', 'exago', 'options')):
        assert os.path.exists(path), \
            f'Attempted to deduce path {path} from ExaGOConfig, but path does not exist!'


def check_preconditions():
    '''Ensure all preconditions for testing python bindings have been met.'''
    check_import()
    check_exago_config()
    check_datafiles()
