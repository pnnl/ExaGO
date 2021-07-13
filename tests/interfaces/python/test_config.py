from test_utilities import *
from check_preconditions import check_preconditions
check_preconditions()


@exago_test
def test_types():
    from exago import config
    assert isinstance(config.libraries(), dict)
    assert isinstance(config.libraries()['pflow'], str)


@exago_test
def test_paths():
    from exago import config
    for k in config.libraries().keys():
        assert os.path.exists(config.libraries()[k]), \
            'Path for library' + k + 'does not exist!'


exago_run_all_tests(__file__)
