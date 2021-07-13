import sys
import os
import inspect
from collections import defaultdict

# Maps filenames to test cases. Each test case will register itself in the entry
# corresponding to it's filename so test cases are grouped by file.
_exago_tests = defaultdict(list)


def exago_test(testcase):
    '''Register a function as a test case'''

    def wrapper():
        name = testcase.__name__
        try:
            testcase()
        except Exception as e:
            print(f'-- Testcase {name}: FAIL')
            print('Got error: ', str(e))
            return 1
        print(f'-- Testcase {name}: PASS')
        return 0

    # Register test case with exago test for `exago_run_all_tests`
    _exago_tests[inspect.getfile(testcase)].append(wrapper)

    return wrapper

def exago_run_all_tests(file_name):
    '''Run all test cases that were registered in the caller's file.
    Invoke with `exago_run_all_tests(__file__)`'''

    fail = 0
    total = 0
    for test in _exago_tests[file_name]:
        fail += test()
        total += 1

    if fail == 0:
        print(f'All {total} test cases passed.')
    else:
        print(f'{fail} / {total} test cases failed.')

    sys.exit(fail)
