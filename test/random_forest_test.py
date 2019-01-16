import unittest
import numpy as np
import src.random_forest as RF

'''
unittest: The Python unit testing framework
unittest supports test automation, sharing of setup and shutdown code for tests, aggregation of tests into collections, and independence of the tests from the reporting framework.
to run the unittest, python -m unittest test.random_forest_test
'''


class RandomForestTest(unittest.TestCase):
    # method called to prepare the test fixture
    def setUp(self):
        pass
    
    # method called immediately after the test method has been called and the result recorded
    def tearDown(self):
        pass
    
    def test_run_random_forest(self):
        RF.run_random_forest()


if __name__ == "__main__":
    unittest.main()