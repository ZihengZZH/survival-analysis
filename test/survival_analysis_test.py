import unittest
import numpy as np
import src.survival_analysis as sa

'''
unittest: The Python unit testing framework
unittest supports test automation, sharing of setup and shutdown code for tests, aggregation of tests into collections, and independence of the tests from the reporting framework.
to run the unittest, python -m unittest test.survival_analysis_test
'''


class SurvivalAnalysisTest(unittest.TestCase):
    # method called to prepare the test fixture
    def setUp(self):
        pass
    
    # method called immediately after the test method has been called and the result recorded
    def tearDown(self):
        pass
    
    def test_survival_analysis_with_all_RNASeq(self):
        # sa.survival_analysis_with_all_RNASeq('rt')
        # sa.survival_analysis_with_all_RNASeq('gbrt')
        pass

    def test_draw_log_p_values(self):
        sa.draw_log_p_values()

if __name__ == "__main__":
    unittest.main()