import numpy as np
import pandas as pd
import src.utility as ut
import src.random_forest as rf
import src.gradient_boost as gb
import src.xgboost_classifier as xgbt
import src.survival_analysis as sa


def tune_hyperparameters_for_classifiers():
    rf.tune_hyperparameters()
    gb.tune_hyperparameters()
    xgbt.tune_hyperparameters()
    

if __name__ == "__main__":
    # run the Random Forest Classifier (pre-trained)
    rf.run_random_forest(load=True, model_no=3)
    # run the Gradient Tree Boosting Classifier (pre-trained)
    gb.run_gradient_boost(load=True, model_no=5)
    # run the XGBoost Classifier (pre-trained)
    xgbt.run_xgboost_classifier(load=True, model_no=2)

    # survival analysis on every gene signatures that were selected by different classifiers
    sa.survival_analysis_with_all_RNASeq('rf')
    sa.survival_analysis_with_all_RNASeq('gbrt')
    sa.survival_analysis_with_all_RNASeq('xgbt')

    # draw the p-values of every gene siganure
    sa.draw_log_p_values()
