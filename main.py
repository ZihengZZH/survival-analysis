import numpy as np
import pandas as pd
import src.utility as ut
import src.random_forest as rf
import src.gradient_boost as gb
import src.survival_analysis as sa


if __name__ == "__main__":
    # rf.tune_hyperparameters()
    # gb.tune_hyperparameters()
    # ut.label_RNASeq_data(lymph=True)
    rf.run_random_forest(load=True, model_no=3)
    gb.run_gradient_boost(load=True, model_no=5)
    sa.survival_analysis_with_all_RNASeq('rf')
    sa.survival_analysis_with_all_RNASeq('gbrt')
    sa.draw_log_p_values()
