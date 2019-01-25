import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from smart_open import smart_open
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

from src.utility import load_data_clinical
from src.utility import load_data_RNASeq


IMPORTANT_FEATURE_RANDOM_FOREST = './results/feature_importance_rf.txt'
IMPORTANT_FEATURE_GRADIENT_BOOST = './results/feature_importance_gbrt.txt'
PLOTS_PATH_RF = './images/random_forest/'
PLOTS_PATH_GBRT = './images/gradient_boost/'
LOG_P_VALUES_PATH = './results/'


def survival_analysis_with_one_RNASeq(model_type, data, feature_list, feature_no):
    # para model_type:
    # para feature_no:
    feature_name = feature_list[feature_no][1]
    T = data[feature_name]
    E = data['label']
    T_E = pd.concat([T, E], axis=1, sort=False)
    T, E = T.tolist(), E.tolist()
    kmf_original = KaplanMeierFitter()
    kmf_original.fit(T, event_observed=E, left_censorship=False)
    # kmf_original.survival_function_
    median = kmf_original.median_

    higher_group, lower_group = [], []
    for idx in T_E.index:
        if T_E.loc[idx, feature_name] >= median:
            higher_group.append((T_E.loc[idx, feature_name], T_E.loc[idx, 'label']))
        else:
            lower_group.append((T_E.loc[idx, feature_name], T_E.loc[idx, 'label']))
    
    assert len(T_E) == len(higher_group) + len(lower_group)
    higher_group, lower_group = np.array(higher_group), np.array(lower_group)

    if len(higher_group) == 0 or len(lower_group) == 0:
        return 0
    
    kmf_higher = KaplanMeierFitter()
    ax = plt.subplot(111)
    kmf_higher.fit(higher_group[:,0], event_observed=higher_group[:,1], label='higher than median')
    ax = kmf_higher.plot(ax=ax)
    print("\nmedian survival time of higher group", kmf_higher.median_)

    kmf_lower = KaplanMeierFitter()
    kmf_lower.fit(lower_group[:,0], event_observed=lower_group[:,1], label='lower than median')
    ax = kmf_lower.plot(ax=ax)
    print("\nmedian survival time of lower group", kmf_lower.median_)

    plt.ylim(0, 1)
    plt.title("life span")

    img_name = "gene_sig_%d.png" % (feature_no+1)
    if model_type == 'rf':
        ax.get_figure().savefig(PLOTS_PATH_RF + img_name)
    elif model_type == 'gbrt':
        ax.get_figure().savefig(PLOTS_PATH_GBRT + img_name)
    ax.clear() # for the next use

    results = logrank_test(higher_group[:,0], lower_group[:,0], higher_group[:,1], lower_group[:,1], alpha=.99)
    
    return -math.log10(results.p_value)


def survival_analysis_with_all_RNASeq(model_type):
    # para model_type:
    data_RNASeq_labels = load_data_RNASeq()
    data_RNASeq_labels = data_RNASeq_labels.drop(columns=['gene'])

    feature_list = [] # list of gene signatures
    # load most important features (index and name)
    if model_type == 'rf':
        for line in smart_open(IMPORTANT_FEATURE_RANDOM_FOREST, 'r', encoding='utf-8'):
            line = line.split()
            feature_list.append((line[2], line[3]))
    elif model_type == 'gbrt':
        for line in smart_open(IMPORTANT_FEATURE_GRADIENT_BOOST, 'r', encoding='utf-8'):
            line = line.split()
            feature_list.append((line[2], line[3]))
    else:
        print("\nPlease indicate the type of model you have trained to produce important features")

    log_p_values = []
    for i in range(len(feature_list)):
        log_p_values.append(survival_analysis_with_one_RNASeq(model_type, data_RNASeq_labels, feature_list, i))

    print(log_p_values)
    save_log_p_values(model_type, log_p_values)
    

def save_log_p_values(model_type, log_p_values):
    # para model_type:
    # para log_p_values:
    filename = LOG_P_VALUES_PATH + 'log_p_values_%s.txt' % model_type
    np.savetxt(filename, log_p_values)
    print("\nlog p-values has been saved to file.")


def draw_log_p_values():
    rf_p_values = np.loadtxt(LOG_P_VALUES_PATH + 'log_p_values_rf.txt')
    gbrt_p_values = np.loadtxt(LOG_P_VALUES_PATH + 'log_p_values_gbrt.txt')
    range_p_values = list(range(1,51))
    plt.plot(range_p_values, rf_p_values, 'r', label='random forest')
    plt.plot(range_p_values, gbrt_p_values, 'b', label='gradient boost')
    plt.legend()
    plt.show()