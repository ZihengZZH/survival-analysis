import numpy as np
import pandas as pd
import sys
import statistics


data_path = './dataset/'
data_path_process= './dataset_proc/'


# load TCGA-BRCA clinical data
def load_data_clinical(proc=True, test=False):
    # para proc: whether or not to load processed clinical data
    # para test: whether or not to print key aspects
    if proc:
        data_clinical = pd.read_csv(data_path_process+'20160128-BRCA-Clinical-processed.txt', sep='\t', index_col=0)
    else:
        data_clinical = pd.read_csv(data_path+'20160128-BRCA-Clinical.txt', sep="\t", header=0, index_col=0).T

    if test:
        # samples
        no_patient, no_feature = data_clinical.shape[0], data_clinical.shape[1]
        # age data
        age = data_clinical['years_to_birth'].dropna().tolist()
        age = [int(x) for x in age]
        age_median, age_min, age_max = statistics.median(age), min(age), max(age)
        # vital status data
        vital = data_clinical['vital_status'].tolist()
        (living, deceased) = (vital.count(0), vital.count(1)) if proc else (vital.count('0'), vital.count('1'))
        # gender data
        gender = data_clinical['gender'].tolist()
        male, female = gender.count('male'), gender.count('female')
        # print key aspects 
        print("#features: %d\t#patients: %d" % (no_feature, no_patient))
        print("#living: %d\t#deceased: %d" %(living, deceased))
        print("median age: %d\tage range: %d-%d" %(age_median, age_min, age_max))
        print("#male: %d\t#female: %d" %(male, female))

    return data_clinical


# load TCGA-BRCA RNASeqGene data
def load_data_RNASeq(proc=True, label=True, lymph=False, test=False, raw_count=False):
    # para proc: whether or not to load processed RNASeq data
    # para label: whether or not to load RNASeq data with labels
    # para test: whether or not to print key aspects
    if proc and label:
        data_RNASeq = pd.read_csv(data_path_process+'20160128-BRCA-RNAseqGene-label.txt', sep='\t', index_col=0)
    elif proc and lymph:
        data_RNASeq = pd.read_csv(data_path_process+'20160128-BRCA-RNAseqGene-lymph.txt', sep='\t', index_col=0)
    elif proc and not label:
        data_RNASeq = pd.read_csv(data_path_process+'20160128-BRCA-RNAseqGene-processed.txt', sep='\t', index_col=0)
    elif raw_count:
        data_RNASeq = pd.read_csv(data_path_process+'20160128-BRCA-RNAseqGene-raw-counts-label.txt', sep='\t', index_col=0)
    else:
        data_RNASeq = pd.read_csv(data_path+'20160128-BRCA-RNAseqGene.txt', sep='\t', header=0, index_col=0).T
        # RPKM index list
        RPKM_index_lst = [x*3+2 for x in list(range(878))] if not raw_count else [x*3 for x in list(range(878))]
        # RPKM partition of the RNASeq
        data_RNASeq_RPKM = data_RNASeq.iloc[RPKM_index_lst]
        # save processed data 
        save_processed_data(data_RNASeq_RPKM, raw_count, data_type='RNASeq')
    
    if test:
        # features (RNA gene)
        no_sample, no_feature = data_RNASeq.shape[0], data_RNASeq.shape[1]
        # print key aspects
        print("#features: %d\t#samples: %d" % (no_feature, no_sample))

    return data_RNASeq


# save processed data to file (suitable for various data types)
def save_processed_data(data, raw_count, data_type=None):
    # para data: the processed data to store
    # para raw_count: whether to load RPKM or raw counts 
    if not data_type:
        print("\nNo data types identified")
        return
    if raw_count and data_type == 'raw_counts':
        data.to_csv(data_path_process+'20160128-BRCA-RNAseqGene-raw-counts-label.txt', sep='\t')
        print("\nProcessed RNASeq data (w/ raw counts & labels) has been successfully written to file ...")
    if raw_count:
        data.to_csv(data_path_process+'20160128-BRCA-RNAseqGene-raw-counts.txt', sep='\t')
        print("\nProcessed RNASeq data (w/ raw counts) has been successfully written to file ...")
    if data_type == 'clinical':
        data.to_csv(data_path_process+'20160128-BRCA-Clinical-processed.txt', sep='\t')
        print("\nProcessed Clinical data has been successfully written to file ...")
    if data_type == 'RNASeq':
        data.to_csv(data_path_process+'20160128-BRCA-RNAseqGene-processed.txt', sep='\t')
        print("\nProcessed RNASeq data has been successfully written to file ...")
    if data_type == 'RNASeq_label':
        data.to_csv(data_path_process+'20160128-BRCA-RNAseqGene-label.txt', sep='\t')
        print("\nProcessed RNASeq with labels data has been successfully written to file ...")
    if data_type == 'RNASeq_lymph':
        data.to_csv(data_path_process+'20160128-BRCA-RNAseqGene-lymph.txt', sep='\t')
        print("\nProcessed RNASeq with #lymph nodes data has been successfully written to file ...")


'''ONLY EXECUTE ONCE'''
# label the RNASeq data with the clinical data
def label_RNASeq_data(lymph=False):
    # para lymph: whether to label RNASeq with lymph data
    if not lymph:
        print("\nLabel the RNASeq data with the vital_status in clinical data ...")
    else:
        print("\nLabel the RNASeq data with the #lymph nodes in clinical data ...")
    # load data first
    data_clinical = load_data_clinical()
    data_RNASeq = load_data_RNASeq(label=False)
    # find the RNASeq dataframe index
    index_RNASeq = data_RNASeq.index.tolist()
    index_RNASeq = [x[:12].lower() for x in index_RNASeq]
    # drop the samples in clinical copy dataframe if there is no corresponding RNASeq data
    data_clinical_copy = data_clinical.copy()
    for index_clinical, _ in data_clinical_copy.iterrows():
        if index_clinical not in index_RNASeq:
            data_clinical_copy = data_clinical_copy.drop([index_clinical])
    
    # find the clinical dataframe index
    data_clinical_copy = data_clinical_copy.sort_index()
    index_Clinical = data_clinical_copy.index.tolist()
    index_Clinical = [index.upper() for index in index_Clinical]
    # drop the samples in RNASeq copy dataframe if there is no corresponding clinical data
    data_RNASeq_copy = data_RNASeq.copy()
    data_RNASeq_copy['label'] = pd.Series(np.zeros(data_RNASeq_copy.shape[0]), index=data_RNASeq_copy.index)
    for index_RNASeq, _ in data_RNASeq_copy.iterrows():
        if index_RNASeq[:12] not in index_Clinical:
            data_RNASeq_copy = data_RNASeq_copy.drop([index_RNASeq])
        # otherwise, label the samples 
        else:
            if not lymph:
                label = data_clinical_copy.loc[index_RNASeq[:12].lower(),:]['vital_status']
            else:
                label = data_clinical_copy.loc[index_RNASeq[:12].lower(),:]['number_of_lymph_nodes']
                label = 0 if pd.isnull(label) else label
                label = 1 if label > 0 else label
            data_RNASeq_copy.loc[index_RNASeq, 'label'] = label
    
    # save the processed data for further analysis
    save_processed_data(data_clinical_copy, False, data_type='clinical')
    if not lymph:
        save_processed_data(data_RNASeq_copy, False, data_type='raw_counts')
    else:
        save_processed_data(data_RNASeq_copy, False, data_type='RNASeq_lymph')


'''
description of clinical data and RNASeqGene data (after pre-processing):

CLINICAL DATA
#features       #samples    #living     #deceased
18              779         658         121  

RNASEQGENGE DATA
#features       #samples
20533           878 
'''
