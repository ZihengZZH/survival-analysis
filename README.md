# survival analysis (breast cancer)

The two datasets that were used in the project are [METABRIC](http://www.cbioportal.org/) and [TCGA](https://portal.gdc.cancer.gov/).


In this project, we mainly use the [TCGA-BRCA](https://portal.gdc.cancer.gov/projects/TCGA-BRCA) dataset and this dataset includes nine disease types:
* Adenomas and Adenocarcinomas
* Adnexal and Skin Appendage Neoplasms
* Basal Cell Neoplasms
* Complex Epithelial Neoplasms
* Cystic, Mucinous and Serous Neoplasms
* Ductal and Lobular Neoplasms
* Epithelial Neoplasms, NOS
* Fibroepithelial Neoplasms
* Squamous Cell Neoplasms


## Table of Contents
1. [TCGA-BRCA](#TCGA-BRCA)
2. [Highlights in Olson2018](#Highlights%20in%20Olson2018)
3. [Survival Analysis](#Survival%20Analysis)
4. [Comparison of classifiers' performance](#Comparison%20of%20classifiers'%20performance)



## TCGA-BRCA
---
The Cancer Genome Atlas Breast Invasive Carcinoma (TCGA-BRCA) data collection is part of a larger effort to build a research community focused on connecting cancer phenotypes to genotypes to subjects from TCGA. The specific subsets of TCGA-BRCA are applied include clinical data, RNASeqGene, methylation data, and mutation data.

### clinical data
---
The TCGA-BRCA clinical data include the followings of 1097 patients:
1. years_to_birth
2. vital_status
3. days_to_death
4. days_to_last_followup
5. tumor_tissue_site
6. pathologic_stage
7. gender
8. race
9. ethnicity 

*vital_status* (diagnosis): the state or condition of being living or deceased; also includes the case where the vital status is unknown.
*days_to_last_followup* (diagnosis): time interval from the date of last followup to the date of initial pathologic diagnosis, represented as a calculated number of days.

| Dataset | #patients | #alive | #deceased | median age* | age range* | gender (M:F) |
| ----- | ----- | ----- | ----- | ----- | ----- | ----- |
| TCGA-BRCA | 1097 | 945** | 152** | 59 | 26-90 | 12:1085 |

\* 15 patients declined to reveal their ages, so age data of 1082 patients are calculated.

\*\* calculated from the vital_status in clinical data, which conforms to the GDC portal, but contradicts to [Wenric2018](https://www.frontiersin.org/articles/10.3389/fgene.2018.00297/full). 

### RNA Sequencing
---
RNA-Seq, also called whole transcriptome shotgun sequencing, uses next-generation sequencing to reveal the presence the quantity of RNA in a biological sample at a given moment. Each sample contrains raw_counts, median_length_normalized, and RPKM for every gene.

#### RPKM (Reads Per Kilobase Million):

This metric attempts to normalize for sequencing depth and gene length. Read counts are normalized for differences in sequencing depth and normalized for gene size. How is how the RPKM is calculated.
1. count up the total reads in a sample and divide that number by 1,000,000 - this is our 'per million' scaling factor
2. divide the read counts by the 'per million' scaling factor, which normalizes for sequencing depth, giving you reads per million (RPM)
3. divide the RPM values by the length of the gene, in kilobases, which gives you RPKM

In this case, 
RPKM = {raw_counts / [sum(raw_counts) / 1000000]} / median_length_nomalized

$RPKM = \frac{raw\_counts \times 1000000}{median\_length\_nomalized \times \sum_{} raw\_counts}$

In the project, only RPKM data are retained. 

#### Label RNA Sequencing data:

After labeling RNA sequencing data with vital_status in clincal data, the shape of processed Clinical data is (779, 18) and the shape of processed RNASeq data is (878, 20533+1), where 1 is the column of labels (0 or 1). Note here that the number of sample in RNASeq data is a little larger than that in Clinical data. This is because for some sample in RNASeq data, two pieces of samples belong to one patient. For instance, TCGA-A7-A13G is the code of the patient, but TCGA-A7-A13G-01A and TCGA-A7-A13G-11A co-exists. The underlying reason awaits to be explored.



## Highlights in [Olson2018](https://psb.stanford.edu/psb-online/proceedings/psb18/olson.pdf)
---
This paper, 'Data-driven advice for applying machine learning to bioinformatics problems', analyses 13 state-of-the-art commonly-used machine learning algorithms on a set of 165 classification problems. The 13 ML algorithms come from the sklearn implementations with individual hyperparameters, and the 165 classification problem datasets come from the Penn Machine Learning Benchmark ([PMLB](https://github.com/EpistasisLab/penn-ml-benchmarks)). 

Authors concluded that both selecting the right ML algorithm and tuning its parameters are vitally importantly for most problems. Furthermore, they listed their top five recommended algorithms on 165 datasets. It is important to note that these algorithms and parameters will not work best on all supervised classification problems, and they should only be used as starting points. Here is the full list of top five recommendations:

1. GradientBoostingClassifier (51/16)
   * loss = 'deviance'
   * learning_rate = 0.1
   * n_estimators = 500
   * max_depth = 3
   * max_features = 'log2'
2. RandomForestClassifier (19/165)
   * n_estimators = 500
   * max_features = 0.25
   * criterion = 'entropy'
3. SVC (16/165)
   * C = 0.01
   * gamma = 0.1
   * kernel = 'poly'
   * degree = 3
   * coef0 = 10.0
4. ExtraTreesClassifier (12/165)
   * n_estimators = 1000
   * max_features = 'log2'
   * criterion = 'entropy'
5. LogisticRegression
   * C = 1.5
   * penalty = 'l1'
   * fit_intercept = True

## Survival Analysis 

Traditionally, survival analysis was developed to measure lifespans of individuals. The analysis can be further applied to not just traditional _births and deaths_, but any duration. Medical professionals might be interested in the time between childbirths, where a birth in this case is the event of having a child, and a death is becoming pregnant again (obviously, we are loose with our definitions of _birth and death_). 

At the time you want to make inferences about durations, it is possible that not all the death events have occured yet. The individuals in a population who have not been subject to the death event are labeled as __right-censored__, i.e., we can not view the rest of their life history due to some external circumstances. All the information we have on these individuals are their current lifetime durations, which is naturally less than their actual lifetimes. There is also __left-censorship__, where an individual's birth event is not seen.

A common mistake data analysts make a choosing to ignore the right-censored individuals. Survival analysis was originally developed to solve this type of problem, to deal with estimation when our data is right-censored. Even in the case where all events have been observed, i.e. no censorship, survival analysis is still a very useful tool to understand durations.

## Comparison of classifiers' performance 

**classification accuracies**

| RNASeq | Random Forest (RFs) | Gradient Tree Boosting (GBDT) | XGBoost |
| --- | --- | --- | --- |
| raw count | <img src="https://github.com/ZihengZZH/survival-analysis/blob/master/results/curve_random_forest.png" width=150> | <img src="https://github.com/ZihengZZH/survival-analysis/blob/master/results/curve_gradient_boost.png" width=150> | <img src="https://github.com/ZihengZZH/survival-analysis/blob/master/results/curve_xgboost.png" width=150> |
| RPKM | <img src="https://github.com/ZihengZZH/survival-analysis/blob/master/results_RPKM/curve_random_forest.png" width=150> | <img src="https://github.com/ZihengZZH/survival-analysis/blob/master/results_RPKM/curve_gradient_boost.png" width=150> | <img src="https://github.com/ZihengZZH/survival-analysis/blob/master/results_RPKM/curve_xgboost.png" width=150> |


**comparison of their performance using survival analysis**

![](https://github.com/ZihengZZH/survival-analysis/blob/master/results/p_values.png)
*p-values of every gene signature, using RNA raw count as features and tested on KM Estimate*
![](https://github.com/ZihengZZH/survival-analysis/blob/master/results_RPKM/p_values.png)
*p-values of every gene signature, using RNA RPKM as features and tested on KM Estimate*



**Statistics of p-values in above figures**

| the number of p-values < 0.01 | Random Forests | Gradient Tree Boosting | XGBoost |
| -- | -- | -- | -- |
| raw counts | 6 | 8 | 7 |
| RPKM | 2 | 12 | 9 |

