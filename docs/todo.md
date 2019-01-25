# TODO

1. WGCDA: Weighted Correlation Network Analysis -- optional
2. investigate either feature selection or not
   * #feature >> #samples (resolved by random forest)
3. ~~Random Forest implementaion on RNA-Seq~~
   * tuning hyperparameters (easier)
4. ~~Gradient Boosting Trees implementation on RNA-Seq~~
   * tuning hyperparameters (harder)
   * WHY GBT training way faster then RF
5. Support Vector Classifier implementaion on RNA-Seq
   * tuning hyperparameters 
6. gene signatures + Cox proportional hazard model
7. comparison between different algorithms
8. co-data ??
9. CorEx: Correlation Explanation -- optional


# TARGET
1. ~~classify RNA-Seq based on vital_status~~
2. access gene ranking obtained from the classification process
   1. predict lymph node metastasis with TCGA data (w/ w/o top 50 genes)
   2. survival p-value for each gene signature
      * sigRF_i = {RF_1, RF_2, ..., RF_50}
      * sigGB_i = {GB_1, GB_2, ..., GB_50}

