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
RNA-Seq, also called whole transcriptome shotgun sequencing, uses next-generation sequencing to reveal the presence the quentity of RNA in a biological sample at a given moment. Each sample contrains raw_counts, median_length_normalized, and RPKM for every gene.

__RPKM (Reads Per Kilobase Million)__:
---
This metric attempts to normalize for sequencing depth and gene length. Read counts are normalized for differences in sequencing depth and normalized for gene size. How is how the RPKM is calculated.
1. count up the total reads in a sample and divide that number by 1,000,000 - this is our 'per million' scaling factor
2. divide the read counts by the 'per million' scaling factor, which normalizes for sequencing depth, giving you reads per million (RPM)
3. divide the RPM values by the length of the gene, in kilobases, which gives you RPKM

In this case, RPKM = {raw_counts / [sum(raw_counts) / 1000000]} / median_length_nomalized

In the project, only RPKM data are retained. 

__Label RNA Sequencing data__:
---
After labeling RNA sequencing data with vital_status in clincal data, the shape of processed Clinical data is (779, 18) and the shape of processed RNASeq data is (878, 20533+1), where 1 is the column of labels (0 or 1). Note here that the number of sample in RNASeq data is a little larger than that in Clinical data. This is because for some sample in RNASeq data, two pieces of samples belong to one patient. For instance, TCGA-A7-A13G is the code of the patient, but TCGA-A7-A13G-01A and TCGA-A7-A13G-11A co-exists. The underlying reason awaits to be explored.