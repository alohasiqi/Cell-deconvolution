# Estimating the effect of tissue- and blood-derived cell reference matrices on deconvolving bulk transcriptomic datasets 

Siqi Sun1, Shweta Yadav1, Mulini Pingili1, Dan Chang1, Jing Wang1*

1 Genomics Research Center, AbbVie, 200 Sidney Street, Cambridge, MA 02139

This the repo for the code in the preprint https://scholar.google.com/citations?view_op=view_citation&hl=en&user=sjt8cbcAAAAJ&citation_for_view=sjt8cbcAAAAJ:Se3iqnhoufwC

## Analyses workflow


## Key Points
### Tissue-derived CRMs showed higher goodness-of-fit compared to blood-derived CRMs for deconvolving bulk tissue transcriptomics.
### All CRMs yield consistent goodness-of-fit for deconvolving bulk blood transcriptomics.
### Tissue-derived CRMs represent more accurate cellular proportion estimates and reveal more treatment-related cell types, and the specific tissue type is relevant to the deconvolution performance.

## Abstract
### Background 
Cell deconvolution is a method used to characterize the composition of the mixed cell population in bulk transcriptomic datasets. Tissue- and blood-derived cell reference matrices (CRMs) are typically used, but their impact on deconvolution has yet to be evaluated. 

### Methods
Recursive feature elimination and random forest methods were applied to build tissue- and blood-derived CRMs using single-cell RNA sequencing (scRNA-seq) datasets of inflammatory bowel disease (IBD) which comprises two subtypes (Crohn’s disease (CD) and ulcerative colitis (UC)). Combining with three published blood-derived CRMs (IRIS, LM22, and ImmunoStates), public bulk transcriptomes datasets and simulated datasets generated from scRNA-Seq were used to evaluate the deconvolution performance by goodness-of-fit and cell fractions correlation. Additionally, two infliximab-treated bulk datasets were used to compare the treatment-related cell types revealed by tissue- and blood-derived CRMs. Lung adenocarcinoma (LUAD) single-cell and TCGA bulk transcriptomic datasets were also used for evaluation.

### Results
For CD and UC tissue bulk datasets, tissue-derived CRMs showed better deconvolution goodness-of-fit scores compared to blood-derived CRMs. Meanwhile, tissue-derived CRMs represented more accurate cellular proportion estimates for most cell types, such as immune and stromal cells in tissue pseudobulk datasets. They also revealed more treatment-related cell types. Additionally, CRMs derived from the scRNA-seq datasets from CD and UC colon samples had similar performance. In contrast, CRMs derived from CD ileum scRNA-seq dataset had better performance in the ileum bulk datasets. All CRMs yielded consistent deconvolution performance when deconvolving blood bulk transcriptomics. The similar results have also been shown using LUAD datasets.  

### Conclusions
Our results emphasize the importance of selecting appropriate CRMs for cell deconvolution, particularly in bulk tissue transcriptomes in immunology and oncology. Such considerations can be extended to encompass other disease implications.
