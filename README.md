
[![DOI](https://zenodo.org/badge/754718392.svg)](https://zenodo.org/doi/10.5281/zenodo.10636308)

![Wellcome Centre for Anti-Infectives Research Logo](static/wcar.png)

![Coverage track for gene Tb927.10.5620](make_images_code/bw_file/full_images/Tb927.10.5620_coverage.png)


# Post-transcriptional reprogramming by thousands of mRNA untranslated regions in trypanosomes
we describe the UTR sequences that control gene expression in the context of a constitutively transcribed genome
and conclude that thousands of UTRs post-transcriptionally reprogram gene expression profiles in trypanosomes.

## Code Repository
We deposited and divided in subsection the code used for the data analysis

-  <a href="alignment_code/">alignment_code</a> <br>
this section contains the scripts to align the fastq files to the reference genome

-  <a href="pca_code/">alignment_code</a> <br>
This section contains the scripts to visualize a PCA on the read counts for the analysed samples <br>
The scripts are implemented in the make_pca.ipynb notebook

-  <a href="peak_finding_code/">peak_finding_code</a> <br>
This section contains the scripts to identify the regions in the genome responsible for positive or negative selection.<br><br>
The notebook <b>b_sample_analysis.ipynb</b> shows the step used for the analysis of the <br>
positive fragments, from peak finding to creating the paper table (Paper_Table_Blasticidine.csv). <br>
This notebook comments on all the steps of the analysis, but for some reason, the github rendering is broken <br>
I suggest cloning the repo and visualising it on your machine. <br><br>
The notebook <b>g_sample_analysis.ipynb</b> shows the step used for the analysis of the <br>
negative fragments, from peak finding to creating the paper table (Paper_Table_Ganciclovir.csv). <br>
