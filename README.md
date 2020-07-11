# A machine learning approach to optimizing cell-free DNA sequencing panels: with an application to prostate cancer

This repository contains supporting material referenced in the publication with the above title. The manuscript can be found both on bioarxiv and in an upcoming publication. 

1. bioarxiv doi: [https://doi.org/10.1101/2020.04.30.069658](https://doi.org/10.1101/2020.04.30.069658)
2. Upcoming publication: TBD


Key Files: 

1.  [comparison_panels/](comparison_panels/)  
Contains variant location compositions for the frequency, union-existing, and orchid panels.

1.  [comparison_panels/orchid_panel.bed](comparison_panels/orchid_panel.bed)  
A bed file containing variants used to generate the hybrid capture probes for the orchid generated panel.
  
1.  [comparison_panels/supporting](comparison_panels/supporting)  
Contains panel supporting reference materal and parsing scripts.
  
1.  [notebooks/panel_generation.ipynb](notebooks/panel_generation.ipynb)  
A jupyter notebook used to generate the orchid panel, includes __Figures 1A,1B,1C, and S1__
  
1.  [notebooks/panel_evaluation.ipynb](notebooks/panel_evaluation.ipynb)  
A jupyter notebook used to assess the orchid panel, includes __Figures 2A,2B,2C,2D, and S2__
  
1.  [notebooks/reanalysis.ipynb](notebooks/reanalysis.ipynb)  
A jupyter notebook used to assess orchid panel _in silico_ performance and to model one-vs-rest feature importance, includes __Figures 1D and 3__  
  
1.  [notebooks/panel_performance.ipynb](notebooks/panel_performance.ipynb)  
A jupyter notebook to assess the number of patient variants captured by the panel, includes __Figures 4A and 4B__


Notes:

Some folders referenced in the above notebooks are not present in this repository due to their large sizes. These are made available upon request and include:

1.  __TN/__  
Contains patient tumor-normal variants

1.  __cfDNA/__  
Contains patient cfDNA detected variants

1.  __SVMv4_Enrichment_20171107__  
Contains the SVM model used to generate the orchid panel
