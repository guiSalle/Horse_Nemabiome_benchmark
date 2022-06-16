# Horse_Nemabiome_benchmark

This repository contains codes and input data used to benchmark the horse nemabiome data using the internal transcribed spacer 2 region, and to test the predictive ability of an amplicon located within the COI mitochondrial gene. 


# Available scripts
These are within the ```scripts_for_sub``` directory.
* ```dada2_ITS_mock_benchmark.sh``` : shell script to process the raw sequencing data and produce the ASV tables using various DADA2 parameters.

* ```dada2_nmb_paired.sh``` : shell script to process the raw sequencing data and produce the ASV tables using various DADA2 parameters. This script relies on the run_dada_paired.R script that is freely available here: https://github.com/qiime2/q2-dada2/blob/master/q2_dada2/assets/run_dada_paired.R 

* ```minimapCOI_benchmark.sh```: shell script to produce genome coverage statistics over the COI amplicon with the minimap software and using various set of parameters.

* ```barcode.ipynb```: python script to compute similarity across sequences and output a COI consensus sequence

* ```NMB_COMP_sub160622_git.R```: R script for pre-processed data analysis

# Data directories

*```Cya_Str_cox1_BOLD_GB_mito.fasta```: COI database extraced from BOLD, GENBANK and mitogenome sequences

*```mock_DADA2```: ASV table counts yielded by the dada2_ITS_mock_benchmark.sh data

*```life_stage```: Contain two folders corresponding to the data processing of amplicon sequencing data applied to different life-stages and processed with the best set of parameters.


<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
