# MuSeq

This data resposity consist of the data files, python code and R code used to analyze Mu transposition across the E. coli genome.

The python code was last run succesfully on Dec, 12 2019 on the TACC computing cluster Lonestar5 using Python 3.7.0

Coding Files:  
MAPS.py is the python code used to identify Mu transposition sites.  
spot_detect.py is python code used to measure the distance between GFP and CFP foci.
run_classic_lasso.R

Data Files:  
Insertion counts are provided in the Mu_counts file  
Sequencing depths are provided in the Sequence_Depth file  
Links to the SRA deposit of FASTQ files are in the SRA_links.txt file
all_interaction_counts_with_norm.csv

LASSO information
run_classic_lasso.R -- R code to run LASSO analysis on processed count data from Mu transposition experiments. Designed to work directly on the data given in all_interaction_counts_with_norm.csv

all_interaction_counts_with_norm.csv -- Bin-level information on Mu transposition rates. Column headers are as follows:
index: line number (not used)
Bin_ID: bin of origin for Mu transposition
Other_bin: Destination of Mu transposition
count_raw: number of observed transposition events from Bin_ID to Other_bin
norm_val: An abundance-based normalization factor for each bin, given as the inverse of the product of the fractional abundance of Bin_ID and Other_bin in untargeted genome sequencing. As described in  STAR Methods of the accompanying paper, the normalization factor is based on the assumption that the background rate of interaction between Bin_ID and Other_bin will be proportional to the product of their abundances. Higher values of norm_val correspond to lower abundance pairings
count_normed: Product of count_raw and norm_val, used only for plotting purposes

If there are any questions please contact dwalker86@gmail.com
