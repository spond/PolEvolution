##POL EVOLUTION

This collection of HyPhy (http://www.hyphy.org) scripts perform sequence analysis
on collections of longitudinal HIV-1 pol (bulk) sequence isolates, sampled from
multiple individuals at multiple time points. More details of the analysis can 
be found in this paper: http://www.ncbi.nlm.nih.gov/pubmed/23840830

This library can 

* Estimate a nucleotide substitution rate for individuals 
* Test whether or not the molecular clock is violated for a particular individuals
* Compare global nucleotide substitution rates between groups of patients
* Estimate dN/dS values for individuals, and groups of individuals
* Use CTL haplotype information to partition the sites for each individual into epitope and non-epitope sets (not the same for all individuals), and compare rates of evolution (and dN/dS) between them in an individual, within and between groups of individiuals
* Use ARV information to perform the same types of analyses for DRAM vs not DRAM sites

##USAGE


Assuming, you have a command line version of HyPhy installed (GUI versions are 
OK as well), run the following command:

	HYPHYMP CPU=1 /path/to/PolEvolution/HBL/RateAnalysis.bf

And enter the following values at the prompts (use HLA information, read HLA information 
from data/hla.csv in the PolEvolution Directory, do Codon analyses, TN93 model, sequence data from data/clean_pol.fas,
split subjects into groups based on Dual Infection status read from data/group.csv, and save the results to results/test.csv)


<pre>


			+------------------------------+
			|Include HLA or DRM information|
			+------------------------------+


	(1):[Neither] No HLA or DRM information
	(2):[HLA] Provide a .csv file with HLA haplotypes for each patient ID
	(3):[DRAM] Provide a .csv file with ARV lists to determine relevant DRM for each patient ID

 Please choose an option (or press q to cancel selection):2

/Volumes/TeraMonkey/Users/sergei/Coding/PolEvolution/CTL/The CTL haplotype file in .CSV format:../data/hla.csv


			+----------------------+
			|Perform codon analyses|
			+----------------------+


	(1):[No] Only perform nucleotide rate analyses [FASTER]
	(2):[Yes] Perform codon-based rate analyses (dN/dS) in addition to nucleotide rate analyses

 Please choose an option (or press q to cancel selection):2


			+-------------------------+
			|Use this nucleotide model|
			+-------------------------+


	(1):[JC69] Jukes Cantor 69 (all equal rates, all equal frequencies)
	(2):[F81] Felsenstein 81 (all equal rates, empirical frequencies)
	(3):[TN93] Tamura and Nei 93 (two transition rates, shared transversion rate, empirical frequencies)
	(4):[GTR] General Time Reversible (6 substitution rates, empirical frequencies)

 Please choose an option (or press q to cancel selection):3

[path]/PolEvolution/HBL/Load the combined HIV polymerase sequence file:../data/clean_pol.fas


			+----------------------------------------------------------------------------------------------------+
			|Split subjects into two groups based on external information (e.g. dual infection, treatment status)|
			+----------------------------------------------------------------------------------------------------+


	(1):[No] Perform a joint analysis on all subjects
	(2):[Yes] Split the subjects into two groups, perform two-group analyses, and compare the rates between groups

 Please choose an option (or press q to cancel selection):2

[path]/PolEvolution/HBL/Partitioning list of the form: pid, [0/1]::../data/group.csv

Descriptive label for 'positive' patients, e.g. Treated:Dually Infected

[path]/PolEvolution/HBL/Save the resulting .CSV file to::../results/test.csv


</pre>



##DATA
	
	
