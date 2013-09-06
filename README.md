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

And enter the following prompts 

	TN93dist data/test.fas data/test.dst 0.05 RESOLVE CSV

Output (all written to stderr, so can be redirected)

Example:

	Read 8 sequences of length 1320
	Will perform 28 pairwise distance calculations
	Progress:      100% (       7 links found)
	Maximum distance = 0.0890019


ARGUMENTS
---------

	
	
