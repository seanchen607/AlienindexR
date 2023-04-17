# AlienindexR 
### an R package to calculate Alien Index (AI) for genes in target species based on homolog search.

[![](https://img.shields.io/badge/release%20version-0.99.1-green.svg)](https://github.com/seanchen607/FScanR)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

‘AlienindexR’ calculates Alien Index (AI) for genes in target species from BLAST results between this species and donor/recipient species. 

## :gear: Install AlienindexR in R (>= 3.5.0)

	library(devtools)
	install_github("seanchen607/AlienindexR")
	library(AlienindexR)

## :file_folder: prepare input files

To acquire the homolog search result in the correct format, use parameter below: 

	# for BLAST
	-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
  
	# for diamond
	-outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

For details, please see [Tutorial](https://seanchen607.github.io/AlienindexR.html).

## :bulb: Run AlienindexR

To calculate Alien Index for shared homologs in the target species with donor/recipient, use codes below to run examples: 

	# Read R documentation
	?AlienindexR
	
	# load input files
	test_donor <- read.table(gzfile(system.file("extdata", "test_0.txt.gz", package = "AlienindexR")), header=FALSE, sep="\t")
	test_recipient <- read.table(gzfile(system.file("extdata", "test_1.txt.gz", package = "AlienindexR")), header=FALSE, sep="\t")
	
	# Calculate Alien Index
	AlienindexR(test_donor, test_recipient)

AlienindexR will output a two-column data frame containing the gene names and their corresponding Alien Indexes.

## :orange_book: What is Alien Index (AI)?

Alienness takes as input the result of a BLASTp search of a whole set of predicted proteins of interest against the protein libraries of donor/recipient species. 
Two individual homolog searches are performed between protein sequences in the target species (as BLAST queries) and protein libraries of donor and recipient species, respectively. 
The blast results in tab-delimited format (12 columns) contain the best blast hit for each query. 

The Alien Index (AI) is computed with the following formula:

	AI = log(best_recipient_evalue + 1e-200) - log(best_donor_evalue + 1e-200)

Where the ***best donor e-value*** is the e-value from the best hit assigned to the taxonomic donor group; 
and ***best recipient e-value*** is the e-value from the best hit assigned to the taxonomic recipient group.

For more details, please visit [Alienness](https://alienness.sophia.inrae.fr/cgi/faq.cgi).

## :pencil2: Authors

Xiao CHEN, PhD

Marine College, Shandong University, Weihai

<https://www.researchgate.net/profile/Xiao_Chen126>

[![Twitter](https://img.shields.io/twitter/url/http/shields.io.svg?style=social&logo=twitter)](https://twitter.com/intent/tweet?hashtags=AlienindexR&url=https://github.com/seanchen607/AlienindexR&screen_name=SC607)

If you use [AlienindexR](https://github.com/seanchen607/AlienindexR) in
published research, please cite the most appropriate paper(s) from this
list:

1.  **X Chen<sup>\*</sup>**, L Yang. *Manuscript in preparation*.

