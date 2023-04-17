##' 'AlienindexR' calculates Alien Index (AI) for genes in target species from BLAST results between this species and donor/recipient species. 
##'
##' The output by BLAST or diamond BLAST will be used as input of 'AlienindexR' and should be in a tabular format with 12 columns. 
##'
##' For BLAST, the output parameter should be: -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'. 
##'
##' For diamond BLAST, the output parameter should be: -outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore. 
##'
##'
##' @title AlienindexR
##' @param blast_output_donor Input file with 12 columns in tab-delimited format
##' @param blast_output_recipient Input file with 12 columns in tab-delimited format
##' @return dataframe
##' @importFrom stats complete.cases
##' @export 
##' @author Xiao Chen
##' @references 1. X Chen, L Yang. Manuscript in preparation.
##' @examples
##' test_donor <- read.table(gzfile(system.file("extdata", c("test_0.tab", package = "AlienindexR")), header=TRUE, sep="\t")
##' test_recipient <- read.table(gzfile(system.file("extdata", c("test_1.tab", package = "AlienindexR")), header=TRUE, sep="\t")
##' AlienindexR(test_donor, test_recipient)

## Main
AlienindexR <- function(blast_output_donor, 
                   blast_output_recipient
                   ) {

    donor <- blast_output_donor
    recipient <- blast_output_recipient
    colnames(donor) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
    colnames(recipient) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
    # shared genes between donors and recipients
    shared.genes = unique(intersect(donor$qseqid, recipient$qseqid))

    # calculate AI for each shared homolog between donors and recipients
    AI.df <- data.frame()
    if (length(shared.genes) > 0) {
        for (gene in shared.genes) {
			donor.sub = donor[donor$qseqid == gene,,drop=F]
			recipient.sub = recipient[recipient$qseqid == gene,,drop=F]
			best_recipient_evalue = min(recipient.sub$evalue)
			best_donor_evalue = min(donor.sub$evalue)
			AI = log(best_recipient_evalue + 1e-200) - log(best_donor_evalue + 1e-200)
			AI.df.sub = as.data.frame(list(AI))
			colnames(AI.df.sub) = "Alien.Index"
			rownames(AI.df.sub) = gene
			AI.df = rbind(AI.df, AI.df.sub)
        }
    } else {
        message("No shared homologs between donors and recipients!")
    }
    
    return(AI.df)
}
