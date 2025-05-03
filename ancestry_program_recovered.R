start_time <- Sys.time()
# PART 1
if (!require("BiocManager"), quietly = TRUE))
 install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("AnnotationHub")
install.packages("pacman")
pacman::p_load(AnnotationHub, tidyr, dplyr, purrr, readr, readxl)

# quickly read notepad WITHOUT TOP TEXT into R
df <- read_tsv('C:\\Users\\22vyaish\\Documents\\HPF\\RLMCC_AncestryDNA_ydna.txt')

# find data only on Y chromosome, or 'chromosome 24'
# AND gets rid of alleles where value is 0

df <- subset(df, df$chromosome == 24)
df <- subset(df, df$allele1!=0)

# PART 2a: USING MY OWN DB
# I want to import my excel as another dataframe,
dbSNP <- read_excel('C:\\Users\\22vyaish\\Documents\\HPF\\dbSNP_IV.xlsx', sheet = "Final_Data")

# join two databases into one with only valid RSIDs
matched <- inner_join(df, dbSNP, by = c("rsid" = "rsid"))

# is there a mutation?
tf < c(TRUE)
matched$tf <- tf
for (i in 1:length(matched$allele1)) {
    if (matched$AltAllele[i] == matched$allele1[i]) {
        matched$tf[i] <- TRUE
    } else {
        matched$tf[i] <- FALSE
    }
    if (matched$AltAllele[i] == matched$allele2[i]) {
        matched$tf[i] <- TRUE
    }
}
for (l in 1:1) {
if (nrow(df2) > 0) {break}
# PART 2b: USING ONLINE DATA
# creates an empty datafrane and names the colunns
df2 <- data.frame(matrix(ncol = 4, nrow = nrow(df)))
names(df2)[1] <- "RSID"
names(df2)[2] <- "Reference_Allele"
names(df2)[3] <- "Alternate_Allele"
names(df2)[4] <- "Origin"

# read code of URL into vector readpage
my_url <- paste0("https://www.ncbi.nlm.nih.gov/snp/", df$rsid)
mypattern = '\"([^<]*)\"'

# run thru for all
k <- 1
for (k in 1:length(my_url)) {

    # read lines into test variable
    test <- readLines(my_url[k])

    # extrapolate 10 lines of test with the alleles
    datalines = grep(mypattern, test[120:130], value=TRUE)

    # credit: https://statistics.berkeley.edu/computing/fags/reading-web-pages-r
    getexpr = function(s,g)substring(s,g,g+attr(g,'match.length')-1)
    gg = gregexpr(mypattern,datalines)
    matches = mapply(getexpr,datalines,gg)
    result = gsub(mypattern, '\\1', matches)
    names(result) = NULL

    # 3rd part is always the REF>ALT
    a_var <- result[3]

    ## allele assignments
    ref_allele <- substr(a_var, 24, 24)
    alt_allele <- substr(a_var, 26, 26)

    df2[k,] <- c(df$rsid[k], ref_allele, alt_allele, NA)
}}
# 1st: subset out the 'bad' data from df2

nucleotide_vars <- c('A', 'C', 'T', 'G')
df2 <- subset(df2, df2$Reference_Allele %in% nucleotide_vars)

# 2nd: subset to find matches between given data and df2
smatch2 <- inner_join(df, df2, by = c("rsid" = "RSID"))

# 3rd: then use RSid results to find potential origins
snp_index_online <- read_excel('C:\\Users\\22vyaish\\Documents\\HPF\\SNP_Index.xlsx', sheet = "Human")
snp_index_online <- snp_index_online %>% drop_na(rs_numbers, Alternate_Names)
phylo_tree <- inner_join(smatch2, snp_index_online, by = c("rsid" = "rs_numbers"))

tf2 <- c(TRUE)
smatch2$tf2 <- tf2
i <- 1
for (i in 1:nrow(smatch2)) {
    if (smatch2$Alternate_Allele[i] == smatch2$allele1[i]) {
        smatch2$tf2[i] <- TRUE
    } else {smatch2$tf2[i] <- FALSE}
    if (smatch2$Alternate_Allele[i] == smatch2$allele2[i]) {
        smatch2$tf2[i] <- TRUE
    }
}
final_df <- subset(smatch2, smatch2$tf2 == TRUE)
matched <- subset(matched, matched$tf == TRUE)
sink("Result.txt")
cat("A note about YDNA assignments:")
cat("\n")
cat("YDNA analysis can only be done on men, since only men have a Y chromosome.
However, since YDNA is only inherited from a genetic father, women can use
YDNA test results from a genetic brother or father. For more information
on Y-DNA analysis, check out this link: https://isogg.org/wiki/Y_chromosome_DNA_tests")
cat("\n")
cat("\n")
cat("Your results:")
cat("\n")
j <- 1
for (j in 1:length(matched$Assignment)) {
    outl <- paste("I can tell you have ancestors from ", matched$Assignment[j], "with ", matched$confidence[j], "confidence. \n")
    cat(outl)
}
cat("\n")
cat("Explore your SNPs using NCBI's database: https://www.ncbi.nlm.nih.gov/snp/
Just search up the RSID, and click on the first result.")
cat("\n")
cat("Your full list of SNPs:")
cat("\n")
cat(final_df$rsid)
end_time <- Sys.time()
cat("\n")
print(paste("Program run length:", end_time - start_time, "seconds."))
sink()
file.show("Result.txt")