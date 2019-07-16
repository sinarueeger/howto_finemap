
## replicate figure7 in the FINEMAP manuscript
## Fine-mapping of 15q21/LIPC region associated with high-density lipo- protein cholesterol.



## FINEMAP: Benner et al. 2016
## https://academic.oup.com/bioinformatics/article/32/10/1493/1743040

## HDL summary stats:
## Data: https://www.ncbi.nlm.nih.gov/pubmed/25961943
## Download: http://diagram-consortium.org/2015_ENGAGE_1KG/
## Chromosome and position (build 37, base-pairs)

## libraries needed
## vroom
## ggGWAS
##
library(ggplot2)
library(dplyr)
theme_set(theme_bw())
library(ggGWAS) # for plotting
library(readr) ## read_delim
library(data.table) ## reading in too

fs::dir_create(here::here("data"))
DIR_DATA <- here::here("data")

## Data ------------------------
## -----------------------------
## Chromosome and position (build 37, base-pairs)
url_summary_stats <- "mccarthy.well.ox.ac.uk/publications/2015/ENGAGE_1KG/HDL_Meta_ENGAGE_1000G.txt.gz"

download.file(url_summary_stats, glue::glue("{DIR_DATA}/hdl_summarystats.txt.gz"))
dat_raw <- vroom::vroom(glue::glue("{DIR_DATA}/hdl_summarystats.txt.gz"),
                        col_select = list(p = `p-value`, n = n_samples, everything()))

## calc p value again
dat_raw$p_sina <- GWAS.utils::z2p(abs(dat_raw$beta/dat_raw$se))
qplot(-log10(p), -log10(p_sina), data = dat_raw %>% sample_n(1000)) +
  geom_abline(intercept = 0, slope = 1)

## check how many p values are 0
dat_raw %>% filter(p_sina == 0)
dat_raw %>% filter(p == 0)

## lets display this with a qqplot, just a check
ggplot(data = dat_raw) +
  stat_gwas_qq_hex(aes(y = p_sina)) +
  geom_abline(intercept = 0, slope = 1)

## we are only interested in the region around
CHR <- 15
BP_SNP <- 58680954 ## centered around https://www.ncbi.nlm.nih.gov/snp/rs2043085
WINDOW <- 1e6
BP_FROM <- BP_SNP - WINDOW/2
BP_TO <- BP_SNP + WINDOW/2

dat <- dat_raw %>%
  filter(between(pos, BP_FROM, BP_TO) & chr == glue::glue("chr{CHR}"))

## lets make a locus plot
ggplot(data = dat) +
  geom_point(aes(pos, -log10(p_sina)))

## prepare tools ----------------------
## --------------------------------------

## path to tools
## download from here:
## www.christianbenner.com
LDSTORE <- here::here("bin", "ldstore_v1.1_MacOSX", "ldstore")
FINEMAP <- here::here("bin", "finemap_v1.3.1_MacOSX", "finemap_v1.3.1_MacOSX")
PLINK2 <- here::here("bin", "plink2")

system_glued <- function(x)
{
  system(glue::glue(x))
}


## prepare reference data ----------------
## ---------------------------------------

## sample info
## get chr 15
## on GRCh37
# samples
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v2.20130502.ALL.ped
# ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.xlsx
url_reference_panel <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

download.file(url_reference_panel,
              glue::glue("{DIR_DATA}/1000genomes_chr15.vcf.gz"))
download.file("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi",
              glue::glue("{DIR_DATA}/1000genomes_chr15.vcf.gz.tbi"))
download.file("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v2.20130502.ALL.ped",
              glue::glue("{DIR_DATA}/samples_info.ped"))

## include only certain samples
## include only certain snps

## vcf to plink
FILE_1KG <- glue::glue("{DIR_DATA}/1000genomes_chr15")

## TODO: select individuals
##
## --keep accepts a space/tab-delimited text file with family IDs in the first column and within-family IDs in the second column, and removes all unlisted samples from the current analysis. --remove does the same for all listed samples.
keep <- data.table::fread( glue::glue("{DIR_DATA}/samples_info.ped")) %>%
  filter(Population %in% c("GBR", "FIN")) %>%
  select(`Family ID`, `Individual ID`) %>%
  mutate(`Family ID` = 0)

readr::write_delim(keep, path = glue::glue("{DIR_DATA}/keep.txt"), delim = " ", col_names = FALSE)


system_glued("{PLINK2} --vcf {FILE_1KG}.vcf.gz --keep {DIR_DATA}/keep.txt --max-alleles 2 --mind 0.01 --geno 0.01 --hwe 1e-6 --make-bed --out {FILE_1KG} --from-bp {BP_FROM} --to-bp {BP_TO} --chr {CHR}")


## prepare master file -------------------
## ---------------------------------------
## has to look like this
#z;ld;snp;config;cred;log;n_samples
#dataset1.z;dataset1.ld;dataset1.snp;dataset1.config;dataset1.cred;dataset1.log;5363
#dataset2.z;dataset2.ld;dataset2.snp;dataset2.config;dataset2.cred;dataset2.log;5363

FILE <- glue::glue("{DIR_DATA}/data-hdl")
sample_size <- median(dat_raw$n) ## look that up in the paper
master <- data.frame(z = glue::glue("{FILE}.z"),
                     ld = glue::glue("{FILE}.ld"),
                     snp = glue::glue("{FILE}.snp"),
                     config = glue::glue("{FILE}.config"),
                     cred = glue::glue("{FILE}.cred"),
                     log = glue::glue("{FILE}.log"),
                     n_samples = sample_size
                     )

readr::write_delim(master, path = "master", delim = ";")



## prepare z data ------------------------
## ---------------------------------------

# create data.z
#rsid chromosome position allele1 allele2 maf beta se
#rs1         10 1 T C 0.35 0.0050 0.0208
#rs2         10 1 A G 0.04 0.0368 0.0761
#rs3         10 1 G A 0.18 0.0228 0.0199
## check if they are in the bim file
snps_1kg <- data.table::fread(glue::glue("{FILE_1KG}.bim"), header = FALSE) %>%
  rename(chr = V1, rsid = V2, pos = V4, allele2 = V5, allele1 = V6)
# Allele 1 (corresponding to clear bits in .bed; usually minor) > other allele
# Allele 2 (corresponding to set bits in .bed; usually major) > reference allele

## then join + rearrange
data_z <- dat %>%
  rename(chromosome = chr, allele1 = reference_allele, allele2 = other_allele) %>%
  inner_join(snps_1kg %>% select(pos, rsid, allele1, allele2)) %>%
  rename(position = pos) %>%
  mutate(maf = 0.1, chr = CHR) %>% ## overwrite chromosome
  select(rsid, chromosome, position, allele1, allele2, maf, beta, se)

## TODO: get proper MAF
data_z <- data_z %>% slice(1:100)

readr::write_delim(
  data_z,
  path = glue::glue("{FILE}.z"),
  delim = " "
)



## --incl-variants 		Extract LD information for variants given in the specified text file.
## The specified file has 5 columns with a header: RSID, position, chromosome, A_allele and B_allele 		Requires --matrix or --table
incl_variants <- data_z %>%
    select(rsid, position, chromosome, allele1, allele2) %>%
    rename(RSID = rsid, A_allele = allele1, B_allele = allele2)

readr::write_delim(incl_variants,
              path = glue::glue("{FILE}-incl-variants"),
              delim = " ")

## generate LD data -----------------------
## ----------------------------------------

if(FALSE) {
  ## sub bcor file
  system_glued(
    "{LDSTORE} --bplink {FILE_1KG} --bcor {FILE_1KG}.bcor --n-threads 1 --accuracy low --ld-thold 0 --ld-n-samples-avail-prop 0"
  )
  #
  ##  --incl-range {BP_FROM}-{BP_TO}


}

## create LD matrix
system_glued(
  "{LDSTORE} --bcor {FILE_1KG}.bcor_1 --matrix {FILE}.ld --incl-variants {FILE}-incl-variants"
)


out <- readr::read_delim("out", " ")
head(out %>% filter(RSID %in% data_z$rsid))

## run FINEMAP ----------------------------
## ----------------------------------------

system_glued(
  "{FINEMAP} --sss --in-files master --dataset 1 --log"
)



## read in results ------------------------
## ----------------------------------------

config <- data.table::fread(glue::glue("{FILE}.config"))
snp <- data.table::fread(glue::glue("{FILE}.snp"))

ggplot(data = snp) +
  ## add log10bf
  geom_point(aes(position, abs(log10bf), color = corr_group)) +
  scale_color_distiller("LD", type = "div", palette = "Spectral", limits = c(0, 1)) +

  ## color the causal variants
  geom_point(aes(position, abs(log10bf)), shape = 2, data = snp %>% filter(group == 1)) +
    labs(
      title = "Locuszoom plot for HDL GWAS",
      subtitle = glue::glue("Summary statistics for chromosome {CHR}, {BP_FROM}-{BP_TO}"),
      caption = glue::glue("Data source: Surakka et al. 2015: {url_summary_stats}\n 1000 genomes ref panel: {url_reference_panel}")
    )


#
# log10bf
#
# prob column contains the marginal Posterior Inclusion Probabilities (PIP). The PIP for the l th SNP is the posterior probability that this SNP is causal.
#
# log10bf column contains the log10 Bayes factors. The Bayes factor quantifies the evidence that the l th SNP is causal with log10 Bayes factors greater than 2 reporting considerable evidence
#
# group column contains the group number that the SNP belongs to
#
# corr_group column contains the correlation with the marginally most significant SNP among SNPs in the same group with this SNP
#
# prob_group column contains the posterior probability that there is at least one causal signal among SNPs in the same group with this SNP
#
# log10bf_group column contains the log10 Bayes factors for quantifying the evidence that there is at least one causal signal among SNPs in the same group with the SNP. Bayes factors greater than 2 report considerable evidence
#
#
