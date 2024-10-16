library(tidyverse)
library(magrittr)
library(adegenet)
library(poppr)

##read in marker info

marker_info<-read_csv("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/OtsTRIR/marker_info.csv")

marker_info %<>%
  mutate(a1_count =  as.numeric(substr(a1_count, 3, nchar(a1_count)))) %>%
  mutate(a2_count =  as.numeric(substr(a2_count, 3, nchar(a2_count)))) %>%
  mutate(ind = str_remove(ind, "^\\./")) %>%
  mutate(ind = str_remove(ind, "\\.genos"))


##read in genotypes##
genos_0.1 <- read_csv("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/OtsTRIR/OtsTRIR_GTs_0.1.csv")


##clean up sample name column##

genos_0.1 %<>%
  mutate(adapter = str_extract(Sample, "[ATCG]{6}-[ATCG]{6}")) %>%
  mutate(sample_simple = str_extract(Sample, "[:upper:][:lower:]{2}[AJCU][RC]\\d{2}\\w{4}_\\d{4}")) %>%
  relocate(Sample, sample_simple, adapter)

##identify duplicate positive controls and look for how similar their genotype calls are (concordance)###

genos_0.1 %<>% 
  group_by(sample_simple) %>%
  slice_max(order_by = `On-Target Reads`, n = 2)

dups <- genos_0.1[genos_0.1$sample_simple %in% genos_0.1$sample_simple[duplicated(genos_0.1$sample_simple)],]
dups <- dups[order(dups$sample_simple),]


##Positive QC samples are the only duplicates##
##8 QC duplicates##

dups_genos <- dups[,9:ncol(dups)] 

rep_info <- matrix(ncol=ncol(dups_genos), nrow=nrow(dups_genos)/2)
colnames(rep_info) <- colnames(dups_genos)
for (j in 1:(nrow(dups_genos)/2)) {
  for (i in 1:ncol(dups_genos)) {
    rep_info[j,i] <- sum(dups_genos[(j*2)-1,i]==dups_genos[(j*2),i])
  }
}

geno_concordance <- as.data.frame(as.matrix(rep_info)) %>%
  rowMeans()
rep_data <- as.data.frame(cbind(dups[c(1:length(geno_concordance))*2,1], geno_concordance))

ggplot(data=rep_data)+geom_histogram(aes(x=geno_concordance))+theme_classic()

range(rep_data$geno_concordance)

##concordance between quality controls and corresponding samples is high (0.997 - 1.0)##

#this writes a new dataset (0.2) by choosing the samples within duplicates and keeping the one with the highest genotyping success
genos_0.2 <- genos_0.1 %>%
  group_by(sample_simple) %>%
  filter(`On-Target Reads` == max(`On-Target Reads`))

##dataset now has 86 unique samples genotyped at 352 loci and one sex marker##

##filter marker info file to remove duplicate positive controls##

marker_info<-marker_info %>%
  filter(ind %in% genos_0.2$Sample)

##look at sex genotyping success##

##for some reason sex genotypes have a comma inserted
##remove it

genos_0.2 %<>%
  mutate(`Ots_SEXY3-1` = case_when(`Ots_SEXY3-1` == "XX," ~ "XX", 
                                   `Ots_SEXY3-1` == "XY," ~ "XY", 
                                   `Ots_SEXY3-1` == "00," ~ "00"
  ))

sex_marker_info <- marker_info %>%
  filter(str_detect(marker, "SEXY")) %>%
  mutate(called_geno = replace_na(called_geno, "00")) %>%
  mutate(called_geno = case_when(called_geno == "A1HOM" ~ "XX",
                                 called_geno == "HET" ~ "XY",
                                 called_geno == "00" ~ "00"))


ggplot(data = sex_marker_info)+geom_point(aes(a2_count, a1_count, color = called_geno))+scale_colour_viridis_d(name = "Called Sex Genotype")+theme_classic()+xlab("Y-specific probe count")+ylab("Theoretical X chromosome count")+geom_abline(aes(intercept = 0, slope = 0.1))+geom_abline(aes(intercept = 0, slope = 5))+geom_abline(aes(intercept = 0, slope = 10))+geom_abline(aes(intercept = 0, slope = 0.2))+xlim(0,max(c(sex_marker_info$a1_count, sex_marker_info$a2_count)))+ylim(0,max(c(sex_marker_info$a1_count, sex_marker_info$a2_count)))+geom_abline(aes(intercept = 0, slope = 1), color = "darkred")

###Sex marker looks good. Leave as is without "correction" script##

##visualize IFI scores, individual call rates, and locus call rates##
ggplot(genos_0.2)+geom_histogram(aes(x=IFI))+geom_vline(aes(xintercept= 2.5), color="red")+theme_classic()

ggplot(genos_0.2)+geom_histogram(aes(x=`%GT`))+geom_vline(aes(xintercept= 90), color="red")+theme_classic()

##calculate locus call rates##

missingness <- (colSums(genos_0.2[,c(9:ncol(genos_0.2))] == "00" | genos_0.2[,c(8:(ncol(genos_0.2)-1))] == "0"))/nrow(genos_0.2) #warning hardcoding: "[,8:398]" is hardcoded to work on the example script using the Omy panel with 390 markers, these values will need to be changed to reflect the genotype columns of the genos r object that YOU are running. This excludes columns with metadata and genotyping results such as "sample name" "ifi" "on-target reads" etc
missing <- as.data.frame(missingness)
missing$marker <- row.names(missing)

ggplot(missing) + geom_histogram(aes(x=missingness))+geom_vline(aes(xintercept= 0.2), color="red")+geom_vline(aes(xintercept= 0.1), color="blue")+theme_classic()+xlab("missingness (loci)")


range(genos_0.2$`%GT`)
##0 individuals have call rates below 70%%
##individual call rates range 95.47 - 99.72%

#remove markers with call rates below 50%
very_bad_markers <- missing[missing$missingness>0.5, 2]
print(paste(length(very_bad_markers), "markers with > 50% missing data"))

##1 marker with call rate below 50%
##Ots_wenYhap_33126##

##remove it##

genos_0.3 <- genos_0.2 %>%
  dplyr::select(-one_of(very_bad_markers))

##Recalculate IFI##

IFI <- marker_info %>%
  filter(marker %in% colnames(genos_0.3)) %>%
  group_by(ind) %>%
  summarize(back_count = sum(a1_count[called_geno == "A2HOM"], na.rm = TRUE)
            + sum(a2_count[called_geno == "A1HOM"], na.rm = TRUE)
            + sum(a1_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a2_count > a1_count)], na.rm = TRUE )
            + sum(a2_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a1_count > a2_count)], na.rm = TRUE ),
            
            hom_ct = sum(a1_count[called_geno == "A1HOM"], na.rm = TRUE)
            + sum(a2_count[called_geno == "A2HOM"], na.rm = TRUE)
            + sum(a2_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a2_count > a1_count)], na.rm = TRUE )
            + sum(a1_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a1_count > a2_count)], na.rm = TRUE ),
            
            ifi2 = (back_count/hom_ct)*100)



IFI$sample <- str_extract(IFI$ind, "[:upper:][:lower:]{2}[AJCU][RC]\\d{2}\\w{4}_\\d{4}")
IFI$adapter <- str_extract(IFI$ind, "[ATCG]{6}-[ATCG]{6}") 


genos_0.3 <- genos_0.3 %>%
  left_join(select(IFI, sample, adapter, ifi2), by = c("sample_simple" = "sample", "adapter" = "adapter")) %>%
  mutate(IFI = ifi2) %>%
  select(-one_of("ifi2"))

range(genos_0.3$IFI)
##0 samples have IFI scores greater than 10##

##recalculate individual call rates##

genos_0.3 %<>% #dont forget to change to the correct dataset
  rowwise %>%
  mutate(`%GT` = 100*(1-sum((cur_data() == "00"), na.rm = TRUE)/(ncol(.)-8))) %>% #change number
  ungroup() #always ungroup after rowwise


missing_ind<- which(genos_0.3$`%GT` < 90)
length(missing_ind)


##0 individuals have call rates below 90%##

genos_0.4<-genos_0.3

#recalculate locus level missingness

missingness3 <- (colSums(genos_0.4[,c(9:(ncol(genos_0.4)))] == "00" | genos_0.4[,c(9:(ncol(genos_0.4)))] == "0"))/nrow(genos_0.4) 
missing3 <- as.data.frame(missingness3)
missing3$marker <- row.names(missing3)


##identify markers missing more than 20% of individuals##

bad_markers <- missing3[missing3$missingness3>0.2, 2]

##1 markers missing more in more than 20% of individuals
##"Ots_RAS1"

#remove this marker##
genos_0.4 <- genos_0.4 %>%
  dplyr::select(-one_of(bad_markers))


##recalculate IFI

IFI <- marker_info %>%
  filter(marker %in% colnames(genos_0.4)) %>%
  group_by(ind) %>%
  summarize(back_count = sum(a1_count[called_geno == "A2HOM"], na.rm = TRUE)
            + sum(a2_count[called_geno == "A1HOM"], na.rm = TRUE)
            + sum(a1_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a2_count > a1_count)], na.rm = TRUE )
            + sum(a2_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a1_count > a2_count)], na.rm = TRUE ),
            
            hom_ct = sum(a1_count[called_geno == "A1HOM"], na.rm = TRUE)
            + sum(a2_count[called_geno == "A2HOM"], na.rm = TRUE)
            + sum(a2_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a2_count > a1_count)], na.rm = TRUE )
            + sum(a1_count[is.na(called_geno) == TRUE & ((a1_count + a2_count)>=10) & (a1_count > a2_count)], na.rm = TRUE ),
            
            ifi2 = (back_count/hom_ct)*100)

IFI$sample <- str_extract(IFI$ind, "[:upper:][:lower:]{2}[AJCU][RC]\\d{2}\\w{4}_\\d{4}")
IFI$adapter <- str_extract(IFI$ind, "[ATCG]{6}-[ATCG]{6}") 


genos_0.4 <- genos_0.4 %>%
  left_join(select(IFI, sample, adapter, ifi2), by = c("sample_simple" = "sample", "adapter" = "adapter")) %>%
  mutate(IFI = ifi2) %>%
  select(-one_of("ifi2"))

# now filter on IFI

high_IFI<-which(genos_0.4$IFI > 2.5)
length(high_IFI)


##0 individuals have IFI scores greater than 2.5

##recalculate marker call rates## 

missingness3 <- (colSums(genos_0.4[,c(9:(ncol(genos_0.4)))] == "00" | genos_0.4[,c(9:(ncol(genos_0.4)))] == "0"))/nrow(genos_0.4) 
missing3 <- as.data.frame(missingness3)
missing3$marker <- row.names(missing3)

##At this point the dataset has 86 Individuals and 350 loci and a sex marker##
##Individual call rates range 95.7 - 100%; Marker call rates range 89 - 100%
##IFI ranges 0.09 - 1.71##


##Evaluate possible paralogs##

#get marker names of markers with 0.1 > missingness > 0.2
miss0.1 <- missing3[missing3$missingness3 > 0.1,]
miss_mod <- miss0.1[miss0.1$missingness3 < 0.2, 2]


hets <- filter(marker_info, called_geno == "HET" | is.na(called_geno))

models <- hets %>%
  filter(marker %in% colnames(genos_0.4)) %>%
  filter(is.na(a1_count) == FALSE & is.na(a2_count) == FALSE) %>%
  group_by(marker) %>%
  group_map(~ lm(a1_count ~ a2_count, data= .))


lms <- lapply(models, coef)
ggplot()+geom_histogram(aes(x = sapply(lms,`[`,2)))+theme_classic()+ggtitle("allele ratios for all NA and HET calls")+geom_vline(aes(xintercept = 1.5), color = "red", linetype = 2)+geom_vline(aes(xintercept = (2/3)), color = "red", linetype = 2)+xlab("allele ratio (a1/a2)")+geom_vline(aes(xintercept = 1), color = "black")


lms_anova <- lapply(models, summary)


# collect info about each bad model
paralog_possible <- which(abs(sapply(lms,`[`,2)) > 1.5) #bad because a positively skewed allele ratio
paralog_possible2 <- which(abs(sapply(lms,`[`,2)) < (2/3)) # bad because a negative skewed allele ratio

paralog_possible3 <- which(sapply(lms_anova, function(x) x$coefficients[,4][2])> 0.01) # bad because too much variance in allele ratio, even if mean ratio is 1

paralog_possible <- c(paralog_possible, paralog_possible2, paralog_possible3)


plots <- marker_info %>%
  filter(marker %in% colnames(genos_0.4)) %>%
  filter(is.na(a1_count) == FALSE & is.na(a2_count) == FALSE) %>%
  group_by(marker) %>%
  do(plots=ggplot(data=.)+geom_point(aes(a1_count, a2_count, color = called_geno))+theme_classic()+geom_abline(aes(slope=1, intercept=0))+geom_abline(aes(slope = 10, intercept=0), color = "green")+geom_abline(aes(slope = 0.1, intercept=0), color = "red")+geom_abline(aes(slope = 0.2, intercept=0), color = "blue")+geom_abline(aes(slope = 5, intercept=0), color = "blue")+coord_equal(ratio=1)+geom_abline(slope = -1, intercept = 10)+ggtitle(unique(.$marker)))

#plot all "bad markers"

#first add the missningness markers to the list to examine
mod_bad_plot_index <- which(plots$marker %in% miss_mod)
paralog_possible <- c(mod_bad_plot_index, paralog_possible)

# then loop through the plots by changing the index until you have looked at all your questionable markers
plots$plots[[paralog_possible[10]]]

to_filt <- c("Ots_MetA")

##remove 1 locus that is possible paralog##
genos_0.5 <- genos_0.4 %>%
  dplyr::select(-one_of(to_filt))


##identify duplicated sampling events with relatedness estimates in coancestry
##format vcf for coancestry##


replacedGT<-as.data.frame(genos_0.5[, c(9:(ncol(genos_0.5)-1))]) ## -1 removes sex ID marker##

delim_gtypes <- replacedGT[, ] %>%
  mutate_all(funs(case_when(. == "AA" ~ "A/A",
                            . == "AT" ~ "A/T", 
                            . == "AC" ~ "A/C", 
                            . == "AG" ~ "A/G", 
                            . == "TT" ~ "T/T", 
                            . == "TA" ~ "T/A", 
                            . == "TC" ~ "T/C", 
                            . == "TG" ~ "T/G",
                            . == "CC" ~ "C/C", 
                            . == "CA" ~ "C/A", 
                            . == "CT" ~ "C/T", 
                            . == "CG" ~ "C/G", 
                            . == "GG" ~ "G/G", 
                            . == "GA" ~ "G/A", 
                            . == "GT" ~ "G/T", 
                            . == "GC" ~ "G/C", 
                            . == "--" ~ "-/-", 
                            . == "A-" ~ "A/-", 
                            . == "T-" ~ "T/-", 
                            . == "C-" ~ "C/-", 
                            . == "G-" ~ "G/-", 
                            . == "-A" ~ "-/A",
                            . == "-T" ~ "-/T", 
                            . == "-C" ~ "-/C", 
                            . == "-G" ~ "-/G", 
                            . == "00" ~ "0/0", 
                            . == "0" ~ "0/0"
  )))



library(splitstackshape)


vrtnames<-c(colnames(delim_gtypes))

##split genotypes into 2 alleles

dfspltgeno<-as.data.frame(cSplit(delim_gtypes, splitCols = vrtnames, sep = "/", type.convert = FALSE))


coancestrylbls<-data.frame(c(genos_0.5$sample_simple))

colnames(coancestrylbls)<-c("Ind")

coancestrylbls %<>%
  mutate(Ind = substr(Ind, 4, 16))


rownames(dfspltgeno)<-as.character(coancestrylbls$Ind)

##replace nucleotides with numeric values for coancestry## 
##A = 1; T = 2; C = 3; G = 4; - = 5 ###

dfspltgeno %<>%
  mutate_all(funs(case_when(. == "A" ~ 1,
                            . == "T" ~ 2,
                            . == "C" ~ 3,
                            . == "G" ~ 4, 
                            . == "-" ~ 5,
                            . == "0" ~ 0, 
  )))



write.table(dfspltgeno, file = "C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/OtsTRIR/OtsTRIR coancestry input.txt", quote = FALSE, sep = "\t", col.names = FALSE, row.names = TRUE)


##read in relatedness estimates##


rel<-read.table("C:/ZSL/Coancestry/OtsTRIR/RelatednessEstimates.txt", sep = ",")

##grab wang estimator##

r<-c(1:3, 6)

rel<-rel[, r]

colnames(rel)<-c("dyad", "Ind1", "Ind2", "wang")

##relatedness ranges -0.5 - 0.6##Thus no duplicates

##identify and remove monomorphic markers##

genos_2.0<-genos_0.5

##remove "dots" from column names##

colnames(genos_2.0) <- gsub("\\.", "_", colnames(genos_2.0))


genos_2.1 <- as.matrix(genos_2.0[,c(9:(ncol(genos_2.0)-1))])

row.names(genos_2.1) <- genos_2.0$sample_simple
genind_1.0 <- df2genind(genos_2.1, sep ="", ploidy=2,NA.char = "00")

##write genalex csv## 
genind2genalex(genind_1.0, filename = "C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/OtsTRIR/Ots TRIR genalex.csv")

##read in genalex file##
GT_poppr<-read.genalex("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/OtsTRIR/Ots TRIR genalex.csv")

locus_info<-locus_table(GT_poppr)

mono<-which(is.na(locus_info[, 4]))

##74 monomorphic markers##
##"Ots11_11925999"         "Ots11_32418659"         "Ots17_1066109_C6"       "Ots17_1345774_C6"       "Ots17_1486479_C6"       "Ots17_1488679_C6"       "Ots17_885364"           "Ots18_32088284"        
#"Ots19_46172133"         "Ots19_46172427"         "Ots28_11202190"         "Ots28_11205423"         "Ots28_11206740"         "Ots2_42405643"          "Ots33_19359879"         "Ots4_40942276"         
#"Ots4_41638710"          "Ots4_64978818"          "Ots5_44795073"          "Ots5_70908626"          "Ots7_51409415"          "Ots7_53631522"          "Ots7_54212944"          "Ots9_28975221"         
#"Ots_123048-521"         "Ots_123921-111"         "Ots_127760-569"         "Ots_129144-472"         "Ots_129870-55"          "Ots_97660-56"           "Ots_ARNT"               "Ots_CCR7"              
##"Ots_CHI06027687_143477" "Ots_CHI06048618_5222"   "Ots_GH2"                "Ots_GPH-318"            "Ots_GST-207"            "Ots_GST-375"            "Ots_GTH2B-550"          "Ots_GnRH-271"          
#"Ots_IGF-I_1-76"         "Ots_IL8R_C8"            "Ots_Ikaros-250"         "Ots_IsoT"               "Ots_LWSop-638"          "Ots_Myc-366"            "Ots_NFYB-147"           "Ots_SL"                
#"Ots_U2362-330"          "Ots_U2567-104"          "Ots_crRAD18937-60"      "Ots_crRAD23631-48"      "Ots_crRAD255-59"        "Ots_crRAD26081-28"      "Ots_crRAD46081-56"      "Ots_crRAD46751-42"     
#"Ots_crRAD55400-59"      "Ots_hnRNPL-533"         "Ots_hsc71-3'-488"       "Ots_nkef-192"           "Ots_nramp-321"          "Ots_pigh-105"           "Ots_stk6-516"           "Ots_tpx2-125"          
#"Ots_txnip-321"          "Ots_u07-20_332"         "Ots_u07-25_325"         "Ots_u07-64_221"         "Ots_u1008-108"          "Ots_u211-85"            "Ots_wenYhap_106664_9"   "Ots_wenYhap_25067_92"  
#"Ots_wenYhap_71572"      "Ots_zn593-346"


##Have to add 8 to column identifiers of monomorphic loci##
mm<-mono + 8

##remove monomorphic markers##

genos_2.0<-genos_2.0[, -mm]



##Dataset now has 407 individuals genotyped at 327 markers and a sex marker

ggplot(genos_2.0)+geom_density(aes(x=`On-Target Reads`))+geom_vline(aes(xintercept=median(`On-Target Reads`)), color = "red") +theme_classic()

ggplot(genos_2.0)+geom_density(aes(x=`%On-Target`))+geom_vline(aes(xintercept=median(`%On-Target`)), color = "red") +theme_classic()


##Depths##

marker_info %>%
  filter(marker %in% colnames(genos_2.0)) %>%
  filter(ind %in% genos_2.0$Sample) %>%
  mutate(sumdepth=a1_count+a2_count) %>%
  summarise(mean=mean(sumdepth, na.rm = TRUE), median=median(sumdepth, na.rm = TRUE), sd=sd(sumdepth, na.rm = TRUE))


marker_info %>%
  filter(marker %in% colnames(genos_2.0)) %>%
  filter(ind %in% genos_2.0$Sample) %>%
  mutate(sumdepth=a1_count+a2_count) %>%
  ggplot + aes(x=sumdepth)+geom_histogram()+theme_classic()+xlab("Mean Depth Per Locus Per Individual")


marker_info %>%
  filter(marker %in% colnames(genos_2.0)) %>%
  filter(ind %in% genos_2.0$Sample) %>%
  mutate(sumdepth=a1_count+a2_count) %>%
  ggplot + aes(x=sumdepth)+geom_histogram()+theme_classic()+xlab("Mean Depth Per Locus Per Individual")+xlim(0,1000)+ggtitle("inset of plot above from 0 to 1000 reads")


save(genos_2.0, file = "C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/OtsTRIR/Ots TRIR genos_2.0.R")


##export adult migration timing markers##

##read marker information for Ots GT-seq panel

library(readxl)

ots_panel_info<-read_excel("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/OtsAC22_23COQR/Ots GT-seq panel/Ots GT-seq panel info.xlsx", sheet = 2)

##subset panel info to those markers in the SFGL panel

ots_SFGL_panel_info<-ots_panel_info[ots_panel_info$Assay %in% colnames(genos_2.0), ]

##subset to run timing markers#
rt_markers<-subset(ots_SFGL_panel_info, ots_SFGL_panel_info$PresumedType == "runtiming")

run_timing<-genos_2.0[, c(1:8, which(colnames(genos_2.0) %in% rt_markers$Assay))]

##put markers in chromosomal order##

run_timing<-run_timing[, c(1:12, 34:36, 13:16, 37, 17:21, 38, 22:23, 39, 24:33)]


##write run timing markers to csv for polarization

write.csv(run_timing, quote = FALSE, row.names = FALSE, file = "C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/OtsTRIR/OtsTRIR_runtiming_genotypes.csv")


##load filtered genos to include genotypic sex##

load("C:/Users/olsenk2/OneDrive - Oregon State University/Desktop/OtsTRIR/Ots TRIR genos_2.0.R")
