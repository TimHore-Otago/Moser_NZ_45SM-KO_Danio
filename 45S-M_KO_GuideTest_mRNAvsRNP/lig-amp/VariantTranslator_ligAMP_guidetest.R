# Variant_Translator 
# script variant for 45S-M_KO_29dpf_01

if (!require(readxl, quietly = TRUE)) {  install.packages("readxl")}
library(readxl)
if (!require(dplyr, quietly = TRUE)) {  install.packages("dplyr")}
library(dplyr)
if (!require(ggplot2, quietly = TRUE)) {  install.packages("ggplot2")}
library(ggplot2)
if (!require(tidyr, quietly = TRUE)) {  install.packages("tidyr")}
library(tidyr)
if (!require(ggpubr, quietly = TRUE)) {  install.packages("ggpubr")}
library(ggpubr)
if (!require(reshape2, quietly = TRUE)) {  install.packages("reshape2")}
library(reshape2)
if (!require(ggrepel, quietly = TRUE)) {  install.packages("ggrepel")}
library(ggrepel)

#set WD
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#import sample info
samples <- read.csv("45S-M_KO_guidetest.csv", header=TRUE) # list includes the PBAT samples not relevant to this script

#read the lists produced by the Variant_caller and Variant_lister shell scripts.
setwd("./Output")
file_names <- list.files(pattern = "^.*.txt$")
df_list <- list()       # Create an empty list to store individual data frames
for (file in file_names) {      # Import each file into a data frame and store it in the list
  file_content <- readLines(file)
  df <- data.frame(Reads = file_content, Source = file, stringsAsFactors = FALSE)
  df_list[[file]] <- df
}
df<- do.call(rbind.data.frame, c(df_list, make.row.names = FALSE))

#combine into a single dataframe
df$Run <- df$Source  #refers to the sequencing run this sample was on
df <- df %>%  separate(Run, into = c("Run"), sep = "-", remove = TRUE)
df_r1 <- as.data.frame(df[grepl("1:N:0", df$Reads),])
colnames(df_r1) <- c("R1","Source","Run")
df_r2 <- as.data.frame(df[grepl("2:N:0", df$Reads),])
colnames(df_r2) <- c("R2","Source","Run")
df_r1 <- df_r1 %>%  separate(R1, into = c("Reads",NA, "R1casVAR"), sep = " ", remove = TRUE)
df_r2 <- df_r2 %>%  separate(R2, into = c("Reads",NA, "R2casVAR"), sep = " ", remove = TRUE)
df <- left_join(df_r2,df_r1,by=c("Reads","Source","Run"))
df$casVARcode <- paste0(df$R1casVAR, "-", df$R2casVAR) # easier on the eye

df$ID <- gsub("_L001_conc_pr_on_rl_Output.txt", "", df$Source)
samples$ID <- gsub("_L001_R1_001.fastq", "", samples$filename)
samples <- samples[grep("_R1", samples$filename), ]
df <- inner_join(df,samples,by="ID") #should not increase obs beyond those in df




#quick clean
rm(list = setdiff(ls(),c( "df"))) #,"samples"

#translate the presence absence code into Indel Variants
df <- df %>%
  mutate(casVAR = case_when(
    casVARcode == "NNNNNYNYN-NNNYYNNNY" ~ "VII", 
    casVARcode == "NNNNNYNNN-NNNYYNNNN" ~ "IV+V",
    casVARcode == "YNNNNYNYN-YNNNNNNNN" ~ "V",
    casVARcode == "YNNNNYNYN-YNNNNNNNY" ~ "V",
    casVARcode == "NYYNNNNNN-NNNNNNNNN" ~ "II+III",
    casVARcode == "NYYNNYNNN-NNNNNNNNN" ~ "II+III",
    casVARcode == "NYNNNNNNN-NNNYYNNNN" ~ "II+III+IV",
    casVARcode == "NNNNNYNYN-NNNNYNNNY" ~ "IV+V",
    casVARcode == "NNYNNYNNN-NNNNYNNNN" ~ "IV+V",
    casVARcode == "NNNNNYNYN-NNNNYNNNN" ~ "IV+V",
    casVARcode == "NNNNNNNNN-NNNYYNYNN" ~ "II+IV", #yes, twice bigish localized indels push it across the gap
    casVARcode == "NYNNNNNNN-NNNYNNNNN" ~ "II+III", #makes sense      
    casVARcode == "NNNNNNNNN-NNNYYNNNN" ~ "III+IV", #makes sense                  
    casVARcode == "YNYNNNNNN-YNNYNNNNN" ~ "III", # III visible on both
    casVARcode == "YNYNNNNNN-YNNNNNNNN" ~ "III", # III onR1
    casVARcode == "YNNNNNNNN-YNNYNNNNN" ~ "III", # III onR2
    casVARcode == "NNNNNNNNN-NNNNNNNNN" ~ "unspecified",
    casVARcode == "NNNNNNNNN-NNNNNNYNN" ~ "VI",
    casVARcode == "NNNNNNNNN-NNNNYNNNN" ~ "IV",
    casVARcode == "NNNNNNNYN-NNNNNNNNY" ~ "VII",
    casVARcode == "NNNNNYNNN-NNNNNNNNN" ~ "V",
    casVARcode == "NYNNNNNNN-NNNNNNNNN" ~ "II",
    casVARcode == "NYNNNNNNN-NNNNNNYNN" ~ "II+VI",
    casVARcode == "NYNNNNNNN-NNNNYNNNN" ~ "II+IV",
    casVARcode == "YNNNNNNNN-YNNNNNNNN" ~ "I",
    casVARcode == "YNNNNNNNN-YNNNNNYNN" ~ "VI", # -> lopsided VI  (escapes smartpipe)
    casVARcode == "YNNNNNNNN-YNNNNNYNY" ~ "VI", # -> the VI del is so big, the perfectness filter for R1 flags both R1 and R2 as fine.(escapes smartpipe)
    casVARcode == "YNNNNYNNN-YNNNNNNNN" ~ "V",  # -> the V del is so big, the perfectness filter for R1 flags both R1 and R2 as fine.(escapes smartpipe)
    casVARcode == "NYNNNNNNN-NNNNNNYNY" ~ "II+VI",
    casVARcode == "NNNNNNNNN-NNNNYNYNN" ~ "IV", # -> IV so big it is flagged as being VI (escapes smartpipe)
    casVARcode == "NNNNNYNNN-NNNNYNNNN" ~ "IV+V",
    casVARcode == "NYNNNNNNN-NNNNYNYNN" ~ "II+IV",
    casVARcode == "NNNNNNNYN-NNNNNNNNN" ~ "VII", #one read too short to be flagged but thats fine
    casVARcode == "NNNNNNNNN-NNNNNNNNY" ~ "VII", #one read too short to be flagged but thats fine
    casVARcode == "YNNNNNNYN-YNNNNNNNY" ~ "VII",
    TRUE ~ NA_character_
  ))
df$casVAR <- ifelse(is.na(df$casVAR), "unspecified", df$casVAR)

# Minimum Read Depth
sample_counts <- table(df$Individual)
factors_to_keep <- names(sample_counts[sample_counts > 149])
df <- df[df$Individual %in% factors_to_keep, ]


# calculate variant proportions /sample
casVARprop <- df %>%  group_by(casVAR,treatment,Individual) %>%  summarize(Count = n()) %>%  group_by(Individual) %>%  mutate(Proportion = Count / sum(Count) * 100) %>%  ungroup()

casVARprop$casVAR <- as.factor(casVARprop$casVAR)

# Assign colors and order to indel variants
color_lookup <- data.frame(
  casVAR =     c("I",          "II",     "III",    "IV",      "II+III","II+III+IV", "II+IV", "III+IV",    "V",     "VI",      "IV+V",   "II+VI",    "VII", "unspecified"),
  color_code = c("#8EE5EE", "#986afc", "#c46afc", "#e129f2", "#7a6afc","#ff99cc", "#b56afc", "#ff6df2", "#ffb90f", "#fcb16a", "#cd950c", "#8b6508", "#dcde73", "#808080")
)
casVARprop <- merge(casVARprop, color_lookup, by = "casVAR", all.x = TRUE)
desired_order<-c("I",	"II",	"II+III",	"II+III+IV",	"II+IV",	"III",	"III+IV",	"IV",	"V",	"VI",	"IV+V",	"II+VI",	"VII",	"unspecified")
casVARprop$casVAR <- factor(casVARprop$casVAR, levels = desired_order)



#make a second data frame containing reference conformity proportion
preads <- aggregate(Count ~ Individual, data = casVARprop, sum)#n paired reads
ref_con <- subset(casVARprop,casVAR=="I")
ref_con <- rename(ref_con,pf_readprs_cVc = Count)
preads<- rename(preads,sum_readprs_cVc = Count)
ref_con <- left_join(ref_con,preads,by=c("Individual"))
ref_con <- rename(ref_con, ref_con = Proportion)


#tidy up
rm(list = setdiff(ls(),c( "ref_con", "casVARprop","df","color_lookup","samples")))
#x_order<-c("RNA1","RNA2","RNA3","RNA4","RNA5","RNA6","RNA7","RNA8",
#           "RNP1","RNP2","RNP3","RNP4","RNP5","RNP6","RNP7","RNP8","RNP9","RNP10","RNP11","RNP12",
#           "Uni1","Uni2","Uni3","Uni4")
#casVARprop$casVAR <- factor(casVARprop$casVAR, levels = desired_order)


# Plotting only reference conformity proportion
ggplot(subset(ref_con), aes(x = Individual, y = ref_con,fill = treatment))+ #, fill = treatment
  geom_col()+
  ylim(0,100)+
  theme(legend.position = "bottom")


# Plotting all indel variants
ggplot(casVARprop, aes(x = Individual, y = Proportion, fill = casVAR))+
  geom_col(position = position_stack(reverse = TRUE)) +
  facet_grid(cols = vars(treatment),scales = "free", space = "free")+
  ylab("Cas9 Indel variant Proportion") +
  scale_fill_manual(values = setNames(color_lookup$color_code, color_lookup$casVAR)) +
  theme_classic() +
  theme(legend.position = "bottom", legend.text = element_text(margin = margin(t = 0, r = 0, b = 0, l = -5)),
        legend.key.size = unit(0.4, "cm"), legend.key.width = unit(0.3, "cm")) +
  guides(fill = guide_legend(nrow = 1))
#theme(axis.title.x = element_blank(), axis.text.x = element_blank())

################################################################ old below

library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(reshape2)
library(ggrepel)

#read Sperm txts
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("./cReg_mRNA_vs_RNP/")
#iSeq037-TM-ZF-AB-RNP-02-n_S6_L001_Rboth_on_rl_Output.txt
file_names <- list.files(pattern = "^iSeq037-TM-ZF-AB.*L001_conc_pr_on_rl_Output.txt$")

df_list <- list()# Create an empty list to store individual data frames
# Import each file into a data frame and store it in the list
for (file in file_names) {
  file_content <- readLines(file)
  df <- data.frame(Reads = file_content, Source = file, stringsAsFactors = FALSE)
  df_list[[file]] <- df
}
# Combine all data frames into a single data frame
temp_df <- do.call(rbind.data.frame, c(df_list, make.row.names = FALSE))
df_r1 <- as.data.frame(temp_df[grepl("1:N:0", temp_df$Reads),])
colnames(df_r1) <- c("R1","Source")
df_r2 <- as.data.frame(temp_df[grepl("2:N:0", temp_df$Reads),])
colnames(df_r2) <- c("R2","Source")
df_r1 <- df_r1 %>%  separate(R1, into = c("Reads",NA, "R1casVAR"), sep = " ", remove = TRUE)
df_r2 <- df_r2 %>%  separate(R2, into = c("Reads",NA, "R2casVAR"), sep = " ", remove = TRUE)
df <- left_join(df_r2,df_r1,by=c("Reads","Source"))
#quick clean
#rm(temp_df,df_list,file_content,file,file_names,df_r1,df_r2,df)



### combine
#df <- rbind(fc,sp)
df$casVARcode <- paste0(df$R1casVAR, "-", df$R2casVAR) # easier on the eye

table(df$casVARcode)
#I_R1|II_R1|IIIr1_R1|IIIr2_R1|IV_R1|V_R1|VI_R1|VIIr1_R1|VIIr2_R1| - |I_R2|II_R2|IIIr1_R2|IIIr2_R2|IV_R2|V_R2|VI_R2|VIIr1_R2|VIIr2_R2

# it will also introduce new variants combos 
df <- df %>%
  mutate(casVAR = case_when(
    casVARcode == "NNNNNYNYN-NNNYYNNNY" ~ "VII", 
    casVARcode == "NNNNNYNNN-NNNYYNNNN" ~ "IV+V",
    casVARcode == "YNNNNYNYN-YNNNNNNNN" ~ "V",
    casVARcode == "YNNNNYNYN-YNNNNNNNY" ~ "V",
    casVARcode == "NYYNNNNNN-NNNNNNNNN" ~ "II+III",
    casVARcode == "NYYNNYNNN-NNNNNNNNN" ~ "II+III",
    casVARcode == "NYNNNNNNN-NNNYYNNNN" ~ "II+III+IV",
    casVARcode == "NNNNNYNYN-NNNNYNNNY" ~ "IV+V",
    casVARcode == "NNYNNYNNN-NNNNYNNNN" ~ "IV+V",
    casVARcode == "NNNNNYNYN-NNNNYNNNN" ~ "IV+V",
    casVARcode == "NNNNNNNNN-NNNYYNYNN" ~ "II+IV", #yes, twice bigish localized indels push it across the gap
    casVARcode == "NYNNNNNNN-NNNYNNNNN" ~ "II+III", #makes sense      
    casVARcode == "NNNNNNNNN-NNNYYNNNN" ~ "III+IV", #makes sense                  
    casVARcode == "YNYNNNNNN-YNNYNNNNN" ~ "III", # III visible on both
    casVARcode == "YNYNNNNNN-YNNNNNNNN" ~ "III", # III onR1
    casVARcode == "YNNNNNNNN-YNNYNNNNN" ~ "III", # III onR2
    casVARcode == "NNNNNNNNN-NNNNNNNNN" ~ "unspecified",
    casVARcode == "NNNNNNNNN-NNNNNNYNN" ~ "VI",
    casVARcode == "NNNNNNNNN-NNNNYNNNN" ~ "IV",
    casVARcode == "NNNNNNNYN-NNNNNNNNY" ~ "VII",
    casVARcode == "NNNNNYNNN-NNNNNNNNN" ~ "V",
    casVARcode == "NYNNNNNNN-NNNNNNNNN" ~ "II",
    casVARcode == "NYNNNNNNN-NNNNNNYNN" ~ "II+VI",
    casVARcode == "NYNNNNNNN-NNNNYNNNN" ~ "II+IV",
    casVARcode == "YNNNNNNNN-YNNNNNNNN" ~ "I",
    casVARcode == "YNNNNNNNN-YNNNNNYNN" ~ "VI", # -> lopsided VI  (escapes smartpipe)
    casVARcode == "YNNNNNNNN-YNNNNNYNY" ~ "VI", # -> the VI del is so big, the perfectness filter for R1 flags both R1 and R2 as fine.(escapes smartpipe)
    casVARcode == "YNNNNYNNN-YNNNNNNNN" ~ "V",  # -> the V del is so big, the perfectness filter for R1 flags both R1 and R2 as fine.(escapes smartpipe)
    casVARcode == "NYNNNNNNN-NNNNNNYNY" ~ "II+VI",
    casVARcode == "NNNNNNNNN-NNNNYNYNN" ~ "IV", # -> IV so big it is flagged as being VI (escapes smartpipe)
    casVARcode == "NNNNNYNNN-NNNNYNNNN" ~ "IV+V",
    casVARcode == "NYNNNNNNN-NNNNYNYNN" ~ "II+IV",
    casVARcode == "NNNNNNNYN-NNNNNNNNN" ~ "VII", #one read too short to be flagged but thats fine
    casVARcode == "NNNNNNNNN-NNNNNNNNY" ~ "VII", #one read too short to be flagged but thats fine
    casVARcode == "YNNNNNNYN-YNNNNNNNY" ~ "VII",
    TRUE ~ NA_character_
  ))

df$casVAR <- ifelse(is.na(df$casVAR), "unspecified", df$casVAR)

#extract source
df <- df %>%  separate(Source, into = c("Front","SampleID"), sep = "_", remove = TRUE)
df <- df %>%  separate(Front, into = c(NA,NA,NA,"line","trtmnt","repl"), sep = "-", remove = TRUE)
df$Sample <- paste(df$trtmnt,df$repl)
##
#casVARprop <- df %>%  group_by(Sire, casVAR, Tissue) %>%  summarize(Count = n()) %>%  group_by(Sire,Tissue) %>%  mutate(Proportion = Count / sum(Count) * 100) %>%  ungroup()
casVARprop <- df %>%  group_by(Sample, casVAR,trtmnt) %>%  summarize(Count = n()) %>%  group_by(Sample) %>%  mutate(Proportion = Count / sum(Count) * 100) %>%  ungroup()
casVARprop$casVAR <- as.factor(casVARprop$casVAR)

# Create a color lookup table with casVAR and corresponding color codes
color_lookup <- data.frame(
  casVAR =     c("I",          "II",     "III",    "IV",      "II+III","II+III+IV", "II+IV", "III+IV",    "V",     "VI",      "IV+V",   "II+VI",    "VII", "unspecified"),
  color_code = c("#8EE5EE", "#986afc", "#c46afc", "#e129f2", "#7a6afc","#ff99cc", "#b56afc", "#ff6df2", "#ffb90f", "#fcb16a", "#cd950c", "#8b6508", "#dcde73", "#808080"))
# Merge the color lookup table with casVARprop based on casVAR values
casVARprop <- merge(casVARprop, color_lookup, by = "casVAR", all.x = TRUE)
desired_order<-c("I",	"II",	"II+III",	"II+III+IV",	"II+IV",	"III",	"III+IV",	"IV",	"V",	"VI",	"IV+V",	"II+VI",	"VII",	"unspecified")
casVARprop$casVAR <- factor(casVARprop$casVAR, levels = desired_order)
#casVARprop$Sire <- sub("Sire", "#", casVARprop$Sire)
#casVARprop$SireTrt <- ifelse(as.numeric(gsub("#", "", casVARprop$Sire)) <= 21, paste0("i", gsub("#", "", casVARprop$Sire)), paste0("u", gsub("#", "", casVARprop$Sire)))
#casVARprop$SireTrt <- ifelse(as.numeric(gsub("#", "", casVARprop$Sire)) <= 21, paste0("", gsub("#", "", casVARprop$Sire)), paste0("", gsub("#", "", casVARprop$Sire)))
#casVARprop$trtmnt<- ifelse(as.numeric(gsub("#", "", casVARprop$Sire)) <= 21, "Injected", "Uninjected")

sum_readprs <- aggregate(Count ~ Sample, data = casVARprop, sum)
sum_readprs <-rename(sum_readprs,sum_readprs=Count)
#pairs_sp <- subset(sum_readprs, Tissue == "Sperm")
#pairs_fc <- subset(sum_readprs, Tissue == "FC")
cVp <- casVARprop
cVp <- left_join(cVp,sum_readprs,by = join_by(Sample))
rm(list = setdiff(ls(),c( "cVp","casVARprop","color_lookup","pairs_sp","pairs_fc" )))
cVp <- cVp %>%  mutate(my_alpha = ifelse(sum_readprs > 150, 1, 0)) #adjust second number to make low read cols transparent


# plots # plots # plots # plots # plots # plots  
ggplot(subset(cVp,sum_readprs > 150), aes(x = Sample, y = Proportion, fill = casVAR)) +
  geom_col(position = position_stack(reverse = TRUE),aes(alpha = I(my_alpha))) +
  scale_y_continuous(expand = c(0, 0), limits = c(-1, NA))+
  ylab("Proportions (Finclips)") +
  scale_fill_manual(values = setNames(color_lookup$color_code, color_lookup$casVAR)) +
  theme_classic() +
  labs(fill = "Indel Variant:")+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "bottom", legend.text = element_text(margin = margin(t = 0, r = 0, b = 0, l = -5)),
        legend.key.size = unit(0.4, "cm"), legend.key.width = unit(0.3, "cm")) +
  guides(fill = guide_legend(nrow = 1))+ #, reverse = TRUE
  facet_grid(cols = vars(trtmnt), scales = "free_x", space = "free_x",switch="both")+
  theme(strip.background = element_blank(),strip.text.x =  element_blank(),strip.placement = "outside")#+
#geom_text_repel(data = pairs_fc, aes(x = SireTrt, y = 20, label = sum_readprs), size = 3, direction = "y", segment.color = NA, inherit.aes = FALSE)


#just pf!
#simply subset casVar==I from casVARprop
#to make on onT filterable merge the result with preads
preads <- aggregate(Count ~ Sample, data = casVARprop, sum)#n paired reads
pf_cVc <- subset(casVARprop,casVAR=="I")
pf_cVc <- rename(pf_cVc,pf_readprs_cVc = Count)
preads<- rename(preads,sum_readprs_cVc = Count)
pf_cVc <- left_join(pf_cVc,preads,by=c("Sample"))
pf_cVc <- rename(pf_cVc, pf_cVc = Proportion)
pf_cVc$casVAR <-NULL
pf_cVc$color_code <-NULL
cVp <- casVARprop
cVp$color_code <-NULL
#rm(list = setdiff(ls(),c( "pf_cVc", "cVp")))
#pf_cVc <- pf_cVc %>%  separate(Source, into = c("Run",NA,NA,"sampleNR","exp","trtmnt","dpf"), sep = "-", remove = TRUE)

pf_cVc$moe_pf_cVc <- (2.326348*sqrt(((pf_cVc$pf_cVc/100)*(1-(pf_cVc$pf_cVc/100)))/pf_cVc$sum_readprs_cVc))*100
cdmTHEME <- theme_classic()+theme(plot.title = element_text(size=12),plot.subtitle = element_text(size=9),panel.border = element_blank(),legend.title = element_blank())
#ggplot(cREG2_29dpf_min, aes(x=uniqueident, y=i49i55_pf_cVc,fill=Tank)) 
ggplot(subset(pf_cVc,sum_readprs_cVc>150), aes(x = Sample, y = pf_cVc,fill=trtmnt))+
  geom_col(width=.7,show.legend = TRUE)+ #color="black"
  #labs(title = "cas9 induced fem-rDNA indels - 29dpf - Whole Fish")+
  scale_fill_manual(values = c("darkorange","darkgoldenrod1","cadetblue3"), name="Treatment", 
                    labels=c("mRNA", "RNP", "Uninjected"))+ #INj, MIX1 Mix2 MIx3 UNI
  xlab(" Individuals")+
  ylab("Ref. Conformity [%]")+
  theme_bw()+
  cdmTHEME+
  #annotate("text", x = cREG2_29dpf_min$uniqueident, y =cREG2_29dpf_min$X.PF_sp_i49and55, label="-")+ #change between pipelines
  theme(axis.text.x=element_blank(), axis.ticks.x= element_blank(),legend.title = element_text(), legend.position = "bottom",legend.text = element_text(margin = margin(t = 0, r = 0, b = 0, l = -5)))+#,legend.position = "bottom" +
  geom_errorbar(aes(x=Sample, ymin=pf_cVc-moe_pf_cVc, ymax=pf_cVc+moe_pf_cVc),linewidth=0.25, colour="#3e403f", alpha=1, width=0.7)


setwd("/Users/timmoser/Desktop/Google_Drive_@Work/C123456_Writeup/06_FemKO/Figures")
#write.csv(pf_cVc,"cR1_29_pf_cVc.csv") #perfectness based on casVAR caller
#write.csv(cVp,"cR1_29_cVp.csv") #casVAR proportions..
#rm(melted)
#ggarrange(left,ent,widths=c(2,0.6))

ggsave("mRNAvsRNP.pdf", 
       path = "/Users/timmoser/Desktop/Google_Drive_@Work/C123456_Writeup/06_FemKO/Figures",
       width = 5,
       height = 6.73,  #7 for one row
       units = c("cm"),
       dpi = 300,
       device = cairo_pdf(family="sans"))


(0.528301887) * log2(0.528301887) + (0.018867925) * log2(0.018867925) + (0.006289308) * log2(0.006289308) + (0.113207547) * log2(0.113207547) + (0.06918239) * log2(0.06918239) + (0.044025157) * log2(0.044025157) + (0.012578616) * log2(0.012578616) + (0.113207547) * log2(0.113207547) + (0.062893082) * log2(0.062893082) + (0.031446541) * log2(0.031446541)



#modstats

df <- subset(pf_cVc,sum_readprs_cVc>150)
aggregate(df$pf_cVc, by = list(df$trtmnt), 
          FUN = function(x) {
            med <- round(median(x), 2)
            n <- length(x)
            mean <- round(mean(x), 2)
            min_val <- round(min(x), 2)
            max_val <- round(max(x), 2)
            se <- round(sd(x)/sqrt(n), 2)
            c(median = med, n = n, mean = mean, min = min_val, max = max_val, se = se)
          })
