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
samples <- read.csv("45S-M_KO_scriptTest.csv", header=TRUE) # list includes the PBAT samples not relevant to this script

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
df$ID <- gsub("_L001_pr_on_rl_Output.txt", "", df$Source)
samples$ID <- gsub("_L001_R1_001.fastq", "", samples$filename)
samples <- samples[grep("lh-amp", samples$protocol), ] # remove PBAT files from list
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


# Plotting only reference conformity proportion
ggplot(subset(ref_con), aes(x = Individual, y = ref_con,fill = treatment))+ #, fill = treatment
  geom_col()+
  ylim(0,100)+
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


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

