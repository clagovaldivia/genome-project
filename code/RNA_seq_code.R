library(readxl)
library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(forcats)
library(gridExtra)
library(cowplot)
library(patchwork)

#Load COG categories description

COG_description<- read_delim("~/Documents/RNA_htseq_2/COG_list", 
                             delim = "\t", escape_double = FALSE, 
                             col_names = FALSE, trim_ws = TRUE)
colnames(COG_description)<- c("COG_category","Description")
COG_description$COG_d<- paste(COG_description$COG_category, COG_description$Description, sep=":")

# List all the count file names for Sample A and Sample B

sampleA_files <- list.files("/home/claudia/Documents/RNA_htseq_2/DNA_egg/1", pattern = ".xlsx", full.names = TRUE)
sampleB_files <- list.files("/home/claudia/Documents/RNA_htseq_2/DNA_egg/2", pattern = "xlsx", full.names = TRUE)

### Load all SRR4342129 files ###

names <-gsub(".xlsx", "", sampleA_files)
merged_list <- list()
merged_29<- NULL
for(i in names){
  base_name<- sub("^([^.]+\\.[^.]+)\\..*$", "\\1" , basename(i))
  filepath <- as.data.frame(read_excel(paste("/home/claudia/Documents/RNA_htseq_2/DNA_egg/1/", basename(i), ".xlsx", sep=""), col_names = TRUE ))
  colnames(filepath)[1]<- "ID"
  assign(base_name, filepath)
  filelength <- as.data.frame(read_delim(paste("/home/claudia/Documents/RNA_htseq_2/RNA_htseq/1/",base_name,".genelengths",sep=""), delim = ";", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, skip = 2))
  filelength<- as.data.frame(lapply(filelength, function(x) gsub(".*=", "", x)))
  colnames(filelength)<- c("ID", "em_target", "em_score", "em_value","em_tcov","genelength")
  filelength<- as.data.frame(lapply(filelength, function(x) gsub("diamond\t", "", x)))
  filecount<- as.data.frame(read.delim(paste("/home/claudia/Documents/RNA_htseq_2/RNA_htseq/1/", base_name , ".txt", sep="") , header=FALSE))
  colnames(filecount)<- c("ID", "count")
  merged<-as.data.frame(merge(filecount, filelength , by = "ID"))
  all_merged<- merge(cbind(merged, "RPKM" = ((merged$count)*1000000000 ) / (as.numeric(merged$genelength)*sum(merged$count))), filepath, by= "ID")
  assign(paste("merged", base_name, sep="_"), all_merged )
  merged_list[[base_name]] <- all_merged
}

merged_list_29 <- lapply(merged_list, function(df) {
  df$evalue <- as.numeric(df$evalue)
  df$RPKM <- as.numeric(df$RPKM)
  return(df)
})

final_merged_29 <- merged_list_29 %>%
  bind_rows() %>%                  # Combine data frames
  group_by(ID) %>%                 # Group by "ID"
  summarize("COG_category" = first(COG_category),     # Keep the first occurrence of COG
            "KEGG_pathway" = first(KEGG_Pathway),
            "KEGG_ko" = first(KEGG_ko),
            "Description" = first(Description), # Keep the first occurrence of KEGG pathway
            RPKM = sum(RPKM, na.rm = TRUE)) %>%    # Sum "RPKM" for non-unique "ID" values 
  arrange(desc(RPKM))

final_final_29 <- final_merged_29 %>%
  group_by(COG_category) %>%               # Group by the "COG" column
  summarize(Total_RPKM = sum(RPKM)) %>%    # Sum "RPKM" for non-unique "ID" values 
  arrange(desc(Total_RPKM))

final_29<- merge(final_final_29, COG_description, by="COG_category")


### Load all SRR4342133 files ###

names <-gsub(".xlsx", "", sampleB_files)
merged_list <- list()
merged_33<- NULL
for(i in names){
  base_name<- sub("^([^.]+\\.[^.]+)\\..*$", "\\1" , basename(i))
  filepath <- as.data.frame(read_excel(paste("/home/claudia/Documents/RNA_htseq_2/DNA_egg/2/", basename(i), ".xlsx", sep=""), col_names = TRUE ))
  colnames(filepath)[1]<- "ID"
  assign(base_name, filepath)
  filelength <- as.data.frame(read_delim(paste("/home/claudia/Documents/RNA_htseq_2/RNA_htseq/2/",base_name,".genelengths",sep=""), delim = ";", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, skip = 2))
  filelength<- as.data.frame(lapply(filelength, function(x) gsub(".*=", "", x)))
  colnames(filelength)<- c("ID", "em_target", "em_score", "em_value","em_tcov","genelength")
  filelength<- as.data.frame(lapply(filelength, function(x) gsub("diamond\t", "", x)))
  filecount<- as.data.frame(read.delim(paste("/home/claudia/Documents/RNA_htseq_2/RNA_htseq/2/", base_name , ".txt", sep="") , header=FALSE))
  colnames(filecount)<- c("ID", "count")
  merged<-as.data.frame(merge(filecount, filelength , by = "ID"))
  all_merged<- merge(cbind(merged, "RPKM" = ((merged$count)*1000000000 ) / (as.numeric(merged$genelength)*sum(merged$count))), filepath, by= "ID")
  assign(paste("merged", base_name, sep="_"), all_merged )
  merged_list[[base_name]] <- all_merged
}

merged_list_33 <- lapply(merged_list, function(df) {
  df$evalue <- as.numeric(df$evalue)
  df$RPKM <- as.numeric(df$RPKM)
  return(df)
})

final_merged_33 <- merged_list_33 %>%
  bind_rows() %>%                  # Combine data frames
  group_by(ID) %>%                 # Group by "ID"
  summarize("COG_category" = first(COG_category),     # Keep the first occurrence of COG
            "KEGG_pathway" = first(KEGG_Pathway),
            "KEGG_ko" = first(KEGG_ko),
            "Description" = first(Description), # Keep the first occurrence of KEGG pathway
            RPKM = sum(RPKM, na.rm = TRUE)) %>%    # Sum "RPKM" for non-unique "ID" values 
  arrange(desc(RPKM))

final_final_33 <- final_merged_33 %>%
  group_by(COG_category) %>%               # Group by the "COG" column
  summarize(Total_RPKM = sum(RPKM)) %>%    # Sum "RPKM" for non-unique "ID" values 
  arrange(desc(Total_RPKM))

final_33<- merge(final_final_33, COG_description, by="COG_category")


# Create a histogram with COG

p1<- final_29[1:20,] %>% 
  mutate(COG_category = reorder(COG_category, Total_RPKM)) %>%  
  ggplot(aes(x=COG_category, y=Total_RPKM, fill=COG_d)) +
  geom_bar(stat = "identity")+
  ggtitle("D1")+
  coord_flip()+
  theme(legend.position = "right")

p2<- final_33[1:20,] %>% 
  mutate(COG_category = reorder(COG_category, Total_RPKM)) %>%  
  ggplot(aes(x=COG_category, y=Total_RPKM, fill=COG_d)) +
  geom_bar(stat = "identity")+
  ggtitle("D3")+
  coord_flip()+
  theme(legend.position = "right")


# Print the combined plot

combined_plot <- p1 / p2
combined_plot
ggsave("~/Documents/RNA_htseq_2/combined_COG_plot.png", combined_plot, width = 20, height = 20)

q1<- final_merged_29[1:15,] %>% 
  mutate(COG_category = reorder(COG_category, RPKM)) %>%  
  ggplot(aes(x=COG_category, y=RPKM, fill=Description)) +
  geom_bar(stat = "identity")+
  ggtitle("D1")+
  coord_flip()+
  theme(legend.direction = "vertical",
        legend.box = "horizontal",
        legend.position = "bottom")

ggsave("~/Documents/RNA_htseq_2/genes_29_plot.png", q1, width = 20, height = 20)


q2<- final_merged_33[1:15,] %>% 
  mutate(COG_category = reorder(COG_category, RPKM)) %>%  
  ggplot(aes(x=COG_category, y=RPKM, fill=Description)) +
  geom_bar(stat = "identity")+
  ggtitle("D3")+
  coord_flip()+
  theme(legend.direction = "vertical",
        legend.box = "horizontal",
        legend.position = "bottom")

ggsave("~/Documents/RNA_htseq_2/genes_33_plot.png", q2, width = 20, height = 20)

