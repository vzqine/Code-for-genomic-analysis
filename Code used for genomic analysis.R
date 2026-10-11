#==================== Install and load required packages ====================
packages <- c("tidyverse", "readxl","openxlsx","ggpubr")
install.packages(setdiff(packages, rownames(installed.packages())))
lapply(packages, library, character.only = TRUE)

#===== Modify to your actual working directory path =====
setwd()
getwd() 
rm(list=ls())
#==================== Read raw data (execute only once) ====================
df_raw <- read_excel("Mutation.all.xlsx",sheet = 1,col_names = T)

# Clean character‑type columns: trim leading/trailing whitespace
# Resolve matching failure caused by extra spaces in intergenic entries
df_raw <- df_raw %>% 
  mutate(across(c(snp_type, annotation_type, classification, Lineage, Group), str_squish))

#==================== Scheme B: retain intergenic regions; filter T0=0, max>=0.05 ====================
df3 <- df_raw%>%
  filter(T0==0) %>%              # Filter out mutations already present in the ancestor
  filter(max>=0.05)              # Filter mutations with maximum frequency below 5%

# Convert dataset to long format
dflong <- df3 %>%
  select(Lineage,Group,T0,T5,T10,T15,T20,position,Class,snp_type) %>%
  pivot_longer(names_to = "Transfer",
               values_to = "Freq",
               cols = starts_with("T"))

## Count non‑zero mutations per Lineage per Transfer (mutation position as unit)
mutsum <- dflong %>%
  filter(Freq!=0) %>%    
  group_by(Lineage,Group,Transfer) %>%
  count(position) %>%
  as.data.frame()

# Extract numeric index for lineage labels (compatible with R_1 / M_1 format)
mutsum$Lineage <- gsub("^[A-Z]_", "", mutsum$Lineage)

# Summarize total mutation number for each lineage at each time point
lineage <- mutsum %>%
  group_by(Lineage,Group,Transfer) %>%
  summarise(total=sum(n),.groups="drop") %>%
  mutate(Lineage = as.factor(Lineage),
         Transfer = factor(Transfer,levels=c("T0","T5","T10","T15","T20"))) %>%
  as.data.frame()

write.xlsx(lineage,"lineage_include_intergenic.xlsx")

# Line plot: dynamics of mutation number along transfers for each lineage
ggplot(lineage,aes(x = Transfer, y = total,group=Lineage,color=Lineage))+
  geom_line()+
  geom_point()+
  facet_wrap(~Group)+
  labs(x = "Transfer",y = "Number of mutations")+
  theme_bw()+
  theme(axis.text = element_text(size=6,color="black"),
        axis.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.position = "right",
        legend.title = element_text(size=8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.key.size = unit(0.9, "lines"),legend.margin = margin(0, 0, 0, 0))+
  ylim(0,25)

ggsave("LineageMutations_include_intergenic.png",width = 85,height = 55,units = "mm")
ggsave("LineageMutations_include_intergenic.pdf",width = 85,height = 55,units = "mm")


#==================== Stacked bar plot: statistics by snp_type ====================
SNP <- dflong %>%
  filter(Freq!=0) %>%
  group_by(Lineage,Group,Transfer) %>%
  count(snp_type) %>%
  as.data.frame()

ave <- SNP %>%
  group_by(Lineage,Group,snp_type) %>%
  summarize(ave = mean(n),.groups="drop") %>%
  as.data.frame()

ggplot(ave, aes(x = Lineage, y = ave, fill = snp_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Lineage",
       y = "Total mutations",
       fill = "Mutation Type") +
  facet_wrap(~Group, scales = "free_x")+
  theme_bw()+
  theme(axis.text.x = element_text(size=6,color="black",angle = 90),
        axis.text.y = element_text(size=6,color="black"),
        axis.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.position = "right",
        legend.title = element_text(size=8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.key.size = unit(0.9, "lines"),legend.margin = margin(0, 0, 0, 0)) +
  ylim(0,20)

ggsave("MutationType_Lineage_include_intergenic.png",width = 95,height = 55,units = "mm")
ggsave("MutationType_Lineage_include_intergenic.pdf",width = 95,height = 55,units = "mm")


#==================== Jaccard index for evolutionary parallelism ====================
df4 <- df_raw %>%
  filter(T0==0) %>%
  filter(max>=0.05)

df5 <- df4 %>%
  select(Lineage,Group,position) %>% 
  group_by(Lineage,Group) %>%
  count(position) %>%
  as.data.frame()

write.xlsx(df5,"position_include_intergenic.xlsx")

# Function to compute Jaccard index; add protection against division‑by‑zero when union length equals zero
calculate_jaccard <- function(df) {
  combinations <- combn(unique(df$Lineage), 2, simplify = FALSE)
  jaccard_indices <- lapply(combinations, function(pair) {
    lineage1 <- pair[1]
    lineage2 <- pair[2]
    positions_lineage1 <- df %>% filter(Lineage == lineage1) %>% distinct(position)
    positions_lineage2 <- df %>% filter(Lineage == lineage2) %>% distinct(position)
    
    intersect_length <- length(intersect(positions_lineage1$position, positions_lineage2$position))
    union_length <- length(union(positions_lineage1$position, positions_lineage2$position))
    jaccard_index <- ifelse(union_length == 0, 0, intersect_length / union_length)
    
    return(data.frame(Lineage1 = lineage1, Lineage2 = lineage2, JaccardIndex = jaccard_index))
  })
  return(bind_rows(jaccard_indices))
}

# Calculate Jaccard index separately for each experimental Group
jaccard_results <- df5 %>%
  group_by(Group) %>%
  nest() %>%
  mutate(jaccard_data = map(data, calculate_jaccard)) %>%
  select(-data) %>%
  unnest(jaccard_data)

write.xlsx(jaccard_results,"JaccardIndex_include_intergenic.xlsx")

# Boxplot of Jaccard index + t‑test between M and R groups
ggplot(jaccard_results, aes(x = Group, y = JaccardIndex, color = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  labs(x = "Group",
       y = "Jaccard similarity coefficient") +
  stat_compare_means(method = "t.test", comparisons = list(c("M", "R")), label = "p.signif") +
  theme_bw() +
  theme(axis.text = element_text(size=8, color="black"),
        axis.title = element_text(size=8),
        legend.position = "none",
        plot.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylim(0,0.35)

ggsave("JaccardIndex_include_intergenic.png",width = 35,height =55,units="mm")
ggsave("JaccardIndex_include_intergenic.pdf",width = 35,height =55,units="mm")