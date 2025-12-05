setwd("F:/Pathoviewer data/Emma_seaweed")

library(car)
library(multcomp)
library(stats)
library(agricolae)
library(ExpDE)
library(multcompView)
library(pgirmess)
library(TukeyC)
library(cluster)
library(e1071)
library(gplots)
library(MASS)
library(leaps)
library(glmnet)
library(lars)
library(rpart)
library("plyr")
library(openxlsx)
library(ggplot2)
library("ggpubr")
library(psych)
library(ggpubr)
library(gmodels)
library(plotrix)
library(wesanderson)
library(RColorBrewer)
library("devtools")
library("FactoMineR")
library("factoextra")
library("dplyr") 
library(Rmisc)

#### Data input ####

data=read.xlsx("F:/Pathoviewer data/Emma_seaweed/Data_Analysis_Report_20251205_1214.xlsx",
               startRow=174, colNames=TRUE, 
               skipEmptyRow=TRUE, rowNames=FALSE)
colnames(data) <- gsub("Fv/Fm", "Fv.Fm", colnames(data))
data <- data[, c("File", "Roi.No","Roi.Mask","Fv.Fm",
                 "Fv/fm.I.(%)","Fv/fm.II.(%)","Fv/fm.III.(%)","Fv/fm.IV.(%)","Fv/fm.V.(%)", 
                 "ChlIdx", "Chl.I.(%)","Chl.II.(%)","Chl.III.(%)","Chl.IV.(%)","Chl.V.(%)", 
                 "AriIdx","NDVI",
                 "Gfp.I.(%)","Gfp.II.(%)","Gfp.III.(%)","Gfp.IV.(%)","Gfp.V.(%)","Comment")]
library(stringr)
### dpi
dpi_pattern <- "(2day|3day|5day|4day)"

data <- data %>%
  mutate(DPI = str_extract(Comment, dpi_pattern)) %>%
  mutate(DPI = str_replace(DPI, "5day", "1day")) #correct the timepoint
data<-data %>%
  mutate(DPI=case_when(
    DPI=="1day" ~ "1 DPI",
    DPI=="2day" ~ "2 DPI",
    DPI=="3day" ~ "3 DPI",
    DPI=="4day" ~ "4 DPI",
    TRUE ~ DPI
  ))


### Strain
treatment_pattern <- "(ep1|ep2|ep3|ep4)"
data <- data %>%
  mutate(Strain = str_extract(Comment, treatment_pattern)) %>%
  mutate(Strain = case_when(
    Strain == "ep1" ~ "Strain 1",
    Strain == "ep2" ~ "Strain 2",
    Strain == "ep3" ~ "Strain 3",
    Strain == "ep4" ~ "Strain 4",
    TRUE ~ Strain
  ))

#Clean remove ROI with all
data <-subset(data, Roi.No != "All")



# Order Treatment names

data$Strain<-factor(data$Strain,
                       levels=c("Strain 1","Strain 2","Strain 3","Strain 4"))


#### Roi.Mask ####
#### Statistics for Roi.Mask ----

summary_stats <- data %>%
  group_by(Strain, DPI) %>%
  summarise(
    n = n(),
    mean = mean(Roi.Mask, na.rm = TRUE),
    sd = sd(Roi.Mask, na.rm = TRUE),
    se = sd / sqrt(n)
  )


data$Strain <- as.factor(data$Strain)


# --- Function to test per DPI ---
analyze_dpi <- function(df, dpi_value) {
  df <- df[df$DPI == dpi_value, ]
  df$Strain <- droplevels(df$Strain)
  
  if (nlevels(df$Strain) < 2) {
    message(paste("⚠️ Skipping:", dpi_value, "— not enough groups"))
    return(NULL)
  }
  
  # --- Normality & homogeneity ---
  aov_model <- aov(Roi.Mask ~ Strain, data = df)
  shapiro_p <- tryCatch(shapiro.test(residuals(aov_model))$p.value, error = function(e) NA)
  levene_p <- tryCatch(car::leveneTest(Roi.Mask ~ Strain, data = df)$`Pr(>F)`[1], error = function(e) NA)
  
  message(paste0("DPI: ", dpi_value,
                 " | Shapiro p = ", round(shapiro_p, 3),
                 " | Levene p = ", round(levene_p, 3)))
  
  # --- Choose test ---
  if (!is.na(shapiro_p) && !is.na(levene_p) && shapiro_p > 0.05 & levene_p > 0.05) {
    # ✅ ANOVA + Tukey
    tukey_res <- TukeyHSD(aov_model)
    letters <- multcompView::multcompLetters4(aov_model, tukey_res)
    letter_df <- data.frame(
      Strain = names(letters$Strain$Letters),
      Letters = letters$Strain$Letters,
      Test = "ANOVA"
    )
  } else {
    # ⚠️ Kruskal–Wallis + Dunn
    kw <- kruskal.test(Roi.Mask ~ Strain, data = df)
    dunn_res <- suppressWarnings(FSA::dunnTest(Roi.Mask ~ Strain, data = df, method = "none"))
    PTD <- dunn_res$res
    
    letter_df <- rcompanion::cldList(
      comparison = PTD$Comparison,
      p.value = PTD$P.adj,
      threshold = 0.05
    ) %>%
      dplyr::rename(Strain = Group, Letters = Letter)
    letter_df$Test <- "Kruskal-Wallis"
  }
  
  # --- Add y-position for label placement ---
  letter_df$y_pos <- tapply(df$Roi.Mask, df$Strain, max, na.rm = TRUE)[letter_df$Strain] * 1.05
  letter_df$DPI <- dpi_value
  return(letter_df)
}

# --- Run manually for specific DPI subsets ---
letter_1dpi <- analyze_dpi(data, "1 DPI")
letter_2dpi <- analyze_dpi(data, "2 DPI")
letter_3dpi <- analyze_dpi(data, "3 DPI")
letter_4dpi <- analyze_dpi(data, "4 DPI")

# --- Combine results ---
letters_df <- rbind(letter_1dpi,letter_2dpi, letter_3dpi,letter_4dpi)
letters_df

# Creating the plot with the filtered data
#### Box plot ####
win.graph(w=5,h=4)
ggplot(data, aes(x=Strain, y=Roi.Mask, fill=Strain)) +
  geom_boxplot(outlier.size=0.8, width=0.5, position=position_dodge(0.6)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 1.5) +
  geom_text(data = letters_df,
            aes(x = Strain, y = y_pos, label = Letters),
            size = 4, fontface = "bold", vjust = 0) +
  xlab("") +
  ylab("Fungal growth (px)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1.5), face = "bold", vjust = 1.5),
    axis.title = element_text(face = "bold"),
    axis.title.y = element_text(vjust = 1.8),
    axis.text.y = element_text(colour = "black"),
    axis.title.x = element_text(vjust = -0.5),
    axis.text.x = element_text(angle=45, size=10, hjust=0.5, vjust=0.6, colour="black"),
    panel.border = element_rect(colour = "black", fill=NA),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(), 
    legend.position = "none", 
    legend.direction = "vertical"
  ) +
  facet_grid(~DPI)
 # scale_fill_manual(values=c("#c8a5e5", "#f5c380"))



savePlot(filename ="F:/Pathoviewer data/Emma_seaweed/Roi.Mask.pdf",
         type = c("pdf"),
         device = dev.cur(),
         restoreConsole = TRUE)


#### Linechart ####
library(dplyr)
library(ggplot2)

# summarize means for line chart
summary_df <- data %>%
  group_by(DPI, Strain) %>%
  summarise(
    mean_growth = mean(Roi.Mask, na.rm = TRUE),
    se = sd(Roi.Mask, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
win.graph(w=5,h=4)
ggplot(summary_df, aes(x = DPI, y = mean_growth, color = Strain, group = Strain)) +
  geom_line(linewidth = 1) +  
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_growth - se, ymax = mean_growth + se), width = 0.1) +
  ylab("Fungal growth (px)") +
  xlab("") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size=10, colour="black"),
    axis.text.y = element_text(size=10, colour="black"),
    axis.title = element_text(face="bold"),
    legend.title = element_blank(),
    panel.border = element_rect(colour="black", fill=NA)
  )

#### Fv/Fm ####
win.graph(w=5,h=4)
ggplot(data, aes(x=Strain, y=Fv.Fm, fill=Strain)) +
  geom_boxplot(outlier.size=0.8, width=0.5, position=position_dodge(0.6)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 1.5) +
  xlab("") +
  ylab("Chlorophyl Fluorescence (Fv/Fm)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1.5), face = "bold", vjust = 1.5),
    axis.title = element_text(face = "bold"),
    axis.title.y = element_text(vjust = 1.8),
    axis.text.y = element_text(colour = "black"),
    axis.title.x = element_text(vjust = -0.5),
    axis.text.x = element_text(angle=45, size=10, hjust=0.5, vjust=0.6, colour="black"),
    panel.border = element_rect(colour = "black", fill=NA),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(), 
    legend.position = "none", 
    legend.direction = "vertical"
  ) +
  facet_grid(~DPI)


#### ChlIdx ####
win.graph(w=5,h=4)
ggplot(data, aes(x=Strain, y=ChlIdx, fill=Strain)) +
  geom_boxplot(outlier.size=0.8, width=0.5, position=position_dodge(0.6)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 1.5) +
  xlab("") +
  ylab("ChlIdx") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1.5), face = "bold", vjust = 1.5),
    axis.title = element_text(face = "bold"),
    axis.title.y = element_text(vjust = 1.8),
    axis.text.y = element_text(colour = "black"),
    axis.title.x = element_text(vjust = -0.5),
    axis.text.x = element_text(angle=45, size=10, hjust=0.5, vjust=0.6, colour="black"),
    panel.border = element_rect(colour = "black", fill=NA),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(), 
    legend.position = "none", 
    legend.direction = "vertical"
  ) +
  facet_grid(~DPI)

#### mARI ####
win.graph(w=5,h=4)
ggplot(data, aes(x=Strain, y=AriIdx, fill=Strain)) +
  geom_boxplot(outlier.size=0.8, width=0.5, position=position_dodge(0.6)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 1.5) +
  xlab("") +
  ylab("AriIdx") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1.5), face = "bold", vjust = 1.5),
    axis.title = element_text(face = "bold"),
    axis.title.y = element_text(vjust = 1.8),
    axis.text.y = element_text(colour = "black"),
    axis.title.x = element_text(vjust = -0.5),
    axis.text.x = element_text(angle=45, size=10, hjust=0.5, vjust=0.6, colour="black"),
    panel.border = element_rect(colour = "black", fill=NA),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(), 
    legend.position = "none", 
    legend.direction = "vertical"
  ) +
  facet_grid(~DPI)

#### NDVI ####
win.graph(w=5,h=4)
ggplot(data, aes(x=Strain, y=NDVI, fill=Strain)) +
  geom_boxplot(outlier.size=0.8, width=0.5, position=position_dodge(0.6)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6), size = 1.5) +
  xlab("") +
  ylab("NDVI") +
  theme_bw() +
  theme(
    plot.title = element_text(size = rel(1.5), face = "bold", vjust = 1.5),
    axis.title = element_text(face = "bold"),
    axis.title.y = element_text(vjust = 1.8),
    axis.text.y = element_text(colour = "black"),
    axis.title.x = element_text(vjust = -0.5),
    axis.text.x = element_text(angle=45, size=10, hjust=0.5, vjust=0.6, colour="black"),
    panel.border = element_rect(colour = "black", fill=NA),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(), 
    legend.position = "none", 
    legend.direction = "vertical"
  ) +
  facet_grid(~DPI)
#### FvFm class ####
colnames(data)[5:9] <- c("Class I", 
                         "Class II",
                         "Class III", 
                         "Class IV",
                         "Class V")

library(tidyr)
library(dplyr)
new_data <- data %>%
  pivot_longer(
    cols = starts_with("Class"),
    names_to = "Class",
    values_to = "Value"
  )
new_data$Value<-as.numeric(new_data$Value)
new_data<-new_data[,c(2,19:22)]
new_data_summary <- new_data %>%
  group_by(Strain, DPI, Class,Roi.No) %>%
  summarise(Value = mean(Value), .groups = "drop")

win.graph(w=4, h=6)
ggplot(new_data_summary, aes(x = Strain, y = Value, fill = Class)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab("") +
  ylab("Proportion (%)") +
  scale_fill_manual(values = c("#ff0000", "#ff8000", "#ffff00", "#bfff00", "#33cc33")) +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, by = 20)) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.title = element_text(size = rel(1.5), face = "bold", vjust = 1.5),
    axis.title = element_text(face = "bold"),
    axis.title.y = element_text(vjust = 1.8),
    axis.text.y = element_text(colour = "black"),
    axis.title.x = element_text(vjust = -0.5),
    axis.text.x = element_text(angle = 45, size = 10, face = "bold",
                               hjust = 0.5, vjust = 0.6, colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  facet_grid(Roi.No~DPI)


savePlot(filename ="F:/Pathoviewer data/Emma_seaweed/Fv.Fm classes.pdf",
         type = c("pdf"),
         device = dev.cur(),
         restoreConsole = TRUE)
