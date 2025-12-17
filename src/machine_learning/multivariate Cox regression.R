library(dplyr)
library(stringr)
library(survival)
setwd("C:/Users/10027/Desktop/data/GSE39582")
GSE39582_exp <- read.csv("GSE39582/GSE39582_exp_data.csv",row.names = 1)
GSE39582_cli <- read.csv("GSE39582/GSE39582_cli.csv",row.names = 1)
GSE39582_cli_in <- GSE39582_cli[c("Sex","age.at.diagnosis..year.","characteristics_ch1.4",
                                  "os.event","os.delay..months.")]
colnames(GSE39582_cli_in) <- c("gender","age","stage","os_event","os_months")
best_gene <- c("ADAP2","ATP6AP2","BRI3","C15orf48","CD163","CD36","CD68","ENG","F13A1","FABP5","FOLR2","LRPAP1","MMP12","MMP14","MRC1","MSR1","NPL","NRP1","PDK4","PMP22","SPP1","VSIG4","NPC1","NPC2")

GSE39582_exp_in <- GSE39582_exp[rownames(GSE39582_exp) %in% best_gene,]
GSE39582_exp_in <- as.data.frame(t(GSE39582_exp_in))
mydata <- cbind(GSE39582_cli_in,GSE39582_exp_in)
mydata$age <- as.numeric(mydata$age)
mydata$os_months <- as.numeric(mydata$os_months)
mydata <- mydata[!is.na(mydata$os_months),]
mydata$os_event <- as.numeric(mydata$os_event)
mydata$gender <- gsub(" ","",mydata$gender)
mydata$stage <- gsub("tnm.stage: ","",mydata$stage)
mydata$stage <- as.numeric(mydata$stage)
table(mydata$gender)

mydata2 <- mydata[,1:5]
mydata2$gene_mean <- rowMeans(mydata[, c(6:29)], na.rm = TRUE)

formula <- as.formula(
  paste("Surv(os_months, os_event) ~", paste(best_gene, collapse = " + "))
)
cox_multi_GEO_gene <- coxph(formula, data = mydata)
cox_multi_GEO_gene_cli <- coxph(Surv(os_months, os_event) ~., data = mydata)
cox_multi_GEO_gene_mean_cli <- coxph(Surv(os_months, os_event) ~., data = mydata2)
summary(cox_multi_GEO_gene_mean_cli)

library(forestplot)
library(dplyr)

cox_res <- summary(cox_multi)

df_forest <- data.frame(
  Variable = rownames(cox_res$coefficients),
  HR = round(cox_res$coefficients[, "exp(coef)"], 2),
  lower = round(cox_res$conf.int[, "lower .95"], 2),
  upper = round(cox_res$conf.int[, "upper .95"], 2),
  P = signif(cox_res$coefficients[, "Pr(>|z|)"], 3)
)

df_forest

tabletext <- cbind(
  c("Variable", df_forest$Variable),
  c("HR (95% CI)", paste0(df_forest$HR, " (", df_forest$lower, "-", df_forest$upper, ")")),
  c("P value", df_forest$P)
)

pdf("加临床信息的均值表达多cox.pdf",onefile=T)
p <- forestplot(
  labeltext = tabletext,
  mean = c(NA, df_forest$HR),
  lower = c(NA, df_forest$lower),
  upper = c(NA, df_forest$upper),
  zero = 1,
  boxsize = 0.3,
  line.margin = 0.2,
  xlab = "Hazard Ratio"
)
print(p)
dev.off()




# TCGA --------------------------------------------------------------------

setwd("../TCGA/")
TCGA_CRC_cli <- read.csv("TCGA_CRC_Clidata.csv",row.names = 1)
TCGA_CRC_cli_in <- TCGA_CRC_cli[c("gender","pathologic_stage","age_at_initial_pathologic_diagnosis")]
TCGA_CRC_cli_in$pathologic_stage[grepl("^Stage I([ABC])?$",  TCGA_CRC_cli_in$pathologic_stage)]   <- 1
TCGA_CRC_cli_in$pathologic_stage[grepl("^Stage II([ABC])?$", TCGA_CRC_cli_in$pathologic_stage)]  <- 2
TCGA_CRC_cli_in$pathologic_stage[grepl("^Stage III([ABC])?$",TCGA_CRC_cli_in$pathologic_stage)]  <- 3
TCGA_CRC_cli_in$pathologic_stage[grepl("^Stage IV([ABC])?$", TCGA_CRC_cli_in$pathologic_stage)]  <- 4
TCGA_CRC_cli_in$pathologic_stage <- as.numeric(TCGA_CRC_cli_in$pathologic_stage)
TCGA_CRC_cli_in$age_at_initial_pathologic_diagnosis <- as.numeric(TCGA_CRC_cli_in$age_at_initial_pathologic_diagnosis)
colnames(TCGA_CRC_cli_in) <- c("gender","stage","age")
# a <- read.csv("TCGA_CRC_Matadata.csv",row.names = 1)
TCGA_CRC_exp <- read.csv("D:/keti/keti_cell_type/single_cancer/TCGA_tpm/TCGA_COAD_TPM.csv",row.names = 1)
TCGA_gene <- c("SPP1","APOE","CD163","TREM2","GPNMB")
TCGA_CRC_exp <- as.data.frame(t(TCGA_CRC_exp))
TCGA_CRC_exp_in <- TCGA_CRC_exp[TCGA_gene,]
TCGA_CRC_exp_in <- as.data.frame(t(TCGA_CRC_exp_in))
TCGA_CRC_exp_in$sample <- substr(rownames(TCGA_CRC_exp_in),1,12)
TCGA_CRC_exp_in$type <- substr(rownames(TCGA_CRC_exp_in),14,16)
TCGA_CRC_exp_in <- TCGA_CRC_exp_in[TCGA_CRC_exp_in$type=="01A",]
# TCGA_CRC_exp_in$sample <- gsub("\\.","-",TCGA_CRC_exp_in$sample)
rownames(TCGA_CRC_exp_in) <- TCGA_CRC_exp_in$sample
TCGA_CRC_exp_in <- TCGA_CRC_exp_in[,-c(6,7)]
TCGA_CRC_exp_in[] <- lapply(TCGA_CRC_exp_in, as.numeric)
# TCGA_CRC_exp_in <- log(TCGA_CRC_exp_in+1)

TCGA_CRC_sur <- read.delim("D:/keti/keti_cell_type/single_cancer/TCGA-survial/TCGA-COAD.survival.tsv")
TCGA_CRC_sur <- TCGA_CRC_sur[,-1]
TCGA_CRC_sur <- TCGA_CRC_sur[!duplicated(TCGA_CRC_sur), ]
rownames(TCGA_CRC_sur) <- TCGA_CRC_sur$X_PATIENT
common_id <- intersect(rownames(TCGA_CRC_cli_in),rownames(TCGA_CRC_exp_in))
TCGA_mydata <- cbind(TCGA_CRC_cli_in[common_id, ],TCGA_CRC_exp_in[common_id, ])
common_id <- intersect(rownames(TCGA_CRC_sur),rownames(TCGA_mydata))
TCGA_mydata <- cbind(TCGA_mydata[common_id, ],TCGA_CRC_sur[common_id, ])
TCGA_mydata <- TCGA_mydata[,-10]

TCGA_mydata2 <- TCGA_mydata[,c(1:3,9,10)]
TCGA_mydata2$gene_mean <- rowMeans(TCGA_mydata[, c(4:8)], na.rm = TRUE)

formula <- as.formula(
  paste("Surv(OS.time, OS) ~", paste(TCGA_gene, collapse = " + "))
)

cox_multi_GEO_gene_mean_cli <- coxph(Surv(os_months, os_event) ~., data = mydata2)

cox_multi_TCGA_gene <- coxph(formula, data = TCGA_mydata)
cox_multi_TCGA_gene_cli <- coxph(Surv(OS.time, OS) ~ ., data = TCGA_mydata)
cox_multi_TCGA_gene_mean_cli <- coxph(Surv(OS.time, OS) ~., data = TCGA_mydata2)
summary(cox_multi_TCGA_gene_mean_cli)


library(forestplot)
library(dplyr)

cox_res <- summary(cox_multi)

df_forest <- data.frame(
  Variable = rownames(cox_res$coefficients),
  HR = round(cox_res$coefficients[, "exp(coef)"], 2),
  lower = round(cox_res$conf.int[, "lower .95"], 2),
  upper = round(cox_res$conf.int[, "upper .95"], 2),
  P = signif(cox_res$coefficients[, "Pr(>|z|)"], 3)
)

df_forest

tabletext <- cbind(
  c("Variable", df_forest$Variable),
  c("HR (95% CI)", paste0(df_forest$HR, " (", df_forest$lower, "-", df_forest$upper, ")")),
  c("P value", df_forest$P)
)

pdf("加临床信息的基因均值多cox.pdf",onefile=T)
p <- forestplot(
  labeltext = tabletext,
  mean = c(NA, df_forest$HR),
  lower = c(NA, df_forest$lower),
  upper = c(NA, df_forest$upper),
  zero = 1,
  boxsize = 0.3,
  line.margin = 0.2,
  xlab = "Hazard Ratio"
)
print(p)
dev.off()


# 保存一下变量 ------------------------------------------------------------------
cox_multi_TCGA_gene <- coxph(formula, data = TCGA_mydata)
cox_multi_TCGA_gene_cli <- coxph(Surv(OS.time, OS) ~ ., data = TCGA_mydata)
cox_multi_TCGA_gene_mean_cli <- coxph(Surv(OS.time, OS) ~., data = TCGA_mydata2)
cox_multi_GEO_gene <- coxph(formula, data = mydata)
cox_multi_GEO_gene_cli <- coxph(Surv(os_months, os_event) ~., data = mydata)
cox_multi_GEO_gene_mean_cli <- coxph(Surv(os_months, os_event) ~., data = mydata2)

save(
  cox_multi_TCGA_gene,
  cox_multi_TCGA_gene_cli,
  cox_multi_TCGA_gene_mean_cli,
  cox_multi_GEO_gene,
  cox_multi_GEO_gene_cli,
  cox_multi_GEO_gene_mean_cli,
  GEO_sample_information,
  TCGA_sample_information,
  file = "cox_models.RData"
)


# 均值分两组 -------------------------------------------------------------------
mydata2$gene_group <- ifelse(
  mydata2$gene_mean > mean(mydata2$gene_mean, na.rm = TRUE),
  "High",
  "Low"
)
TCGA_mydata2$gene_group <- ifelse(
  TCGA_mydata2$gene_mean > mean(TCGA_mydata2$gene_mean, na.rm = TRUE),
  "High",
  "Low"
)

df_prop <- as.data.frame(
  prop.table(
    table(TCGA_mydata2$gene_group, TCGA_mydata2$stage),
    margin = 1
  )
)
colnames(df_prop) <- c("gene_group","stage","Freq")
ggplot(df_prop, aes(x = gene_group, y = Freq, fill = factor(stage))) +
  geom_bar(stat = "identity") +
  ylab("Proportion") +
  xlab("Gene group") +
  labs(fill = "Stage") +
  theme_classic()


df_prop <- as.data.frame(
  prop.table(
    table(mydata2$gene_group, mydata2$stage),
    margin = 1
  )
)
colnames(df_prop) <- c("gene_group","stage","Freq")
ggplot(df_prop, aes(x = gene_group, y = Freq, fill = factor(stage))) +
  geom_bar(stat = "identity") +
  ylab("Proportion") +
  xlab("Gene group") +
  labs(fill = "Stage") +
  theme_classic()

GEO_sample_information<- mydata2
TCGA_sample_information<- TCGA_mydata2
