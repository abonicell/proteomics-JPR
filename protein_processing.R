# Load required libraries
library(dplyr)
library(tidyr)
library(factoextra)
library(readr)
library(ggpubr)
library(ggplot2)
library(viridis)
library(FactoMineR)
library(cowplot)

# --- Load data --------------------------------------------------------------
Cir <- read.csv("/Users/andreabonicelli/Documents/GitHub/proteomics-JPR/data/Cir.csv")

# Prepare Cir data matrix (subset columns 16 to 19)
Cir_Var <- Cir[, 16:19]
Cir_matrix <- as.matrix(Cir_Var)

# Plot Cir data boxplots by Burial
Cir %>%
  dplyr::select(Burial, IRSF:Am.P) %>%
  gather(Measure, Value, -Burial) %>%
  ggplot(aes(x = factor(Burial), y = Value, color = Burial, fill = Burial)) +
  facet_wrap(~Measure, scales = "free_y") +
  xlab("Burial condition") + ylab("Ratio") +
  scale_fill_manual(values = alpha(c("#00188f", "#ec008c"), .2)) +
  scale_color_manual(values = c("#00188f", "#ec008c")) +
  theme_bw() + geom_boxplot(width = 0.4, lwd = .5) +
  theme(axis.title.x = element_blank()) +
  stat_compare_means()

# --- PCA protein -------------------------------------------------------------
df <- read_csv("/Users/andreabonicelli/Documents/GitHub/proteomics-JPR/data/Matrix 75.csv")

Var <- df[, 5:94]

# Shapiro tests for normality
do.call(rbind, lapply(Var, function(x) shapiro.test(x)[c("statistic", "p.value")]))

# T-tests by Burial for all proteins
lapply(df[, 5:94], function(x) t.test(x ~ df$Burial, var.equal = TRUE))

# Correlations with Age and PMI
cor(Var, df$Age)
cor(Var, df$PMI)

# Assign rownames for PCA variables
rownames(Var) <- df[[1]]

# PCA analysis
res.pca <- PCA(Var, scale.unit = TRUE, ncp = 5, graph = FALSE)
summary(res.pca)
dimdesc(res.pca, axes = 1:2, proba = 0.05)

# Eigenvalues and scree plot
eig.val <- get_eigenvalue(res.pca)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 30)) +
  theme_bw() + scale_color_viridis(option = "E")

# PCA variable contributions and coordinates
var <- get_pca_var(res.pca)
head(var$coord)
head(var$cos2)
head(var$contrib)

# Plot PCA variables colored by contribution
fviz_pca_var(
  res.pca,
  geom = c("point", "text"),
  gradient.cols = c("#00188f", "#ec008c"),
  col.var = "contrib",
  repel = TRUE,
  legend.title = "Contribution"
) + theme_bw()

# Top contributing variables to PCs
fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 15) + theme_bw()
fviz_contrib(res.pca, choice = "var", axes = 2, top = 15) + theme_bw()

# Select top 10 variables by contribution
fviz_pca_var(
  res.pca,
  select.var = list(contrib = 10),
  geom = c("point", "text"),
  col.var = "contrib",
  gradient.cols = c("#00188f", "#ec008c"),
  repel = TRUE,
  legend.title = "Contribution"
) + theme_bw() +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18))

# PCA individuals colored by Burial
fviz_pca_ind(
  res.pca,
  geom.ind = c("point", "text"),
  col.ind = df$Burial,
  palette = c("#00188f", "#ec008c"),
  legend.title = "Burial Conditions",
  addEllipses = FALSE,
  repel = TRUE
) + theme_bw() +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18)) +
  theme(legend.position = "right")

# Combine plots example (replace B and A with your plots)
# plot_grid(B, A, labels = c('A', 'B'), ncol = 2)

# ANOVA tests for proteins by Burial
do.call(rbind, lapply(Var, function(x) shapiro.test(x)[c("statistic", "p.value")]))
lapply(df[, 3:92], function(x) anova(lm(x ~ df$Burial)))

# Boxplot of selected proteins by Burial
df %>%
  dplyr::select(
    Burial,
    APOA1_HUMAN, ASPN_HUMAN, CLUS_HUMAN, COCA1_HUMAN, COMA1_HUMAN,
    ENPL_HUMAN, EZRI_HUMAN, FINC_HUMAN, FLNC_HUMAN, LUM_HUMAN,
    LYOX_HUMAN, MFGM_HUMAN, NUCB2_HUMAN, OSTCN_HUMAN, PGS1_HUMAN,
    POSTN_HUMAN, PPBT_HUMAN, TENA_HUMAN, TETN_HUMAN, TTHY_HUMAN
  ) %>%
  gather(Measure, Value, -Burial) %>%
  ggplot(aes(x = factor(Burial), y = Value, color = Burial, fill = Burial)) +
  facet_wrap(~Measure, scales = "free_y") +
  theme(legend.position = "top") +
  xlab("Deposition") + ylab("Relative abundance") +
  scale_fill_manual(values = alpha(c("#00188f", "#ec008c"), .2)) +
  scale_color_manual(values = c("#00188f", "#ec008c")) +
  theme_bw() + geom_boxplot(width = 0.4, lwd = .5) +
  theme(axis.title.x = element_blank()) +
  stat_compare_means()

# --- PTMs PCA ---------------------------------------------------------------
PTM <- read.csv("/Users/andreabonicelli/Documents/GitHub/proteomics-JPR/data/PTM.csv")

Var_PTM <- PTM[, 3:8]
rownames(Var_PTM) <- PTM[[1]]

# PCA on PTM data
res.pca.ptm <- PCA(Var_PTM, scale.unit = TRUE, ncp = 5, graph = FALSE)
summary(res.pca.ptm)
dimdesc(res.pca.ptm, axes = 1:2, proba = 0.05)
eig.val.ptm <- get_eigenvalue(res.pca.ptm)

fviz_eig(res.pca.ptm, addlabels = TRUE) + scale_color_viridis(option = "E")

var.ptm <- get_pca_var(res.pca.ptm)

# PCA contributions plot
fviz_contrib(res.pca.ptm, choice = "var", axes = 1:2, top = 15) + theme_bw()
fviz_pca_var(
  res.pca.ptm,
  col.var = "contrib",
  geom = c("point", "text"),
  repel = TRUE,
  gradient.cols = c("#00188f", "#ec008c"),
  legend.title = "Contribution"
) + theme_bw() +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18))

# Description of PCA dimensions
res.desc <- dimdesc(res.pca.ptm, axes = c(1, 2), prob = 0.05)
print(res.desc$Dim.1)
print(res.desc$Dim.2)

# PCA individuals colored by Burial for PTM data
fviz_pca_ind(
  res.pca.ptm,
  geom.ind = c("point", "text"),
  col.ind = PTM$Burial,
  palette = c("#00188f", "#ec008c"),
  legend.title = "Burial Conditions",
  addEllipses = FALSE,
  repel = TRUE
) + theme_bw() +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18)) +
  theme(legend.position = "right")

# ANOVA on PTM variables by Burial
lapply(PTM[, 6:17], function(x) anova(lm(x ~ PTM$Burial)))

# PTM Boxplots by Burial
PTM %>%
  dplyr::select(Burial, LUM_FNALQYLR_2_Deamidated_NQ:TSP1_QHVVSVEEALLATGQWK_1_Deamidated_NQ) %>%
  gather(Measure, Value, -Burial) %>%
  ggplot(aes(x = factor(Burial), y = Value, color = Burial, fill = Burial)) +
  facet_wrap(~Measure, scales = "free_y") +
  theme(legend.position = "top") +
  xlab("Deposition") + ylab("Relative abundance") +
  scale_fill_manual(values = alpha(c("#00188f", "#ec008c"), .2)) +
  scale_color_manual(values = c("#00188f", "#ec008c")) +
  theme_bw() + geom_boxplot(width = 0.4, lwd = .5) +
  theme(axis.title.x = element_blank()) +
  stat_compare_means()

# --- Error plots for Age and PMI --------------------------------------------

# Example data for Age and PMI (replace with your real error data)
error_Age <- data.frame(
  Age = 1:10,
  Error = runif(10, 0.1, 0.3)
)

error_PMI <- data.frame(
  PMI = 1:10,
  Error = runif(10, 0.1, 0.3)
)

# Plot Error vs Age
p1 <- ggplot(error_Age, aes(x = Age, y = Error)) +
  geom_point(color = "#00188f") +
  geom_smooth(method = "lm", color = "#ec008c") +
  theme_bw() +
  ggtitle("Error vs Age")

# Plot Error vs PMI
p2 <- ggplot(error_PMI, aes(x = PMI, y = Error)) +
  geom_point(color = "#00188f") +
  geom_smooth(method = "lm", color = "#ec008c") +
  theme_bw() +
  ggtitle("Error vs PMI")

# Combine error plots side-by-side
plot_grid(p1, p2, labels = c("A", "B"), ncol = 2)
