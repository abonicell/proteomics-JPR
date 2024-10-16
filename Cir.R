library(dplyr)
library(tidyr)
library(ggfortify)
library(devtools)
library(ggbiplot)
library(factoextra)
library(gridExtra)
library(grid)
library(readr)
library(ggpubr)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)

# load data ---------------------------------------------------------------
Cir <- read.csv("/Users/andreabonicelli/Documents/GitHub/proteomics-JPR/data/Cir.csv")

Cir_matrix <- as.matrix(Cir)
Cir_Var <- subset(Cir[16:19])
Cir_matrix <- as.matrix(Cir_Var)

Cir %>%
  
  dplyr::select(Burial, IRSF:Am.P) %>%
  gather(Measure, Value, -Burial) %>%
  ggplot(aes(
    x = factor(Burial),
    y = Value,
    color = Burial,
    fill = Burial
  )) +
  facet_wrap( ~ Measure, scales = "free_y") +
  xlab("Burial condition") + ylab("Ratio") +
  scale_fill_manual(values = alpha(c("#00188f", "#ec008c"), .2)) +
  scale_color_manual(values = c("#00188f", "#ec008c")) +
  theme_bw() + geom_boxplot(width = 0.4, lwd = .5) + rremove("xlab")


# PCA protein -------------------------------------------------------------
df <- read_csv("/Users/andreabonicelli/Documents/GitHub/proteomics-JPR/data/Matrix 75.csv")

do.call(rbind, lapply(Var, function(x)
  shapiro.test(x)[c("statistic", "p.value")]))

lapply(df[, c(5:94)], function(x)
  t.test(x ~ df$Burial, var.equal = TRUE))

#df <- df %>% slice(-c(23,24))

Var <- df[5:94]

cor(Var, df$Age)
cor(Var, df$PMI)

df_names <- data.frame(df)
rownames(Var) <- df_names[, 1]

library(viridis)
library("FactoMineR")

res.pca <- PCA(Var,
               scale.unit = TRUE,
               ncp = 5,
               graph = FALSE)

summary.PCA(res.pca)

dimdesc(res.pca, axes = 1:2, proba = 0.05)

eig.val <- get_eigenvalue(res.pca)
eig.val

fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 30)) + theme_bw() +
  scale_color_viridis(option = "E")

var <- get_pca_var(res.pca)

var

# Coordinates
head(var$coord)
# Cos2: quality on the factor map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)

# Color by cos2 values: quality on the factor map
fviz_pca_var(
  res.pca,
  geom = c("point", "text"),
  gradient.cols = c("#00188f", "#ec008c"),
  col.var = "contrib",
  repel = TRUE,
  legend.title = "Contribution"
) + theme_bw()
# Visualize variable
fviz_contrib(res.pca,
             choice = "var",
             axes = 1:2,
             top = 15) + theme_bw()

fviz_contrib(res.pca,
             choice = "var",
             axes = 2,
             top = 15) + theme_bw()

#First 7 variables are correlated above 0.8
fviz_pca_var(
  res.pca,
  select.var = list(contrib = 10),
  geom = c("point", "text"),
  col.var = "contrib",
  gradient.cols = c("#00188f", "#ec008c"),
  repel = TRUE,
  legend.title = "Contribution"
) + theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18))


fviz_pca_ind(
  res.pca,
  geom.ind = c("point", "text"),
  col.ind = df$Burial,
  palette = c("#00188f", "#ec008c"),
  legend.title = "Burial Conditions",
  addEllipses = FALSE,
  repel = TRUE
) + theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18)) +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  theme(legend.position = "right")

library(cowplot)

plot_grid(
  B,
  A ,
  respect = TRUE,
  labels = c('A', 'B'),
  ncol = 3,
  nrow = 1
)

do.call(rbind, lapply(Var, function(x)
  shapiro.test(x)[c("statistic", "p.value")]))
lapply(df[, c(3:92)], function(x)
  anova(lm(x ~ df$Burial)))

df %>%
  
  dplyr::select(
    Burial,
    APOA1_HUMAN,
    ASPN_HUMAN,
    CLUS_HUMAN,
    COCA1_HUMAN,
    COMA1_HUMAN,
    ENPL_HUMAN,
    EZRI_HUMAN,
    FINC_HUMAN,
    FLNC_HUMAN,
    LUM_HUMAN,
    LYOX_HUMAN,
    MFGM_HUMAN,
    NUCB2_HUMAN,
    OSTCN_HUMAN,
    PGS1_HUMAN,
    POSTN_HUMAN,
    PPBT_HUMAN,
    TENA_HUMAN,
    TETN_HUMAN,
    TTHY_HUMAN
  ) %>%
  gather(Measure, Value, -Burial) %>%
  ggplot(aes(
    x = factor(Burial),
    y = Value,
    color = Burial,
    fill = Burial
  )) + facet_wrap( ~ Measure, scales = "free_y") + theme(legend.position =
                                                           "top") +
  xlab("Deposition") + ylab("Relative abundance") +
  scale_fill_manual(values = alpha(c("#00188f", "#ec008c"), .2)) +
  scale_color_manual(values = c("#00188f", "#ec008c")) +
  theme_bw() + geom_boxplot(width = 0.4, lwd = .5) + rremove("xlab")



# PTMs --------------------------------------------------------------------
#PCA
PTM <- read.csv("/Users/andreabonicelli/Documents/GitHub/proteomics-JPR/data/PTM.csv")

do.call(rbind, lapply(Var_PTM, function(x)
  shapiro.test(x)[c("statistic", "p.value")]))

df_names <- data.frame(PTM)
rownames(PTM) <- df_names[, 1]


#PTM <- PTM %>% slice(-c(17,18))
#PTM <- PTM[ which(PTM$State=='Modified'), ]

Var_PTM <- PTM[3:8]

df_names <- data.frame(PTM)
rownames(Var_PTM) <- df_names[, 1]

res.pca <- PCA(Var_PTM,
               scale.unit = TRUE,
               ncp = 5,
               graph = FALSE)

summary.PCA(res.pca)

dimdesc(res.pca, axes = 1:2, proba = 0.05)

eig.val <- get_eigenvalue(res.pca)
eig.val

fviz_eig(res.pca, addlabels = TRUE) +
  scale_color_viridis(option = "E")

get_pca_var(res.pca)
var

# Coordinates
head(var$coord)
# Cos2: quality on the factor map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)

# Visualize variable
fviz_contrib(res.pca,
             choice = "var",
             axes = 1:2,
             top = 15) + theme_bw()

#First 7 variables are correlated above 0.8
fviz_pca_var(
  res.pca,
  col.var = "contrib",
  geom = c("point", "text"),
  repel = TRUE,
  gradient.cols = c("#00188f", "#ec008c"),
  legend.title = "Contribution"
) + theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18))

#identify the most significantly associated variables with a given principal component
res.desc <- dimdesc(res.pca, axes = c(1, 2), prob = 0.05)
# Description of dimension 1
res.desc$Dim.1
# Description of dimension 2
res.desc$Dim.2

groups <- as.factor(PTM$Burial[1:28])

fviz_pca_ind(
  res.pca,
  geom.ind = c("point", "text"),
  col.ind = PTM$Burial,
  palette = c("#00188f", "#ec008c"),
  legend.title = "Burial Conditions",
  addEllipses = FALSE,
  repel = TRUE
) + theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18)) +
  theme(legend.position = "right")


lapply(PTM[, c(6:17)], function(x)
  anova(lm(x ~ PTM$Burial)))

PTM %>%
  
  dplyr::select(Burial,
                "LUM_FNALQYLR_2_Deamidated_NQ":"TSP1_QHVVSVEEALLATGQWK_1_Deamidated_NQ") %>%
  gather(Measure, Value, -Burial) %>%
  ggplot(aes(
    x = factor(Burial),
    y = Value,
    color = Burial,
    fill = Burial
  )) +
  facet_wrap( ~ Measure, scales = "free_y") +
  theme_bw() + theme(legend.position = "top") +
  xlab("Deposition") + ylab("Modification %") +
  scale_fill_manual(values = alpha(c("#00188f", "#ec008c"), .2)) +
  scale_color_manual(values = c("#00188f", "#ec008c")) +
  theme_bw() + geom_boxplot(width = 0.4, lwd = .5) + rremove("xlab")


# changing with both ------------------------------------------------------
A <- ggerrorplot(
  df,
  x = "Age",
  y = "ALBU_HUMAN",
  ylab = "Relative abundance",
  title = "ALBU_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))
B <- ggerrorplot(
  df,
  x = "PMI",
  y = "ALBU_HUMAN",
  ylab = "Relative abundance",
  title = "ALBU_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

C <- ggerrorplot(
  df,
  x = "Age",
  y = "ASPN_HUMAN",
  ylab = "Relative abundance",
  title = "ASPN_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))
D <- ggerrorplot(
  df,
  x = "PMI",
  y = "ASPN_HUMAN",
  ylab = "Relative abundance",
  title = "ASPN_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

E <- ggerrorplot(
  df,
  x = "Age",
  y = "CLC11_HUMAN",
  ylab = "Relative abundance",
  title = "CLC11_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))
F <- ggerrorplot(
  df,
  x = "PMI",
  y = "CLC11_HUMAN",
  ylab = "Relative abundance",
  title = "CLC11_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

G <- ggerrorplot(
  df,
  x = "Age",
  y = "FETUA_HUMAN",
  ylab = "Relative abundance",
  title = "FETUA_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))
H <- ggerrorplot(
  df,
  x = "PMI",
  y = "FETUA_HUMAN",
  ylab = "Relative abundance",
  title = "FETUA_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

I <- ggerrorplot(
  df,
  x = "Age",
  y = "FMOD_HUMAN",
  ylab = "Relative abundance",
  title = "FMOD_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))
L <- ggerrorplot(
  df,
  x = "PMI",
  y = "FMOD_HUMAN",
  ylab = "Relative abundance",
  title = "FMOD_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

M <- ggerrorplot(
  df,
  x = "Age",
  y = "MIME_HUMAN",
  ylab = "Relative abundance",
  title = "MIME_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))
N <- ggerrorplot(
  df,
  x = "PMI",
  y = "MIME_HUMAN",
  ylab = "Relative abundance",
  title = "MIME_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

O <- ggerrorplot(
  df,
  x = "Age",
  y = "NUCB1_HUMAN",
  ylab = "Relative abundance",
  title = "NUCB1_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))
P <- ggerrorplot(
  df,
  x = "PMI",
  y = "NUCB1_HUMAN",
  ylab = "Relative abundance",
  title = "NUCB1_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

ggarrange(A,
          C,
          E,
          G,
          I,
          M,
          O,
          B,
          D,
          F,
          H,
          L,
          N,
          P,
          ncol = 4,
          nrow = 4)

# changing with PMI only --------------------------------------------------
A <- ggerrorplot(
  df,
  x = "PMI",
  y = "ANT3_HUMAN",
  ylab = "Relative abundance",
  title = "ANT3_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

L <- ggerrorplot(
  df,
  x = "PMI",
  y = "CHAD_HUMAN",
  ylab = "Relative abundance",
  title = "CHAD_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

B <- ggerrorplot(
  df,
  x = "PMI",
  y = "B2MG_HUMAN",
  ylab = "Relative abundance",
  title = "B2MG_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

C <- ggerrorplot(
  df,
  x = "PMI",
  y = "CSPG2_HUMAN",
  ylab = "Relative abundance",
  title = "CSPG2_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

D <- ggerrorplot(
  df,
  x = "PMI",
  y = "G3P_HUMAN",
  ylab = "Relative abundance",
  title = "G3P_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

E <- ggerrorplot(
  df,
  x = "PMI",
  y = "IGL1_HUMAN",
  ylab = "Relative abundance",
  title = "IGL1_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

F <- ggerrorplot(
  df,
  x = "PMI",
  y = "KNG1_HUMAN",
  ylab = "Relative abundance",
  title = "KNG1_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

G <- ggerrorplot(
  df,
  x = "PMI",
  y = "PCOC1_HUMAN",
  ylab = "Relative abundance",
  title = "PCOC1_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

H <- ggerrorplot(
  df,
  x = "PMI",
  y = "PGBM_HUMAN",
  ylab = "Relative abundance",
  title = "PGBM_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

I <- ggerrorplot(
  df,
  x = "PMI",
  y = "RCN3_HUMAN",
  ylab = "Relative abundance",
  title = "RCN3_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

ggarrange(A, L, B, C, D, E, F, G, H, I, ncol = 3, nrow = 4)

# changing with age only --------------------------------------------------
A <- ggerrorplot(
  df,
  x = "Age",
  y = "A4_HUMAN",
  ylab = "Relative abundance",
  title = "A4_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

B <- ggerrorplot(
  df,
  x = "Age",
  y = "CFAB_HUMAN",
  ylab = "Relative abundance",
  title = "CFAB_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

C <- ggerrorplot(
  df,
  x = "Age",
  y = "CO9_HUMAN",
  ylab = "Relative abundance",
  title = "CO9_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

D <- ggerrorplot(
  df,
  x = "Age",
  y = "COBA2_HUMAN",
  ylab = "Relative abundance",
  title = "COBA2_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

E <- ggerrorplot(
  df,
  x = "Age",
  y = "NUCB1_HUMAN",
  ylab = "Relative abundance",
  title = "NUCB1_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

F <- ggerrorplot(
  df,
  x = "Age",
  y = "COCA1_HUMAN",
  ylab = "Relative abundance",
  title = "COCA1_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

G <- ggerrorplot(
  df,
  x = "Age",
  y = "COMA1_HUMAN",
  ylab = "Relative abundance",
  title = "COMA1_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

H <- ggerrorplot(
  df,
  x = "Age",
  y = "EZRI_HUMAN",
  ylab = "Relative abundance",
  title = "EZRI_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

I <- ggerrorplot(
  df,
  x = "Age",
  y = "KAZD1_HUMAN",
  ylab = "Relative abundance",
  title = "KAZD1_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

L <- ggerrorplot(
  df,
  x = "Age",
  y = "MFGM_HUMAN",
  ylab = "Relative abundance",
  title = "MFGM_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

M <- ggerrorplot(
  df,
  x = "Age",
  y = "OMD_HUMAN",
  ylab = "Relative abundance",
  title = "OMD_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

N <- ggerrorplot(
  df,
  x = "Age",
  y = "PDIA1_HUMAN",
  ylab = "Relative abundance",
  title = "PDIA1_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

O <- ggerrorplot(
  df,
  x = "Age",
  y = "PEDF_HUMAN",
  ylab = "Relative abundance",
  title = "PEDF_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

P <- ggerrorplot(
  df,
  x = "Age",
  y = "PROC_HUMAN",
  ylab = "Relative abundance",
  title = "PROC_HUMAN"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

ggarrange(A,
          B,
          C,
          D,
          E,
          F,
          G,
          H,
          I,
          L,
          M,
          N,
          O,
          P,
          ncol = 4,
          nrow = 4)



I <- ggerrorplot(
  PTM,
  x = "Age",
  y = "LUM_FNALQYLR_2_Deamidated_NQ",
  ylab = "Deamidation (%)",
  title = "LUM_FNALQYLR_2_Deamidated_NQ"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

L <- ggerrorplot(
  PTM,
  x = "Age",
  y = "SPRC_MRDWLKNVLVTLYER_7_Deamidated_NQ",
  ylab = "Deamidation (%)",
  title = "SPRC_MRDWLKNVLVTLYER_7_Deamidated_NQ"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

M <- ggerrorplot(
  PTM,
  x = "Age",
  y = "OFL_EIDYIQYLR_6_Deamidated_NQ",
  ylab = "Deamidation (%)",
  title = "OFL_EIDYIQYLR_6_Deamidated_NQ"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

N <- ggerrorplot(
  PTM,
  x = "Age",
  y = "OMD_LQDIPYNIFNLPNIVELSVGHNKLK_10_Deamidated_NQ",
  ylab = "Deamidation (%)",
  title = "OMD_LQDIPYNIFNLPNIVELSVGHNKLK_10_Deamidated_NQ"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

O <- ggerrorplot(
  PTM,
  x = "Age",
  y = "LYOX_VSVNPSYLVPESDYTNNVVR_16_Deamidated_NQ",
  ylab = "Deamidation (%)",
  title = "LYOX_VSVNPSYLVPESDYTNNVVR_16_Deamidated_NQ"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

P <- ggerrorplot(
  PTM,
  x = "Age",
  y = "TSP1_QHVVSVEEALLATGQWK_1_Deamidated_NQ",
  ylab = "Deamidation (%)",
  title = "TSP1_QHVVSVEEALLATGQWK_1_Deamidated_NQ"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

ggarrange(I, L, M, N, O, P, ncol = 3, nrow = 2)

I <- ggerrorplot(
  PTM,
  x = "PMI",
  y = "LUM_FNALQYLR_2_Deamidated_NQ",
  ylab = "Deamidation (%)",
  title = "LUM_FNALQYLR_2_Deamidated_NQ"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

L <- ggerrorplot(
  PTM,
  x = "PMI",
  y = "SPRC_MRDWLKNVLVTLYER_7_Deamidated_NQ",
  ylab = "Deamidation (%)",
  title = "SPRC_MRDWLKNVLVTLYER_7_Deamidated_NQ"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

M <- ggerrorplot(
  PTM,
  x = "PMI",
  y = "OFL_EIDYIQYLR_6_Deamidated_NQ",
  ylab = "Deamidation (%)",
  title = "OFL_EIDYIQYLR_6_Deamidated_NQ"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

N <- ggerrorplot(
  PTM,
  x = "PMI",
  y = "OMD_LQDIPYNIFNLPNIVELSVGHNKLK_10_Deamidated_NQ",
  ylab = "Deamidation (%)",
  title = "OMD_LQDIPYNIFNLPNIVELSVGHNKLK_10_Deamidated_NQ"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

O <- ggerrorplot(
  PTM,
  x = "PMI",
  y = "LYOX_VSVNPSYLVPESDYTNNVVR_16_Deamidated_NQ",
  ylab = "Deamidation (%)",
  title = "LYOX_VSVNPSYLVPESDYTNNVVR_16_Deamidated_NQ"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

P <- ggerrorplot(
  PTM,
  x = "PMI",
  y = "TSP1_QHVVSVEEALLATGQWK_1_Deamidated_NQ",
  ylab = "Deamidation (%)",
  title = "TSP1_QHVVSVEEALLATGQWK_1_Deamidated_NQ"
) +
  theme_bw() + theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 14))

ggarrange(I, L, M, N, O, P, ncol = 3, nrow = 2)


ggscatter(
  df,
  x = "Age",
  y = "COCA1_HUMAN",
  ylab = "Relative abundance",
  title = "COCA1_HUMAN",
  add = "reg.line",
  shape = "Burial",
  color = "PMI",
  cor.coef = TRUE,
  add.params = list(color = "black", linetype = "dashed"),
  cor.method = "spearman"
) + gradient_color(c("#00188f", "#ec008c")) +  theme_bw() +
  ggforce::facet_zoom(xlim = c(60, 90), zoom.size = 0.7)

df_zinc <- subset(df, df$Burial == "Zinc coffin")
Z1 <- ggscatter(
  df_zinc,
  x = "Age",
  y = "COCA1_HUMAN",
  ylab = "Relative abundance",
  title = "COCA1_HUMAN, zinc coffin",
  add = "reg.line",
  color = "PMI",
  cor.coef = TRUE,
  add.params = list(color = "black", linetype = "dashed"),
  cor.method = "spearman"
) + gradient_color("RdYlBu") + theme_calc()


df_wood <- subset(df, df$Burial == "Wood coffin")
Z2 <- ggscatter(
  df_wood,
  x = "Age",
  y = "COCA1_HUMAN",
  ylab = "Relative abundance",
  title = "COCA1_HUMAN, wood coffin",
  add = "reg.line",
  color = "PMI",
  cor.coef = TRUE,
  add.params = list(color = "black", linetype = "dashed"),
  cor.method = "spearman"
) + gradient_color("RdYlBu") + theme_bw() +
  ggforce::facet_zoom(xlim = c(60, 90), zoom.size = 1)

ggarrange(Z1,
          Z2,
          ncol = 1,
          nrow = 2,
          common.legend = TRUE)

