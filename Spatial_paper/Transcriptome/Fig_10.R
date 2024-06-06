library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggpubr)

df <- data.table( read.csv("../../df_060322_Annotated.csv") )

Acutus_Chl <- read.csv("Acutus_Chlcorr.csv")
Prop_Chl <- read.csv("Propinquus_Chlcorr.csv")

names(Acutus_Chl)[1] <- "Gene"
names(Prop_Chl)[1] <- "Gene"

Acutus_Chl$Gene <- gsub("Cluster-", "Acutus_", Acutus_Chl$Gene)
Prop_Chl$Gene <- gsub("Cluster-", "Prop_", Prop_Chl$Gene)

chl_df <- rbind(Acutus_Chl, Prop_Chl) %>% dplyr::select(Gene, logFC, t, adj.P.Val)
names(chl_df)[2:4] <- c("Chl_LFC", "Chl_t", "Chl_p")

df <- join(df, chl_df, type="left", by="Gene")

df <- df[!is.na(Chl_LFC) & !is.na(logFC),]

png(res=300, file="F:/Antarctic_copepods/Chl/Cor_fig.png", units='in', height=4, width=6)
ggplot(df[DE=="DE"], aes(x=logFC, y=Chl_LFC, color=Species)) + geom_point(alpha = 0.1) +
    geom_smooth(method='lm', se=F) + stat_cor(aes(label = ..r.label..), method="spearman", cor.coef.name="rho") +
    theme_bw() + facet_grid(.~Species) + labs(y="Chl response", x="Starvation response") + theme(legend.position = "none")
dev.off()


