#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(edgeR))
options(stringsAsFactors=FALSE)
theme_set(theme_classic())

# helper variables
cols.times <- c(Baseline="white", "4 DPI"="#fecc5c", 
                "7 DPI"="#fd8d3c", Terminal="#e31a1c", 
                "10 DPI"="black", ">10 DPI"="grey")
cols.reg <- c("HeV-R"="#e31a1c", "No difference"="grey", "HeV-g2"="#377eb8")

# helper functions
run.de.strain <- function(timepoint, counts.mat, meta.mat, psig=0.05, fsig=1) {
  # subset by time-point and re-level agent column
  meta.mat <- meta.mat %>%
             filter(Timepoint==timepoint) %>%
             mutate(Agent=factor(Agent, 
                                 levels=c("HeV-g2", "HeV-R"),
                                 labels=c("g2", "Redlands")))
  counts.mat <- counts.mat[ , rownames(meta.mat)]
  
  # define model by "agent" column and define contrasts
  mm <- model.matrix(~0 + meta.mat$Agent)
  colnames(mm) <- levels(meta.mat$Agent)
  cont <- makeContrasts(Redlands-g2, levels=mm)
  
  # prep counts. Remove genes where <3 samples are above threshold (20),
  # normalize, and fit a linear model
  fit <- counts.mat[rowSums(counts.mat > 20) > 3, ] %>% 
         DGEList() %>%
         calcNormFactors() %>% 
         voom(design=mm) %>% # normalize
         lmFit(design=mm)
  
  # run differential expression, annotate, and return data frame
  degenes <- contrasts.fit(fit, cont) %>%
             eBayes() %>%
             topTable(n=Inf) %>%
             rownames_to_column("Gene") %>%
             mutate(psig=(adj.P.Val < psig),
                    fsig=(abs(logFC) > fsig),
                    padj=adj.P.Val,
                    Significant=(psig & fsig),
                    Regulation="No difference",
                    Timepoint=timepoint) %>%
             select(Gene, logFC, padj, Significant, Regulation, Timepoint)
  degenes[degenes$Significant & degenes$logFC > 0, "Regulation"] <- "HeV-R"
  degenes[degenes$Significant & degenes$logFC < 0, "Regulation"] <- "HeV-g2"
  degenes$Regulation <- factor(degenes$Regulation, 
                               levels=c("HeV-R", "No difference", "HeV-g2"))
  return(degenes)
}

run.de.timepoint <- function(counts.mat, meta.mat, model.mat, contrast, 
                             psig=0.05, fsig=1) {
  # prep counts. Remove genes where <3 samples are above threshold (20),
  # normalize, and fit a linear model
  fit <- counts.mat[rowSums(counts.mat > 20) > 3, ] %>% 
         DGEList() %>%
         calcNormFactors() %>% 
         voom(design=model.mat) %>% # normalize
         lmFit(design=model.mat)
  
  # run differential expression, annotate, and return data frame
  degenes <- contrasts.fit(fit, contrast) %>%
             eBayes() %>%
             topTable(n=Inf) %>%
             rownames_to_column("Gene") %>%
             mutate(psig=(adj.P.Val < psig),
                    fsig=(abs(logFC) > fsig),
                    padj=adj.P.Val,
                    Significant=(psig & fsig),
                    Regulation="No difference") %>%
             select(Gene, logFC, padj, Significant, Regulation)
  degenes$Regulation <- "None"
  degenes[degenes$Significant & degenes$logFC > 0, "Regulation"] <- "Up"
  degenes[degenes$Significant & degenes$logFC < 0, "Regulation"] <- "Down"
  degenes$Regulation <- factor(degenes$Regulation, 
                               levels=c("Up", "None", "Down"))
  return(degenes)
}

## inputs ----------------------------------------------------------------------
# load metadata and plot samples
meta <- read.csv("samplesheet.csv", na.strings="") %>%
        filter(QC.pass) %>%
        mutate(Timepoint=factor(Timepoint, levels=names(cols.times)))
rownames(meta) <- str_replace(meta$ID, "-", ".")
meta %>%
  ggplot(aes(DPI, NHP, group=NHP, fill=Timepoint)) +
  geom_line() +
  geom_point(pch=21, size=3) +
  scale_fill_manual(values=cols.times) +
  facet_wrap(~Agent, ncol=2, scales="free_y") +
  labs(x="Days postinfection",
       y=element_blank())
ggsave("analysis/samples.png",
       units="cm", width=15, height=8)

# thresholded counts
counts <- read.csv("analysis/counts-thresholded.csv", 
                   row.names=1) %>%
          as.matrix()
# remove controls
counts <- counts[!str_detect(rownames(counts), "^NEG|^POS"), ]

# align rows and columns
x <- intersect(rownames(meta), colnames(counts))
counts <- counts[ , x]
meta <- meta[x, ]
rm(x)

## QC: PCA ---------------------------------------------------------------------
# all samples
# filter genes with no counts above the threshold
dds <- counts[rowSums(counts > 20) > 3, ] %>%
       DGEList() %>%
       calcNormFactors() %>% 
       cpm()
# save CPM for CIBERSORT
dds %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  write.table("analysis/cibersort-input.tsv", sep="\t")
# log transform and format
dds <- dds %>%
       log2() %>% 
       t()

# run PCA
pca <- dds %>%
        prcomp()
pcs <- summary(pca)$importance["Proportion of Variance", 1:2]
pcs <- round(100*pcs)
pcs <- paste0(c("PC1 (", "PC2 ("), pcs, "%)")
pca <- pca$x %>%
       as.data.frame() %>%
       rownames_to_column("ID") %>%
       mutate(ID=str_replace(ID, "\\.", "-")) %>%
       select(ID, PC1, PC2) %>%
       left_join(meta, by="ID")
pca %>%
  ggplot(aes(PC1, PC2, shape=Agent, fill=Timepoint)) +
  geom_hline(yintercept=0, linetype=2, col="grey") +
  geom_vline(xintercept=0, linetype=2, col="grey") +
  geom_point(size=5) +
  scale_fill_manual(values=cols.times) +
  scale_shape_manual(values=c("HeV-R"=21, "HeV-g2"=22)) +
  guides(fill=guide_legend(override.aes=list(shape=21))) +
  labs(x=pcs[1],
       y=pcs[2])
ggsave("analysis/pca-all.png",
       units="cm", width=15, height=8)
rm(pca, pcs)

# clean up
rm(dds)

## differential expression analyses: by viral strain ---------------------------
# put 10 DPI and terminal on the same level
m <- meta %>%
     mutate(Timepoint=factor(Timepoint,
                             levels=c("Baseline", "4 DPI", "7 DPI",
                                      "Terminal", "10 DPI"),
                             labels=c("Baseline", "4 DPI", "7 DPI",
                                      "Terminal/10 DPI", "Terminal/10 DPI")))

# run DE between HeV-g2 and HeV-R and save outputs
# do NOT run on 7 DPI because there are only 2x HeV-R 
dds <- c("Baseline", "4 DPI", "Terminal/10 DPI") %>%
       lapply(run.de.strain, counts.mat=counts, meta.mat=m) %>%
       do.call(rbind, .) %>%
       mutate(Timepoint=factor(Timepoint, 
                               levels=c("Baseline", "4 DPI", 
                                        "Terminal/10 DPI")))
write.csv(dds, "analysis/results-strain.csv", row.names=FALSE)

# total DE genes per timepoint
dds %>%
  filter(Significant) %>%
  group_by(Timepoint, Regulation) %>%
  summarise(Genes=n(),
            .groups="drop") %>% 
  # add in missing zeros
  right_join(expand.grid(Timepoint=levels(dds$Timepoint),
                         Regulation=c("HeV-R", "HeV-g2")),
             by=c("Timepoint", "Regulation")) %>%
  replace_na(list(Genes=0)) %>%
  # update group labels
  mutate(Regulation=factor(Regulation, 
                           labels=c("Higher in HeV-R", 
                                    "Higher in HeV-g2"))) %>%
  ggplot(aes(Timepoint, Genes)) +
  geom_col(col="black", fill="black") +
  geom_text(aes(label=Genes), nudge_y=5) +
  scale_y_continuous(limits=c(0, 100)) +
  facet_wrap(~Regulation) +
  labs(x=element_blank(),
       y="Total significantly DE genes")
ggsave("analysis/degenes-strain.png",
       units="cm", width=15, height=10)

# volcano for terminal/10 DPI
top.de <- dds %>%
          filter(Timepoint=="Terminal/10 DPI",
                 Significant) %>%
          group_by(Regulation) %>%
          top_n(n=12, wt=abs(logFC)) %>% 
          ungroup()
dds %>%
  filter(Timepoint=="Terminal/10 DPI") %>%
  ggplot(aes(logFC, -log10(padj), col=Regulation, size=Regulation)) +
  geom_point() +
  ggrepel::geom_text_repel(data=top.de, aes(label=Gene), 
                           col="black", size=3, force=25, max.iter=1e6) +
  xlim(-4, 4) +
  scale_color_manual(values=cols.reg) +
  scale_size_manual(values=c(1, 0.5, 1)) +
  labs(x="Fold change (log2)",
       y="FDR-adjusted p-value",
       col="Higher in",
       size="Higher in")
ggsave("analysis/volcano-strain.png",
       units="cm", width=15, height=10) 

# clean up
rm(m, dds, top.de)

## differential expression analyses: 4 DPI vs. baseline ------------------------
# g2
m <- meta %>%
     filter(Agent=="HeV-g2",
            Timepoint %in% c("Baseline", "4 DPI")) %>%
     mutate(Timepoint=droplevels(Timepoint))
c <- counts[, rownames(m)]
# define model 
mm <- model.matrix(~0 + m$Timepoint)
colnames(mm) <- c("Baseline", "DPI4")
dds <- run.de.timepoint(c, m, mm, makeContrasts(DPI4-Baseline, 
                                              levels=colnames(mm))) %>%
       mutate(Agent="HeV-g2")

# Redlands
m <- meta %>%
     filter(Agent=="HeV-R",
            Timepoint %in% c("Baseline", "4 DPI")) %>%
     mutate(Timepoint=droplevels(Timepoint))
c <- counts[, rownames(m)]
# define model 
mm <- model.matrix(~0 + m$Timepoint)
colnames(mm) <- c("Baseline", "DPI4")
dds <- run.de.timepoint(c, m, mm, makeContrasts(DPI4-Baseline, 
                                                levels=colnames(mm))) %>%
       mutate(Agent="HeV-R") %>%
       # add g2 results
       rbind(dds)

# Venn diagram of DE genes
dds <- list("HeV-R"=dds %>%
                    filter(Agent=="HeV-R", 
                           Significant) %>%
                    select(Gene) %>%
                    unlist(),
            "HeV-g2"=dds %>%
                     filter(Agent=="HeV-g2", 
                            Significant) %>%
                     select(Gene) %>%
                     unlist())
VennDiagram::venn.diagram(dds, "analysis/venn-4dpi.png", imagetype="png", 
                          units="in", width=2.5, height=2.5,
                          disable.logging=TRUE,  
                          fontfamily="sans", cat.fontfamily="sans", 
                          margin=0.05, cat.pos=c(-140, 140), cat.dist=0.06, 
                          fill=cols.reg[c("HeV-R", "HeV-g2")], alpha=0.5)

# clean up
rm(c, dds, m, mm)

## differential expression analyses: 10/term vs. baseline ----------------------
# g2
m <- meta %>%
     filter(Agent=="HeV-g2",
            Timepoint %in% c("Baseline", "10 DPI")) %>%
     mutate(Timepoint=droplevels(Timepoint))
c <- counts[, rownames(m)]
# define model 
mm <- model.matrix(~0 + m$Timepoint)
colnames(mm) <- c("Baseline", "DPI10")
dds <- run.de.timepoint(c, m, mm, makeContrasts(DPI10-Baseline, 
                                                levels=colnames(mm))) %>%
       mutate(Agent="HeV-g2")

# Redlands
m <- meta %>%
     filter(Agent=="HeV-R",
            Timepoint %in% c("Baseline", "Terminal")) %>%
     mutate(Timepoint=droplevels(Timepoint))
c <- counts[, rownames(m)]
# define model 
mm <- model.matrix(~0 + m$Timepoint)
colnames(mm) <- c("Baseline", "Terminal")
dds <- run.de.timepoint(c, m, mm, makeContrasts(Terminal-Baseline, 
                                                levels=colnames(mm))) %>%
       mutate(Agent="HeV-R") %>%
       # add g2 results
       rbind(dds)

# Venn diagram of DE genes
dds <- list("HeV-R"=dds %>%
              filter(Agent=="HeV-R", 
                     Significant) %>%
              select(Gene) %>%
              unlist(),
            "HeV-g2"=dds %>%
              filter(Agent=="HeV-g2", 
                     Significant) %>%
              select(Gene) %>%
              unlist())
VennDiagram::venn.diagram(dds, "analysis/venn-10dpiterm.png", imagetype="png", 
                          units="in", width=2.5, height=2.5,
                          disable.logging=TRUE,  
                          fontfamily="sans", cat.fontfamily="sans", 
                          margin=0.05, cat.pos=c(-40, 30), cat.dist=0.06, 
                          fill=cols.reg[c("HeV-R", "HeV-g2")], alpha=0.5)

# clean up
rm(c, dds, m, mm)

## differential expression analyses: >10 DPI vs. baseline ----------------------
# filter metadata to HeV-g2 only, baseline vs. >10 DPI
m <- meta %>%
     filter(Agent=="HeV-g2",
            Timepoint %in% c(">10 DPI", "Baseline")) %>%
     mutate(Timepoint=factor(Timepoint, levels=unique(Timepoint)))
# filter counts
c <- counts[ , rownames(m)]
mm <- model.matrix(~0 + m$Timepoint)
colnames(mm) <- c("Baseline", "Convalescent")
dds <- run.de.timepoint(c, m, mm, makeContrasts(Convalescent-Baseline, 
                                                levels=colnames(mm)))
write.csv(dds, "analysis/results-recovery.csv", row.names=FALSE)

# total DE genes
dds %>%
  filter(Significant) %>%
  group_by(Regulation) %>%
  summarise(Genes=n(),
            .groups="drop") %>% 
  # add in missing zeros
  right_join(data.frame(Regulation=c(">10 DPI", "Baseline")),
             by="Regulation") %>%
  replace_na(list(Genes=0)) %>%
  ggplot(aes(Regulation, Genes)) +
  geom_col(col="black", fill="black") +
  geom_text(aes(label=Genes), nudge_y=0.5) +
  scale_y_continuous(limits=c(0, 5), breaks=0:5) +
  labs(x="Higher expression",
       y="Significantly DE genes")
ggsave("analysis/degenes-recovery.png",
       units="cm", width=10, height=10)

# volcano plot
dds %>%
  ggplot(aes(logFC, -log10(padj), col=Regulation)) +
  geom_point() +
  scale_y_continuous(limits=c(0, 5), breaks=0:5) +
  xlim(-5, 5) +
  scale_color_manual(values=cols.reg) +
  labs(x="Fold change (log2)",
       y="FDR-adjusted p-value",
       col="Higher in") +
  theme(legend.position="bottom")
ggsave("analysis/volcano-recovery.png",
       units="cm", width=10, height=10) 

# clean up
rm(m, dds, c, mm)

## DCQ -------------------------------------------------------------------------
# read cibersort output. Remove columns with all zeros plus rows with p > 0.05
dcq <- read.csv("analysis/cibersort-output.csv") %>%
       filter(P.value < 0.05) %>% 
       select(-P.value, -Correlation, -RMSE) %>%
       column_to_rownames("Mixture")

# remove cell types with no samples >10%
dcq <- dcq[, apply(dcq, 2, max) > 0.1]
dcq$ID <- str_replace(rownames(dcq), "\\.", "-")

# move to "long" format and add in metadata
dcq <- dcq %>%
       reshape2::melt(id.vars="ID",
                      variable.name="Celltype",
                      value.name="Abundance") %>%
       left_join(meta, by="ID") %>%
       # re-format celltype names
       # update time point levels to group 10 DPI and terminal together
       mutate(Celltype=factor(Celltype, levels=levels(Celltype),
                              labels=c("Memory B-cells", 
                                       "CD8 T-cells",
                                       "Naive CD4 T-cells",
                                       "Memory CD4 T-cells",
                                       "γδ T-cells", 
                                       "NK cells",
                                       "Monocytes", 
                                       "Macrophages",
                                       "Mast cells", 
                                       "Neutrophils")),
              Timepoint=factor(Timepoint, labels=c("Baseline", "4 DPI", "7 DPI",
                                                   "Terminal/10 DPI",
                                                   "Terminal/10 DPI", ">10 DPI"))) 

# remove the terminal sample from g2
dcq <- dcq[dcq$ID != "nhpD7239-dpi10", ]

# wilcoxon tests for each comparison: no significance after 
dcq %>%
  # remove >10 DPI and 7 DPI
  filter(Timepoint != ">10 DPI",
         Timepoint != "7 DPI") %>%
  select(Celltype, Timepoint, Agent, Abundance) %>%
  group_by(Celltype, Timepoint) %>%
  rstatix::wilcox_test(Abundance ~ Agent) %>%
  # adjust P-value 
  rstatix::adjust_pvalue(method="bonferroni") %>%
  arrange(p.adj)

dcq %>%
  ggplot(aes(Timepoint, Abundance)) +
  geom_boxplot(aes(fill=Agent), outliers=FALSE, alpha=0.5) +
  geom_point(aes(fill=Agent, group=Agent), pch=21, 
             position=position_dodge(width=0.75)) +
  scale_fill_manual(values=cols.reg) +
  facet_wrap(~Celltype, ncol=5) +
  labs(y="Relative abundance") +
  theme(legend.position=c(0.1, 0.9), 
        axis.text.x=element_text(angle=45, hjust=1)) 
ggsave("analysis/cibersort-output.png",
       units="cm", width=20, height=15)

## done! -----------------------------------------------------------------------
sessionInfo()
