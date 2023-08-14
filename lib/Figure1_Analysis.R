
# Add to any ggplot2 plot to make it look better
improve_plot <- function() {

  library('ggplot2')
  t <- theme(axis.text = element_text(color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(linetype="dashed"))
  return(t)
}


# This function `run_pcoa` performs principal coordinate analysis (PCoA) 
# on a given phylogenetic tree and abundance matrix, and then annotates 
# the results using metadata.
#
# Args:
#   tree_str: A phylogenetic tree, usually in Newick format or as a tree object.
#   df_metaphlan: A dataframe of microbial abundances. Rows represent microbes 
#                 and columns represent samples.
#   df_meta: A dataframe of metadata for the samples. It should at least contain 
#            a 'sample' column to match with the columns of `df_metaphlan`.
#

run_pcoa <- function(tree_str, df_metaphlan, df_meta) {

    library('ape')
    library('phyloseq')
    library('tidyverse')
    library('vegan')

    tree <- read.tree(text=tree_str)
    if (str_detect(tree$tip.label[1],fixed("|")))
        tree$tip.label <- sapply(tree$tip.label, function(x) str_split(x,fixed("|"))[[1]][7])

    otu <- otu_table(as.matrix(df_metaphlan), TRUE)
    phylo <- phyloseq(otu, tree)
    dists <- UniFrac(phylo)
    pcoa.results <- pcoa(dists)
    # Run significance testing
    permanova <- adonis2(formula = t(dists) ~ Detected, data = df_meta)
    print(permanova)
    return(pcoa.results)
}


# This function `plot_pcoa` visualizes the results of a principal coordinate analysis (PCoA) 
# using microbial abundance data and accompanying metadata.
#
# Args:
#   pcoa.results: An object containing the PCoA results, typically produced by a function like `pcoa`.
#   df_metaphlan: A dataframe of microbial abundances. Rows represent microbes 
#                 and columns represent samples.
#   df_meta: A dataframe of metadata for the samples. It should at least contain 
#            a 'sample' column to match with the columns of `df_metaphlan`.
#   outfile_plot: A string indicating the path and name of the file where the plot 
#                 should be saved.
#   add_loading: A logical value. If TRUE, loadings (correlations of taxa to the PCoA axes) 
#                will be added to the plot. Defaults to FALSE.
#   taxa_oi: A vector of taxonomic names of interest for loading vectors.

plot_pcoa <- function(pcoa.results, df_metaphlan, df_meta, outfile_plot, add_loading=FALSE,taxa_oi=NULL) {
    library('extrafont')
    df_plot <- pcoa.results$vectors %>% data.frame %>% select(Axis.1, Axis.2)
    df_plot$sample <- str_replace_all(rownames(df_plot), "^X","")
    df_plot <- left_join(df_plot, df_meta)

    p <- ggplot(df_plot, aes(x = Axis.1, y = Axis.2)) + geom_point(aes(fill=Detected), pch=21, size=3, alpha=0.65) + scale_fill_manual(labels = c("GA3-", "GA3+"), values = c("steelblue", "#c33b3b")) + theme_bw(16) + theme(axis.text=element_text(color="black"),
        panel.grid.major.x = element_line(linetype="dashed"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(linetype="dashed")) + xlab("PCoA1") + ylab("PCoA2")

    # Add loading vectors if requested
    if (add_loading) {
        U <- calculate_loading_vector(pcoa.results, df_metaphlan, taxa_oi)
        p <- p + geom_segment(size=0.5, color="black", aes(x = 0, y = 0, xend = U[1], yend = U[2]), arrow = arrow(length = unit(0.25, "cm")), alpha=1.0) + annotate("text",x = U[1]+0.065, y = U[2]+0.01, label=taxa_oi)
    }

    pdf(outfile_plot, width=4.5, height=3.25, family="Arial")
    print(p)
    dev.off()

}

# This function `calculate_loading_vector` computes the loadings (correlations of taxa to the PCoA axes)
# based on the methodology used in `biplot.pcoa`.
#
# Source: Inspired by the `biplot.pcoa` function.
#
# Args:
#   pcoa.results: An object containing the PCoA results, typically produced by a function like `pcoa`.
#   df_metaphlan: A dataframe of microbial abundances. Rows represent microbes 
#                 and columns represent samples.
#   taxa_oi: A vector of taxonomic names of interest. This function will calculate 
#            loading vectors specifically for these taxa.

calculate_loading_vector <- function(pcoa.results, df_metaphlan, taxa_oi) {

    x <- pcoa.results
    Y <- t(df_metaphlan)
    plot.axes <- c(1,2)

    pr.coo <- x$vectors
    n <- nrow(Y)
    points.stand <- scale(pr.coo[, plot.axes])
    S <- cov(Y, points.stand)
    U <- S %*% diag((x$values$Eigenvalues[plot.axes]/(n - 
        1))^(-0.5))
    colnames(U) <- colnames(pr.coo[, plot.axes])
    return(U[taxa_oi,])

}

# This function `prep_tree_for_philr` prepares a given phylogenetic tree and species abundance
# dataframe for phylogenetic ILR (Isometric Log Ratio) transformation using the `philr` method.
#
# Args:
#   tree_str: A character string or object representing the phylogenetic tree, often in Newick format.
#   df_species: A dataframe where rows represent different samples and columns represent species 
#               abundances. Each entry gives the abundance of a species in a particular sample.

prep_tree_for_philr <- function(tree_str, df_species) {
    library('phyloseq')
    library('ape')
    library('stringr')

    tree <- read.tree(text=tree_str)
    if (str_detect(tree$tip.label[1],fixed("|")))
        tree$tip.label <- sapply(tree$tip.label, function(x) str_split(x,fixed("|"))[[1]][7])

    tree_sub <- keep.tip(phyloseq(tree), tree$tip[tree$tip%in%rownames(df_species)])
    tree_sub <- makeNodeLabel(tree_sub)
    return(tree_sub)

}

# This function `run_philr` performs a phylogenetic Isometric Log-Ratio (ILR) transformation on a given
# species abundance dataframe using a provided sub-tree. The `philr` method is used for the transformation.
#
# Args:
#   df_species: A dataframe where rows represent different samples and columns represent species 
#               abundances. Each entry gives the abundance of a species in a particular sample.
#   tree_sub: A phylogenetic sub-tree object, typically derived from the main phylogenetic tree,
#             which is used for the ILR transformation process.

run_philr <- function(df_species, tree_sub) {
    library('philr')

    pseudo <- min(df_species[df_species > 0])

    df_species_mod <- df_species + pseudo / 2

    gp.philr <- philr(t(df_species_mod), tree_sub, 
                      part.weights='enorm.x.gm.counts', 
                      ilr.weights='blw.sqrt')
    return(gp.philr)
}

# `find_significant_nodes` function identifies the significant nodes in a phylogenetic sub-tree based on 
# the results of a given phylogenetic ILR transformed data and additional metadata. This is useful for pinpointing 
# specific branches or nodes in the tree that might be driving observed differences between samples or groups.
#
# Args:
#   tree_sub: A phylogenetic sub-tree object. This tree typically represents a subset of the main phylogenetic tree.
# 
#   gp.philr: The results of the phylogenetic Isometric Log-Ratio (ILR) transformation, usually performed using `philr`.
#             This data can be used to understand the changes in abundance patterns across the tree hierarchy.
#
#   df_meta: A dataframe containing metadata for each sample in the study. This data might include details like 
#            sample origin, treatment group, time point, etc., and is used for grouping or categorizing samples 
#            during the significance analysis.
#
#   df_tax: A dataframe detailing the taxonomy associated with each species or taxonomic unit. It helps in interpreting 
#           the results by linking tree nodes to their corresponding taxonomic identities.

find_significant_nodes <- function(tree_sub, gp.philr, df_meta, df_tax) {

    tax <- df_tax
    
    df_sig <- NULL                                 
    for (i in 1:ncol(gp.philr)) {
        vals.pos <- gp.philr[,i][df_meta$Detected == TRUE]
        vals.neg <- gp.philr[,i][df_meta$Detected == FALSE]
        df_sig <- rbind(df_sig,data.frame(colnames(gp.philr)[i],mean(vals.pos), mean(vals.neg), wilcox.test(vals.pos, vals.neg)$p.value, name.balance(tree_sub, tax,colnames(gp.philr)[i])))
    }
                                              
    colnames(df_sig) <- c("Node","GA3+","GA3-","pvalue","balance")
    df_sig$bonferroni <- p.adjust(df_sig$pvalue, "bonferroni")
    df_sig$fdr <- p.adjust(df_sig$pvalue, "fdr")
    return(df_sig)
}

