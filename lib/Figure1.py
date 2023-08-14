
import Figure1_Data
from Figure1_SpeciesToGenusTree import SpeciesToGenusTree

from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()

import pandas as pd

from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace


def make_abundance_scatter_plot(outfile):
    """
    Makes Figure 1 B and C

    """

    ro.globalenv['outfile'] = outfile

    ro.r('''

    library('ggplot2')
    library('extrafont')
    source("lib/Figure1_Analysis.R")

    t <- improve_plot()


    p <- ggplot(df_ab_sub,aes(x = log10(fraction_bfr), y=log10(fraction_t6 + 1e-7))) + geom_point(size=3, pch=21, fill="steelblue") + theme_bw(16) + theme(text=element_text(size=16, color="black", family="Arial")) + ylab("Log10 GA3 Structural genes") + xlab("Log10 B. fragilis") + t
    
    pdf(outfile, width=4,height=3.5, family="Arial")
    print(p)
    dev.off()

    ''')

def make_fraction_ga3_minus_barplot(df_ga3, outfile):
    """
    Makes Figure 1A
    """
    ro.globalenv['df_ga3'] = df_ga3
    ro.globalenv['outfile'] = outfile
    ro.r('''
        library('tidyverse')

        source("lib/Figure1_Analysis.R")

        r1 <- df_ga3 %>% filter(type == "Infant") %>% select(-dataset, -type, -fraction_minus) %>% apply(2, sum)
        r2 <- df_ga3 %>% filter(type == "Adult") %>% select(-dataset,-type,-fraction_minus) %>% apply(2, sum)

        print(bind_rows(r1, r2) %>% fisher.test())

        t <- improve_plot()
        df_ga3$dataset <- factor(df_ga3$dataset, levels=c("DIABIMMUNE","TEDDY","HMP","MetaHIT", "Yachida"))

        pdf(outfile, width=3.7, height=5.6)
        p <- ggplot(df_ga3, aes(x = dataset, y = fraction_minus)) + geom_bar(stat="identity", aes(fill=type), color="black") + theme_bw(16) +  scale_y_continuous(labels = scales::percent) + ylab("Percentage GA3-") + xlab("") + scale_fill_brewer(palette = "Set1", name="") + t + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        print(p)
        dev.off()
    ''')


def do_pcoa(tree_str, df_metaphlan, df_meta, outfile_plot, add_loading=False, taxa_oi=None):
    """
    Makes Figure 1 D and E
    """
    ro.globalenv['df_metaphlan'] = df_metaphlan
    ro.globalenv['df_meta'] = df_meta
    ro.globalenv['tree_str'] = tree_str
    ro.globalenv['outfile_plot'] = outfile_plot
    ro.globalenv['add_loading'] = add_loading
    if taxa_oi is not None:
        ro.globalenv['taxa_oi'] = taxa_oi

    ro.r('''
        if (!exists("taxa_oi"))
            taxa_oi <- NULL
        source("lib/Figure1_Analysis.R")
        pcoa.results <- run_pcoa(tree_str, df_metaphlan, df_meta)
        plot_pcoa(pcoa.results, df_metaphlan, df_meta, outfile_plot, add_loading, taxa_oi)

    ''')

def do_philr(tree_str, df_metaphlan, df_taxa, df_meta, outfile_fig):
    """
    Executes a phiLR (Phylogenetic Isometric Log-Ratio Transform) analysis on the given data using R's 'philr' package.

    This function uses R's philr package to perform a phiLR transformation on the input metagenomic data. 
    The function then identifies significant nodes in the phylogenetic tree using the transformed data.

    Parameters:
    - tree_str (str): A string representation of the phylogenetic tree.
    - df_metaphlan (DataFrame): A pandas DataFrame containing metagenomic data, typically derived from MetaPhlAn.
    - df_taxa (DataFrame): A taxa table formatted for the philr R package.
    - df_meta (DataFrame): Metadata associated with the samples in df_metaphlan.
    - outfile_fig (str): Path where the resulting figure should be saved.
    """
    ro.globalenv['df_metaphlan'] = df_metaphlan
    ro.globalenv['df_taxa'] = df_taxa


    ro.globalenv['tree_str'] = tree_str
    ro.globalenv['df_meta'] = df_meta
    ro.globalenv['outfile_fig'] = outfile_fig

    ro.r('''
        library('tidyverse')
        source("lib/Figure1_Analysis.R")
        save(tree_str, df_metaphlan, df_meta, df_taxa, file = "~/Desktop/debug.RData")

        tree_sub <- prep_tree_for_philr(tree_str, df_metaphlan)
        gp.philr <- run_philr(df_metaphlan, tree_sub)
        df_sig <- find_significant_nodes(tree_sub, gp.philr, df_meta, df_taxa)        
        tree_sub_str <- write.tree(tree_sub)
    ''')
    tree_sub_str = ro.globalenv['tree_sub_str']
    df_sig = ro.globalenv['df_sig']
    return df_sig, tree_sub_str


def make_tree_plot(df_f, tree_sub_str, outfile_treefig):
    """
    Makes Figure 1 F
    """

    
    def layout(node):
        if node.is_leaf():
            faces.add_face_to_node(AttrFace("name", fstyle="italic", fsize=14,ftype="Arial"), node, column=0)
        if node.name in nodes_oi:
            node.name = node.name.replace("Node","")
            faces.add_face_to_node(AttrFace("name", fsize=8,fstyle="bold", ftype="Arial"), node, column=1, position="float")

    t = Tree(tree_sub_str[0], format=1)

    nodes_oi = df_f.query('fdr < 0.05')["Node"].values

    line_width=1.25
    nst1 = NodeStyle()
    nst1["fgcolor"] = "#c33b3b"
    nst1["size"] = 16
    nst1["vt_line_width"] = line_width
    nst1["hz_line_width"] = line_width

    nst2 = NodeStyle()
    nst2["size"] = 0
    nst2["vt_line_width"] = line_width
    nst2["hz_line_width"] = line_width

    for n in t.traverse():
        if n.name in nodes_oi:
            n.set_style(nst1)
        else:
            n.set_style(nst2)
        if "s__" in n.name:
            n.name = n.name.replace("s__","")
            n.name = n.name.replace("Bacteroides_","B. ")
            
    ts = TreeStyle()
    ts.layout_fn = layout
    ts.show_leaf_name = False
    t.render(outfile_treefig, w=200, tree_style = ts)


def make_taxa_table(df_metaphlan):
    """
    Generates a taxa table suitable for use with the philr R package.


    Parameters:
    - df_metaphlan (DataFrame): A pandas DataFrame containing metagenomic data, usually derived from MetaPhlAn.

    """
    df_taxa = pd.DataFrame(pd.DataFrame([x.split("|") for x in df_metaphlan.index]))
    df_metaphlan.index = [x.split("|")[-1] for x in df_metaphlan.index]
    df_taxa.index = df_metaphlan.index
    df_taxa.columns = ["ta{}".format(i+1) for i in range(len(df_taxa.columns))]
    return df_taxa, df_metaphlan


if __name__ == "__main__":
    counts_cutoff = 100
    # Figure 1A
    df_plot = Figure1_Data.make_fraction_ga3_minus_dataframe(counts_cutoff)
    make_fraction_ga3_minus_barplot(df_plot, "Figure1A.pdf")

    # Figure 1B-C
    Figure1_Data.send_gene_abundance_data_to_R("diabimmune")
    make_abundance_scatter_plot("Figure1B.pdf")
    
    Figure1_Data.send_gene_abundance_data_to_R("HMP")
    make_abundance_scatter_plot("Figure1C.pdf")

    # Figure 1D-E
    df_genus_diabimmune, df_meta_diabimmune, tree_str = Figure1_Data.load_diabimmune_metaphlan()
    tree_condenser = SpeciesToGenusTree(tree_str, df_genus_diabimmune.index)
    tree_str_genus = tree_condenser.write_tree()
    do_pcoa(tree_str_genus, df_genus_diabimmune, df_meta_diabimmune, "Figure1D.pdf", add_loading=True, taxa_oi="g__Bacteroides")

    df_genus_hmp, df_species_hmp, df_bacteroides_hmp, df_meta_hmp, tree_str = Figure1_Data.load_hmp_metaphlan()
    
    tree_condenser = SpeciesToGenusTree(tree_str, df_genus_hmp.index)
    tree_str_genus = tree_condenser.write_tree()
    do_pcoa(tree_str_genus, df_genus_hmp, df_meta_hmp, "Figure1E.pdf", add_loading=True, taxa_oi="g__Bacteroides")

    # Figure 1F
    df_taxa_species, df_species_hmp = make_taxa_table(df_species_hmp)
    df_sig, tree_sub_str = do_philr(tree_str, df_bacteroides_hmp, df_taxa_species, df_meta_hmp, "outfile_philr.pdf")
    make_tree_plot(df_sig, tree_sub_str, "Figure1F.svg")
