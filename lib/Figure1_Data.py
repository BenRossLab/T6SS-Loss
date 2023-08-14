
from pathlib import Path

import re
from Bio import SeqIO

import numpy as np
import pandas as pd

from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()

def load_count_tables(infile_bfr_counts, infile_t6_counts, counts_cutoff, split_sample_names=False):
    """
    Load the bowtie count tables for B. fragilis markers and GA3 T6SS genes.

    Parameters:
    - infile_bfr_counts: Path to the B. fragilis count table file.
    - infile_t6_counts: Path to the GA3 T6SS count table file.
    - counts_cutoff: Count threshold for including data.
    - split_sample_names: Flag to split sample names based on certain delimiter.
    """

    df_bfr = pd.read_csv(infile_bfr_counts)
    df_bfr = df_bfr[["Tag","total"]].rename(columns={"Tag":"sample"})
    if split_sample_names:
        df_bfr["sample"] = [x.split("_")[0] for x in df_bfr["sample"]]

    df_bfr = df_bfr.query('total >= {}'.format(counts_cutoff))

    samples_use = df_bfr["sample"].values

    df_t6 = pd.read_csv(infile_t6_counts)
    df_t6 = df_t6[["Tag","total"]].rename(columns={"Tag":"sample"})
    if split_sample_names:
        df_t6["sample"] = [x.split("_")[0] for x in df_t6["sample"]]

    samples_use = [x for x in samples_use if x in df_t6["sample"].values]
    df_t6 = df_t6.query('sample in @samples_use')
    assert df_t6.shape[0] == len(samples_use)
    return df_bfr, df_t6


def load_seq_list_from_fasta(infile_fasta):
    """
    Load sequences from a FASTA file.

    Parameters:
    - infile_fasta: Path to the input FASTA file.
    """

    seq_list = []
    with open(infile_fasta) as f:
        for seq in SeqIO.parse(f, "fasta"):
            seq_list.append(seq)
    return seq_list


def identify_ga3_minus(df_bfr, df_t6, fraction_include = 0.1, infile_bfr = "Data/MetaPhlanMarkers_RepGenome.fa", infile_t6 = "Data/GA3_structural.fasta"):
    """
    Identify samples that are missing GA3 T6SS genes.

    Parameters:
    - df_bfr: Dataframe containing B. fragilis markers count data.
    - df_t6: Dataframe containing GA3 T6SS genes count data.
    - fraction_include: Fraction of expected number of counts, below which we consider T6SS to be missing (default = 0.1).
    - infile_bfr: Path to the input B. fragilis markers FASTA file (default = "Data/MetaPhlanMarkers_RepGenome.fa").
    - infile_t6: Path to the input GA3 T6SS genes FASTA file (default = "Data/GA3_structural.fasta").
    """

    seqs_bfr = load_seq_list_from_fasta(infile_bfr)
    total_bfr_len = np.sum([len(x) for x in seqs_bfr])

    seqs_t6 = load_seq_list_from_fasta(infile_t6)
    total_t6_len = np.sum([len(x) for x in seqs_t6])

    df_bfr.index = df_bfr['sample']
    detected = []
    expected_T6SS_vals = []
    for (i, dat) in df_t6.iterrows():
        counts_bfr = df_bfr.loc[dat["sample"],"total"]

        expected_T6SS = counts_bfr * total_t6_len / total_bfr_len
        expected_T6SS_vals.append(expected_T6SS)
        if dat["total"] > expected_T6SS * fraction_include:
            detected.append(True)
        else:
            detected.append(False)

    df_t6["Detected"] = detected
    df_t6["expected_T6SS"] = expected_T6SS_vals

    return df_t6


def group_by_individual(df_t6):
    """
    Group GA3 T6SS gene presence by individual.

    Parameters:
    - df_t6: Dataframe containing GA3 T6SS gene count data.
    """

    r = []
    for (indv, dat) in df_t6.groupby("Individual"):
        if len(dat["Detected"].unique()) > 1:
            continue
        r.append([indv, dat["Detected"].unique()[0]])
    df_r = pd.DataFrame(r, columns = ["individual","Detected"])
    return df_r



def load_yachida_data(counts_cutoff = 100, ga3_minus_cutoff = 0.1):
    """
    Load the data from the Yachida et al. study.

    Parameters:
    - counts_cutoff: Cutoff value for including samples passed to load_count_tables()
    - ga3_minus_cutoff: fraction cutoff passed to identify_ga3_minus()
    """

    infile_yachida = "Data/Yachida_metadata.xlsx"
    df_names = pd.read_excel(infile_yachida)
    df_names["Subject_ID"] = [int(re.search("human feces_([0-9]+)",x).group(1)) for x in df_names["sample_title"].values]

    df_meta = pd.read_excel("Data/41591_2019_458_MOESM3_ESM.xlsx", sheet_name="Table_S2-1", skiprows=2)
    subj_oi = df_meta.query('Group == "Healthy"')["Subject_ID"].values.tolist()
    df_names = df_names.query('Subject_ID in @subj_oi')

    df_bfr, df_t6 = load_count_tables(infile_bfr_counts = "Data/Yachida2019_counts_BfrMarkers_combined.csv", infile_t6_counts = "Data/Yachida2019_counts_T6Struct_combined.csv", counts_cutoff=counts_cutoff)
    df_t6 = identify_ga3_minus(df_bfr, df_t6, ga3_minus_cutoff)

    return df_t6, ["Yachida", "Adult",np.sum(~df_t6["Detected"]), np.sum(df_t6["Detected"]), np.sum(~df_t6["Detected"]) / df_t6.shape[0]]


def load_hmp_data(counts_cutoff = 100, ga3_minus_cutoff = 0.1):
    """
    Load the data from HMP.

    Parameters:
    - counts_cutoff: Cutoff value for including samples passed to load_count_tables()
    - ga3_minus_cutoff: fraction cutoff passed to identify_ga3_minus()
    """

    df_bfr, df_t6 = load_count_tables(infile_bfr_counts = "Data/HMP_counts_BfrMarkers_combined.csv", infile_t6_counts = "Data/HMP_counts_GA3Struct_combined.csv", counts_cutoff=counts_cutoff)
    df_t6 = identify_ga3_minus(df_bfr, df_t6, ga3_minus_cutoff)

    df_meta = pd.read_csv("Data/hmp_combined_metadata.csv", sep="\t")
    df_meta = df_meta.rename(columns={"SRS":"sample"})

    df_t6_full = df_t6.merge(df_meta)
    df_t6_indv = group_by_individual(df_t6_full)
    return df_t6, ["HMP", "Adult",np.sum(~df_t6_indv["Detected"]), np.sum(df_t6_indv["Detected"]), np.sum(~df_t6_indv["Detected"]) / df_t6_indv.shape[0]]


def load_teddy(counts_cutoff=100, ga3_minus_cutoff=0.1):
    """
    Load the data from TEDDY.

    Parameters:
    - counts_cutoff: Cutoff value for including samples passed to load_count_tables()
    - ga3_minus_cutoff: fraction cutoff passed to identify_ga3_minus()
    """

    dataset = "TEDDY"
    df_bfr, df_t6 = load_count_tables("Data/{}_counts_BfrMarkers_combined.csv".format(dataset), "Data/{}_counts_GA3Struct_combined.csv".format(dataset), counts_cutoff=counts_cutoff)
    df_t6 = identify_ga3_minus(df_bfr, df_t6, ga3_minus_cutoff)

    infile_meta = "Data/TEDDY_10_Percent_Rand_Selected_220113.txt"
    df_teddy_meta = pd.read_csv(infile_meta, sep=" ")
    df_t6_full = df_t6.merge(df_teddy_meta[["BioSample","dbGaP_Subject_ID","collinterval"]].rename(columns={"collinterval":"age","dbGaP_Subject_ID":"Individual","BioSample":"sample"}), how="left")
    df_t6_full_indv = group_by_individual(df_t6_full)
    return df_t6_full, ["TEDDY","Infant",np.sum(~df_t6_full_indv["Detected"]), np.sum(df_t6_full_indv["Detected"]), np.sum(~df_t6_full_indv["Detected"]) / df_t6_full_indv.shape[0]]


def load_metahit(counts_cutoff=100, ga3_minus_cutoff=0.1):
    """
    Load the data from MetaHIT.

    Parameters:
    - counts_cutoff: Cutoff value for including samples passed to load_count_tables()
    - ga3_minus_cutoff: fraction cutoff passed to identify_ga3_minus()
    """

    dataset = "MetaHit"
    df_bfr, df_t6 = load_count_tables("Data/{}_counts_BfrMarkers_combined.csv".format(dataset), "Data/{}_counts_GA3Struct_combined.csv".format(dataset), counts_cutoff=counts_cutoff)
    df_t6 = identify_ga3_minus(df_bfr, df_t6, ga3_minus_cutoff)
    return df_t6, ["MetaHIT","Adult",np.sum(~df_t6["Detected"]), np.sum(df_t6["Detected"]), np.sum(~df_t6["Detected"]) / df_t6.shape[0]]


def load_data_diabimmune(counts_cutoff=100, ga3_minus_cutoff=0.1):
    """
    Load the data from DIABIMMUNE.

    Parameters:
    - counts_cutoff: Cutoff value for including samples passed to load_count_tables()
    - ga3_minus_cutoff: fraction cutoff passed to identify_ga3_minus()
    """

    days_to_months = 365 / 12
    
    dataset = "Yassour"
    df_bfr_yassour, df_t6_yassour = load_count_tables("Data/{}_counts_BfrMarkers_combined.csv".format(dataset), "Data/{}_counts_GA3Struct_combined.csv".format(dataset), counts_cutoff=counts_cutoff)
    df_t6_yassour = identify_ga3_minus(df_bfr_yassour, df_t6_yassour, ga3_minus_cutoff)
    
    df_t6_yassour["age"] = [float(x.split("_")[1]) for x in df_t6_yassour["sample"]]
    df_t6_yassour["Individual"] = [str(x.split("_")[0]) for x in df_t6_yassour["sample"]]
    
    dataset = "Kostic"
    df_bfr_kostic, df_t6_kostic = load_count_tables("Data/{}_counts_BfrMarkers_combined.csv".format(dataset), "Data/{}_counts_GA3Struct_combined.csv".format(dataset), counts_cutoff=counts_cutoff)
    df_t6_kostic = identify_ga3_minus(df_bfr_kostic, df_t6_kostic, ga3_minus_cutoff)

    infile_meta = "Data/Kostic_metadata.xlsx"
    df_meta_kostic = pd.read_excel(infile_meta, engine='openpyxl')

    df_meta_kostic_sub = df_meta_kostic[["Gid_shotgun","Age_at_Collection","Subject_ID"]].rename(columns={"Gid_shotgun":"sample","Age_at_Collection":"age", "Subject_ID":"Individual"})
    df_meta_kostic_sub["age"] = df_meta_kostic_sub["age"] / days_to_months
    df_t6_kostic = df_t6_kostic.merge(df_meta_kostic_sub, on="sample")
    
    dataset = "Vatanen"
    df_bfr_vatanen, df_t6_vatanen = load_count_tables("Data/{}_counts_BfrMarkers_combined.csv".format(dataset), "Data/{}_counts_GA3Struct_combined.csv".format(dataset), counts_cutoff=counts_cutoff)
    df_t6_vatanen = identify_ga3_minus(df_bfr_vatanen, df_t6_vatanen, ga3_minus_cutoff)

    infile_meta = "Data/VatanenMetadata.csv"
    df_meta_vatanen = pd.read_csv(infile_meta)

    df_meta_vatanen_sub = df_meta_vatanen[["gid_wgs","age_at_collection","subjectID"]].rename(columns={"gid_wgs":"sample","age_at_collection":"age", "subjectID":"Individual"})
    df_meta_vatanen_sub["age"] = df_meta_vatanen_sub["age"] / days_to_months
    df_t6_vatanen = df_t6_vatanen.merge(df_meta_vatanen_sub, on="sample")

    df_t6_full = df_t6_vatanen.append(df_t6_kostic).append(df_t6_yassour)
    df_t6_full_indv = group_by_individual(df_t6_full)
    return df_t6_full, ["DIABIMMUNE","Infant",np.sum(~df_t6_full_indv["Detected"]), np.sum(df_t6_full_indv["Detected"]), np.sum(~df_t6_full_indv["Detected"]) / df_t6_full_indv.shape[0]]


def load_species_tree(infile_tree = "Data/mpa_v30_CHOCOPhlAn_201901_species_tree.nwk"):
    """
    Load a species tree from a file provided by metaphlan
    https://github.com/biobakery/MetaPhlAn/tree/master/metaphlan/utils

    Parameters:
    - infile_tree: Path to the input species tree file.
    """

    with open(infile_tree) as f:
        tree_str = next(f).rstrip("\n")

    tree_str = re.sub("GCA_[0-9]+\|","",tree_str)
    return tree_str


def remove_missing_species(df_metaphlan, tree_str):
    """
    Remove species from the DataFrame that are not present in the provided tree.

    Parameters:
    - df_metaphlan: DataFrame containing Metaphlan data.
    - tree_str: Newick string representation of the species tree.
    """

    missing = [x for x in df_metaphlan.index if x not in tree_str]
    return df_metaphlan.loc[[x not in missing for x in df_metaphlan.index],:]


def combine_dfs(df_list):
    """
    Merges a list of Metaphlan dataframes and fills in missing species with zeros.

    Parameters:
    - df_list: List of dataframes to be merged.
    """

    df_full = df_list[0]
    for i in range(1, len(df_list)):
        df_full = df_full.merge(df_list[i],left_index=True,right_index=True,how="outer").fillna(0)
    return df_full


def load_combined_metaphlan_file(infile_excel, simplify_taxa = True):
    """
    Load a combined Metaphlan file from an Excel format and process the data.

    Parameters:
    - infile_excel: Path to the input Excel file.
    - simplify_taxa: Flag indicating whether to simplify taxa names.
    """

    df = pd.read_excel(infile_excel, engine='openpyxl')
    if simplify_taxa:
        df.index = [x.split("|")[-1] for x in df["clade"]]
    else:
        df.index = df["clade"]
    df = df.drop("clade", axis=1)
    df.columns = [re.sub("_metaphlan.+","",x) for x in df.columns]

    return df


def load_diabimmune_metaphlan():
    """
    Load Metaphlan genus abundance data from DIABIMMUNE,
    selecting samples with sufficient T6SS.
    """

    df_t6_full,_ = load_data_diabimmune()
    tree_str = load_species_tree()

    df_genus = load_combined_metaphlan_file("Data/Vatanen_metaphlan_combined_g_relative_abundance.xlsx")
    df_genus_vatanen = remove_missing_species(df_genus, tree_str)

    df_genus = load_combined_metaphlan_file( "Data/Kostic_metaphlan_combined_g_relative_abundance.xlsx")
    df_genus_kostic = remove_missing_species(df_genus, tree_str)

    df_genus = load_combined_metaphlan_file("Data/Yassour_metaphlan_combined_g_relative_abundance.xlsx")
    df_genus_yassour = remove_missing_species(df_genus, tree_str)

    df_genus_diabimmune = combine_dfs([df_genus_vatanen, df_genus_kostic, df_genus_yassour])
    df_genus_diabimmune = df_genus_diabimmune / df_genus_diabimmune.sum(0)
    
    df_t6_full = df_t6_full.query('(age >= 12) & (age <= 24)')
    
    df_genus_diabimmune = df_genus_diabimmune[df_t6_full["sample"].values]
    return df_genus_diabimmune, df_t6_full, tree_str


def load_hmp_metaphlan():
    """
    Load Metaphlan genus abundance data from HMP,
    selecting samples with sufficient T6SS.
    """

    df_t6, _ = load_hmp_data() 
    tree_str = load_species_tree()

    df_species = load_combined_metaphlan_file("Data/HMP_metaphlan_combined_s_relative_abundance.xlsx", simplify_taxa=False)

    df_species = remove_missing_species(df_species, tree_str)
    df_species = df_species[df_t6["sample"].values]
    df_species = df_species / df_species.sum(0)

    df_genera = load_combined_metaphlan_file("Data/HMP_metaphlan_combined_g_relative_abundance.xlsx", simplify_taxa=False)
    df_genera = df_genera[df_t6["sample"].values]
    df_genera = df_genera / df_genera.sum(0)
    df_genera.index = [x.split("|")[-1] for x in df_genera.index]

    df_bacteroides = df_species.query('index.str.contains("Bacteroides")', engine="python")
    df_bacteroides.index = [x.split("|")[-1] for x in df_bacteroides.index]

    return df_genera, df_species, df_bacteroides, df_t6[["sample","Detected"]], tree_str


def make_fraction_ga3_minus_dataframe(counts_cutoff):
    """
    Calculate the values for Figure 1A.

    Parameters:
    - counts_cutoff: Count cutoff for B. fragilis positive samples
    """

    _, row_yachida = load_yachida_data(counts_cutoff=counts_cutoff)
    _, row_teddy = load_teddy(counts_cutoff=counts_cutoff)
    _, row_diabimmune = load_data_diabimmune(counts_cutoff=counts_cutoff)
    _, row_hmp = load_hmp_data(counts_cutoff=counts_cutoff)
    _, row_metahit = load_metahit(counts_cutoff=counts_cutoff)
    df_plot = pd.DataFrame([row_teddy, row_diabimmune, row_hmp, row_metahit, row_yachida], columns=["dataset","type","ga3_minus","ga3_plus","fraction_minus"])
    return df_plot


def load_normalized_counts(dataset = "Vatanen", min_bfr=100):
    """
    Load normalized counts for the specified dataset.

    Parameters:
    - dataset: Name of the dataset to load (default: "Vatanen").
    - min_bfr: Minimum count for B. fragilis to consider a sample (default: 100).
    """

    indir = Path("Data/")
    infile = indir / "{}_counts_BfrMarkers_combined.csv".format(dataset)
    df_counts_bfr = pd.read_csv(infile)
    df_counts_bfr = df_counts_bfr.query('total >= {}'.format(min_bfr))
    infile = indir / "{}_counts_GA3Struct_combined.csv".format(dataset)
    df_counts_t6 = pd.read_csv(infile)

    infile_metaphlan = indir / "{}_metaphlan_combined_s_estimated_number_of_reads_from_the_clade.xlsx".format(dataset)
    df_metaphlan = pd.read_excel(infile_metaphlan)
    df_metaphlan.columns = [x.split("_")[0] for x in df_metaphlan.columns]

    df_counts_bfr = df_counts_bfr[["Tag","total"]].merge(df_metaphlan.iloc[0,1:], left_on="Tag",right_index=True)
    df_counts_bfr.columns = ["Sample","counts_bfr", "total_counts"]
    df_counts_bfr["fraction_bfr"] = df_counts_bfr["counts_bfr"] / df_counts_bfr["total_counts"]
    
    df_counts_t6 = df_counts_t6[["Tag","total"]].merge(df_metaphlan.iloc[0,1:], left_on="Tag",right_index=True)
    df_counts_t6.columns = ["Sample","counts_t6", "total_counts"]
    df_counts_t6["fraction_t6"] = df_counts_t6["counts_t6"] / df_counts_t6["total_counts"]

    df_ab = df_counts_bfr.merge(df_counts_t6.drop(columns=["total_counts"]), on=["Sample"])
    return df_ab

def send_gene_abundance_data_to_R(data="diabimmune"):
    """
    Send gene abundance data to R for plotting Figure1BC.

    Args:
    - data: Name of the data to send (default: "diabimmune").
    """
    if data == "diabimmune":
        df_ab_vatanen = load_normalized_counts("Vatanen")
        df_ab_kostic = load_normalized_counts("Kostic")
        df_ab_yassour = load_normalized_counts("Yassour")
        df_ab = df_ab_vatanen.append(df_ab_kostic).append(df_ab_yassour)
    elif data == "HMP":
        df_ab = load_normalized_counts("HMP")
    elif data == "TEDDY":
        df_ab = load_normalized_counts("TEDDY")

    df_ab_sub = df_ab[["Sample","fraction_bfr","fraction_t6", "counts_bfr", "counts_t6"]]
    df_ab_sub["Sample"] = df_ab_sub["Sample"].astype(str)
    df_ab_sub["fraction_bfr"] = df_ab_sub["fraction_bfr"].astype(float)
    df_ab_sub["fraction_t6"] = df_ab_sub["fraction_t6"].astype(float)
    ro.globalenv['df_ab_sub'] = df_ab_sub
