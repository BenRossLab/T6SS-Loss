
from ete3 import Tree, CircleFace, RectFace, TreeStyle,NodeStyle, faces, AttrFace, TextFace

import tempfile
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm

import Analysis_DistToBase

from collections import defaultdict, Counter

from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()


# Allows rendering, https://github.com/etetoolkit/ete/issues/296
import os
os.environ['QT_QPA_PLATFORM']='offscreen'



def run_pastml(infile_tree = "../Data/RAxML_bestTree.simplifyaln", infile_characters = "../Data/ga3_struct_blast.csv", col_oi = "found"):
    outfile_pastml = tempfile.NamedTemporaryFile()
    outdir_pastml = tempfile.TemporaryDirectory()
    
    cmd = ["/home/averster/miniconda3/envs/t6ss/bin/pastml", "--tree", infile_tree, "--data", infile_characters, "--columns", col_oi, "--data_sep", ",", "-o", outfile_pastml.name, "--work_dir", outdir_pastml.name]
    subprocess.call(cmd)
    
    infile_tree_pastnames = Path(outdir_pastml.name) / "named.tree_{}.nwk".format(Path(infile_tree).with_suffix("").name)
    t = Tree(str(infile_tree_pastnames), format=1)
    df_reconstruction, found_dict = load_pastml_results(outfile_pastml.name)
    print(df_reconstruction)
    return df_reconstruction, found_dict, t


def load_pastml_results(infile_reconstruction):
    df_reconstruction = pd.read_csv(infile_reconstruction, sep = "\t")
    found_dict = {}
    for (i, dat) in df_reconstruction.iterrows():
        found_dict[dat["node"]] = dat["found"]
    return df_reconstruction, found_dict


def prep_tree_for_ape(tree_str):
    t = Tree(tree_str)

    i = 1
    for node in t.traverse("levelorder"):
        if node.name == "":
            node.name = "n{}".format(i)
            i += 1

    # Don't name the root
    t.name = ""
    tree_str = t.write(format=1)
    return t, tree_str


def do_reconstruction_ape(infile_tree, df_ga3):
    with open(infile_tree) as f:
        tree_str = f.read().rstrip("\n")

    t, tree_str = prep_tree_for_ape(tree_str)
    t.prune([x for x in t.iter_leaves() if x.name in df_ga3["name"].values])

    ro.globalenv['df_ga3'] = df_ga3
    ro.globalenv['tree_str'] = tree_str

    ro.r('''
    
        library('ape')
        rownames(df_ga3) <- df_ga3$name

        tree <- ape::read.tree(text=tree_str)
        tree <- root(tree, tree$tip.label[1])
        tree <- multi2di(tree)

        tree <- drop.tip(tree, tree$tip.label[!tree$tip.label%in%rownames(df_ga3)])
        stopifnot(all(tree$tip.label%in%rownames(df_ga3)))

        ga3 <- df_ga3[tree$tip.label,"found"]
        names(ga3) <- tree$tip.label
        ans <- ace(ga3, tree, type = "d", model="ARD")

        reconstruction <- apply(ans$lik.anc,1, function(x) colnames(ans$lik.anc)[which.max(x)])
        names(reconstruction) <- tree$node.label

        df_results <- data.frame(node=c(names(ga3),names(reconstruction)), found=c(ga3,reconstruction)) 
    ''')

    df_reconstruction = ro.globalenv['df_results']

    found_dict = {}
    for (i, dat) in df_reconstruction.iterrows():
        if dat["found"] == "TRUE":
            found_dict[dat["node"]] = True
        else:
            assert dat["found"] == "FALSE"
            found_dict[dat["node"]] = False

    return df_reconstruction, found_dict, t



def infer_events(t, found_dict):
    event_dict = {}
    for node in t.traverse("postorder"):
        parent = node.up
        if parent is None:
            continue
        if (found_dict[parent.name]) & (not found_dict[node.name]):
            event_dict[node.name] = "loss"
        elif (not found_dict[parent.name]) & (found_dict[node.name]):
            event_dict[node.name] = "gain"
        else:
            event_dict[node.name] = "standard"
    return event_dict

def prep_tree_for_plotting(t, event_dict, found_dict):
    
    nstyle_blank = NodeStyle()
    nstyle_blank["size"] = 0
    
    star_red = faces.ImgFace("/home/averster/Documents/RossLab/T6SSLoss/Figures/star_red.png")
    star_green = faces.ImgFace("/home/averster/Documents/RossLab/T6SSLoss/Figures/star_green.png")
    
    circle_size = 8
    
    for node in t.traverse("postorder"):
        if node.name == "root":
            continue
        node.set_style(nstyle_blank)
        if found_dict[node.name] & (node.is_leaf()):
            node.add_face(CircleFace(circle_size, "black","#616161"), column=0)
        elif (node.is_leaf()):
            node.add_face(CircleFace(circle_size, "black","white"), column=0)

        if node.name == "":
            continue

        if event_dict[node.name] == "gain":
            node.add_face(star_green, column=0, position="float")
        if event_dict[node.name] == "loss":
            node.add_face(star_red, column=0, position="float")

            
def add_node_names(t):
    i = 1
    for node in t.traverse("postorder"):
        if node.is_leaf():
            name_dict[node.name] = i
            node.name = i
            i+=1
    return name_dict


def add_node_names_dict(t, name_dict):
    for node in t.traverse("postorder"):
        if node.is_leaf():
            node.name = name_dict[node.name]


def layout(node):
    if node.is_leaf():
        N = AttrFace("name", fsize=10)
        faces.add_face_to_node(N, node, 1, position="branch-right")


def plot_tree_figure(infile_tree, df_ga3, outfile, name_dict=None, scale_use=80000, add_names = False):
    #df_reconstruction, found_dict, t = run_pastml(infile_tree)
    #print("First")
    #print(df_reconstruction)

    df_reconstruction, found_dict, t = do_reconstruction_ape(infile_tree, df_ga3)

    event_dict = infer_events(t, found_dict)
    prep_tree_for_plotting(t, event_dict, found_dict)
    
    if add_names:
        add_node_names_dict(t, name_dict)
        t.set_outgroup( t&"1" )
    
    ts = TreeStyle()
    ts.mode = "c" # draw tree in circular mode
    ts.show_leaf_name = False
    ts.scale = scale_use
    if add_names:
        ts.layout_fn = layout

    t.render(outfile, w=400, units="mm", tree_style=ts)

def load_metadata():
    df_meta = pd.read_csv("../Data/Bfragilis_metadata.csv",index_col=0)
    df_meta = df_meta[["strain"]]
    df_meta["genome"] = df_meta.index
    df_meta["number"] = np.arange(start=1,stop=df_meta.shape[0]+1,step=1)

    name_dict = {}
    for (i,dat) in df_meta.iterrows():
        name_dict[dat["genome"]] = str(dat["number"])
    return df_meta, name_dict

def bootstrap_count_gains(df_ga3, infile_bootstrap =  "Data/RAxML_bootstrap.fragilis_tree_markersv3_simplify_boot"):
    results = []
    with open(infile_bootstrap,"r") as f:
        for line in tqdm(f,total=100):
            line = line.rstrip("\n")
            tree_str_file = tempfile.NamedTemporaryFile()
            tree_str_file.write(line.encode())
            df_reconstruction, found_dict, t = do_reconstruction_ape(tree_str_file.name, df_ga3)
            event_dict = infer_events(t, found_dict)
            vals = Counter(event_dict.values())
            loss = 0
            if 'loss' in vals:
                loss=vals['loss']
            gain = 0
            if 'gain' in vals:
                gain=vals['gain']
            results.append([gain,loss])
    df_results = pd.DataFrame(results, columns=["gain","loss"])
    return df_results


if __name__ == "__main__":
    np.random.seed(7)
    infile_ga3 = "Data/ga3_in_tree_annotated.xlsx"
    infile_bootstrap = "/home/averster/Documents/RossLab/T6SSLoss/Data/RAxML_bootstrap.fragilis_tree_markersv3_simplify_boot"
    
    df_ga3 = pd.read_excel(infile_ga3)
    df_ga3 = df_ga3[df_ga3['found'].notna()]
    df_ga3 = df_ga3[["name","n_found", "n_total","found"]]
    df_ga3["found"] = df_ga3["found"] > 0

    Analysis_DistToBase.calculate_parameters(df_ga3, use_median=True)
    df_params_bootstrap = Analysis_DistToBase.calculate_distances_bootstrap(df_ga3, use_median=True)
    print("Fraction of bootstraps with a significant effecct {}, median pvalue {}".format(np.sum(df_params_bootstrap["pval"] < 0.05) / df_params_bootstrap.shape[0], np.median(df_params_bootstrap["pval"].values)))

    df_meta, name_dict = load_metadata()
    plot_tree_figure("Data/RAxML_bestTree.fragilis_tree_core98_simplify", df_ga3, "fragilis_tree_t6ss_loss_core98.pdf")
    plot_tree_figure("Data/RAxML_bestTree.fragilis_tree_core95_simplify", df_ga3, "fragilis_tree_t6ss_loss_core95.pdf")
    plot_tree_figure("Data/RAxML_bestTree.fragilis_tree_markersv3_simplify", df_ga3, "fragilis_tree_t6ss_loss_metaphlanv3.pdf", scale_use=45000)
    df_results = bootstrap_count_gains(df_ga3)
    print("means")
    print(df_results.mean(0))
    print("standard deviation")
    print(df_results.std(0))
