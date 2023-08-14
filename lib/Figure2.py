import numpy as np
import pandas as pd
import Figure2_Analysis
from ete3 import CircleFace, TreeStyle,NodeStyle, faces, AttrFace



def prep_tree_for_plotting(t, event_dict, found_dict, infile_red_star = "Data/star_red.png", infile_green_star = "Data/star_green.png"):
    """
    Prepare an ete3 tree for plotting by adding custom faces to indicate gains and losses

    Parameters:
    - t: The ete3 tree to prepare for plotting.
    - event_dict: A dictionary that maps leaf node names to event labels.
    - found_dict: A dictionary that maps leaf node names to a boolean (found or not found).
    - infile_red_star: Path to the image file for a red star symbol.
    - infile_green_star: Path to the image file for a green star symbol.

    """

    nstyle_blank = NodeStyle()
    nstyle_blank["size"] = 0
    
    star_red = faces.ImgFace(infile_red_star)
    star_green = faces.ImgFace(infile_green_star)
    
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


def add_node_names_dict(t, name_dict):
    """
    Update leaf node names in a tree using a dictionary.

    Parameters:
    - t: The ete3 tree to modify.
    - name_dict: A dictionary mapping original leaf node names to new names.
    """

    for node in t.traverse("postorder"):
        if node.is_leaf():
            # Update the leaf node name using the dictionary
            node.name = name_dict[node.name]


def layout(node):
    """
    Custom layout function for ete3 tree visualization.

    Parameters:
    - node: A node in the ete3 tree.
    """
    if node.is_leaf():
        N = AttrFace("name", fsize=10)
        faces.add_face_to_node(N, node, 1, position="branch-right")


def plot_tree_figure(infile_tree, df_ga3, outfile, name_dict=None, scale_use=80000, add_names = False):
    """
    Plot a tree figure with optional additional information.

    Parameters:
    - infile_tree (str): Path to the input tree file.
    - df_ga3 (pandas.DataFrame): DataFrame containing GA3 +/- data.
    - outfile (str): Path to save the output tree figure.
    - name_dict (dict, optional): Dictionary to map names to labels.
    - scale_use (int, optional): Scaling factor for the tree layout.
    - add_names (bool, optional): If True, add names to the figure.
    """
    df_reconstruction, found_dict, t = Figure2_Analysis.do_reconstruction_ape(infile_tree, df_ga3)

    event_dict = Figure2_Analysis.infer_events(t, found_dict)
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


if __name__ == "__main__":
    np.random.seed(7)
    infile_ga3 = "Data/ga3_in_tree_annotated.xlsx"
    infile_bootstrap = "Data/RAxML_bootstrap.fragilis_tree_markersv3_simplify_boot"
    infile_tree = "Data/RAxML_bestTree.fragilis_tree_markersv3_simplify"
    
    # Load up the GA3 annotation data
    df_ga3 = pd.read_excel(infile_ga3)
    df_ga3 = df_ga3[df_ga3['found'].notna()]
    df_ga3 = df_ga3[["name","n_found", "n_total","found"]]
    df_ga3["found"] = df_ga3["found"] > 0

    Figure2_Analysis.calculate_parameters(df_ga3, infile_tree, use_median=True)
    df_params_bootstrap = Figure2_Analysis.calculate_distances_bootstrap(df_ga3, infile_bootstrap, use_median=True)
    print("Fraction of bootstraps with a significant effect {}, median pvalue {}".format(np.sum(df_params_bootstrap["pval"] < 0.05) / df_params_bootstrap.shape[0], np.median(df_params_bootstrap["pval"].values)))

    plot_tree_figure(infile_tree, df_ga3, "Figure2B.pdf", scale_use=45000)
    df_results = Figure2_Analysis.bootstrap_count_gains_losses(df_ga3)
    print("means of bootstrap losses and gains")
    print(df_results.mean(0))
    print("standard deviation of bootstrap losses and gains")
    print(df_results.std(0))
