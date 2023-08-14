from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()

import tempfile
from collections import defaultdict, Counter
import pandas as pd
import numpy as np
from ete3 import Tree
from tqdm import tqdm



def R_tree_reconstruction_randomization(tree_str, df_ga3, n_reps=1000):
    """
    Performs ancestral state reconstruction on the main tree and conducts a randomization test using the ace function in R.

    The function uses ancestral character estimation (ace) from R's 'ape' package to reconstruct ancestral states 
    on a provided phylogenetic tree. Furthermore, a randomization test is executed on the reconstructed states to 
    gauge the robustness of the ancestral state inferences.

    Parameters:
    - tree_str (str): A string representation of the phylogenetic tree on which ancestral state reconstruction will be performed.
    - df_ga3 (DataFrame): A pandas DataFrame containing the GA3 +/- data for the known species.
                          Used in ancestral state reconstruction.
    - n_reps (int, optional): The number of randomization replicates to perform. Default is 1000.

    """

    ro.globalenv['tree_str'] = tree_str
    ro.globalenv['df_ga3'] = df_ga3
    ro.globalenv['n_reps'] = n_reps
    
    ro.r('''
    
        set.seed(7)

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

        df_results <- data.frame(names=names(c(ga3,reconstruction)), reconstruction=c(ga3,reconstruction)) 

        Q <- ans$rates
        for (i in 1:n_reps) {
            # https://lukejharmon.github.io/pcm/chapter7_introdiscrete
            # https://lukejharmon.github.io/pcm/chapter8_fitdiscrete/
            # Diagonals are ignored

            randomization <- rTraitDisc(tree, model = matrix(c(0, Q, 0), 2), ancestor=TRUE, states=c("TRUE","FALSE"))
            stopifnot(all(names(c(ga3,reconstruction)) == names(randomization)))
            df_results[,paste0("randomization",i)] <- randomization
        }
                      
       ''')
    
    return ro.globalenv['df_results']


def convert_r_results(df_results):

    """
    Converts the results from ancestral state reconstruction into gain or loss events.

    Given a DataFrame that contains results from the ancestral state reconstruction (from the ace function 
    in R's 'ape' package), this function identifies and categorizes branches in the tree as having experienced 
    a gain or loss event.

    Parameters:
    - df_results (DataFrame): A DataFrame containing the results from the ancestral state reconstruction.
    """

    found_dict = defaultdict(dict)
    for (i, dat) in df_results.iterrows():
        for c in df_results.columns:
            if c == "names":
                continue
            if dat[c] == "TRUE":
                found_dict[c][dat["names"]] = True
            else:
                assert dat[c] == "FALSE"
                found_dict[c][dat["names"]] = False
    return found_dict


def quantify_distance_to_base(t, event_dict, root_distance, node2tips, use_median=False):
    """
    Quantify the distance of specific events in the tree to the root/base.
    
    Parameters:
    - t (Tree): A tree structure, typically from a library like ETE Toolkit.
    - event_dict (dict): A dictionary where keys are nodes and values are gains/losses at each node.
    - root_distance (dict): A dictionary that holds precalculated distances from each node to the tree root.
    - node2tips (dict): A dictionary with each node as key and its associated leaf nodes (tips) as values.
    - use_median (bool, optional): Whether to use the median instead of mean. Default is False.
    """

    dists_loss = []
    for node in t.traverse("postorder"):
        if (node.name == "") | (node.name == "root"):
            continue
        if event_dict[node.name] == "loss":
            dists_loss.append(dist_to_tips(node, root_distance, node2tips))
    if use_median:
        return np.median(dists_loss)
    return np.mean(dists_loss)


def prep_tree_for_ape(tree_str):
    """
    Prepare a Newick tree string for use with the ape package in R. Specifically it adds internal node names.

    Parameters:
    - tree_str: A Newick tree string.

    """
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


def infer_events(t, found_dict):
    """
    Infers gain/loss events in a phylogenetic tree based on ancestral state reconstruction.
        
    Parameters:
    - t (Tree): An ete3 phylogenetic tree.
    - found_dict (dict): A dictionary where keys are nodes in the tree and values 
                         represent the character state of each node.
    """
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


def calculate_distance_from_events(t, found_dict, root_distance, node2tips, use_median=False):
    """
    Computes the distance from phylogenetic events (e.g., gains or losses) to the tree base for all given randomizations.
    
    This function is designed to calculate empirical p-values by assessing the distribution of distances to the base for 
    observed phylogenetic events in comparison with a null model or randomized scenario. This helps in statistically 
    determining the significance of the observed distances.

    Parameters:
    - t: The phylogenetic tree structure, typically an object representing the tree topology and node relationships.
    - found_dict (dict): A dictionary containing phylogenetic events of interest (e.g., gains, losses).
    - root_distance (dict): A pre-calculated dictionary that maps each node in the tree to its distance from the root. 
                            From cache_distances()
    - node2tips (dict): A dictionary mapping nodes to the tips or leaves they are associated with.
    - use_median (bool, optional): If set to True, the function will use the median distance instead of the mean. 
                                   Default is False.
    """

    event_dict = infer_events(t, found_dict["reconstruction"])
    dist_reconstruction = quantify_distance_to_base(t, event_dict, root_distance, node2tips, use_median)

    dists = []
    n_randomizations = max([int(x.replace("randomization","")) for x in found_dict.keys() if "randomization" in x])
    for i in range(1, n_randomizations+1):
        event_dict = infer_events(t, found_dict["randomization{}".format(i)])
        dists.append(quantify_distance_to_base(t, event_dict, root_distance, node2tips, use_median))
    dists = np.array(dists)
    # Nans crop up in simulations where there are no losses
    dists = dists[~np.isnan(dists)]
    pval = np.sum(dists < dist_reconstruction) / len(dists)
    if pval == 0:
        pval = 1 / (len(dists) + 1)
    return dists, dist_reconstruction, np.sum(dists < dist_reconstruction), pval



def calculate_distances_bootstrap(df_ga3, infile_bootstrap, use_median=False):
    """
    Computes distances from phylogenetic events to the tree base for each bootstrap replicate.
    
    This function systematically calls 'calculate_distance_from_events' across all bootstrap replicates. By doing so, it 
    assesses the distribution of distances to the base for observed phylogenetic events in each replicate, aiding in 
    the evaluation of the robustness and reliability of the observed events across different phylogenetic reconstructions.

    Parameters:
    - df_ga3 (DataFrame): A pandas DataFrame containing the GA3 +/- data for the known species.
    - infile_bootstrap (str): Path to the file containing bootstrap replicates of the phylogenetic tree.
    - use_median (bool, optional): If set to True, the function will use the median distance instead of the mean. 
                                   Default is False.
    """
    results = []
    dists_all = []
    with open(infile_bootstrap,"r") as f:
        for line in tqdm(f,total=1000):
            line = line.rstrip("\n")
            t, tree_str = prep_tree_for_ape(line)
            t.prune([x for x in t.iter_leaves() if x.name in df_ga3["name"].values])
            # Drop tips without annotations
            
            df_results = R_tree_reconstruction_randomization(tree_str, df_ga3)
            found_dict = convert_r_results(df_results)
            
            root_distance = cache_distances(t)
            node2tips = t.get_cached_content()

            dists, dist_reconstruction, n_lower, pval = calculate_distance_from_events(t, found_dict, root_distance, node2tips, use_median)
            dists_all.append(dists)
            results.append([np.median(dists), dist_reconstruction, n_lower, pval])
    df_out = pd.DataFrame(results, columns=["random_median","data_median","n_lower","pval"])
    return df_out


def calculate_parameters(df_ga3, infile_tree, use_median=False):
    """
    Computes distances from phylogenetic events to the tree base for the main phylogenetic tree.

    This function leverages 'calculate_distance_from_events' to calculate the distances from events to the base
    of the main phylogenetic tree. It helps to quantify parameters or metrics associated with observed phylogenetic
    events in the main tree, which can be contrasted with those from bootstrap replicates to assess reliability and 
    significance.
    Parameters:
    - df_ga3 (DataFrame): A pandas DataFrame containing the GA3 +/- data for the known species.
    - infile_tree (str): Path to the file containing the main phylogenetic tree.
    - use_median (bool, optional): If set to True, the function will use the median distance instead of the mean. 
                                   Default is False.
    """
    with open(infile_tree) as f:
        tree_str = f.read().rstrip("\n")

    t, tree_str = prep_tree_for_ape(tree_str)
    t.prune([x for x in t.iter_leaves() if x.name in df_ga3["name"].values])

    df_results = R_tree_reconstruction_randomization(tree_str, df_ga3)
    found_dict = convert_r_results(df_results)
    
    root_distance = cache_distances(t)
    node2tips = t.get_cached_content()
    dists, dist_reconstruction, n_lower, pval = calculate_distance_from_events(t, found_dict, root_distance, node2tips, use_median)
    print("distance {}, median randomized distance {}, mean randomized distance {}, n_lower{}, pval {}".format(dist_reconstruction, np.median(dists), np.mean(dists), n_lower, pval))



# Adapted from: https://www.biostars.org/p/97409/
def cache_distances(tree):
    """
    Precalculate and cache distances of all nodes to the root.
    
    Parameters:
    - tree (Tree): The input tree structure.
    """

    node2rootdist = {tree:0}
    for node in tree.iter_descendants('preorder'):
        node2rootdist[node] = node.dist + node2rootdist[node.up]
    return node2rootdist


def dist_to_tips(node, root_distance, node2tips):
    """
    Compute the mean distance from a given node to its associated tips (leaf nodes).
    
    Parameters:
    - node (Node): The node of interest.
    - root_distance (dict): A dictionary holding distances from each node to the tree root.
    - node2tips (dict): A dictionary with each node as key and its associated leaf nodes (tips) as values.
    """
    return np.mean([root_distance[tip]-root_distance[node]+node.dist for tip in node2tips[node]])


def do_reconstruction_ape(infile_tree, df_ga3):
    """
    Perform ancestral state reconstruction using R's ape package via rpy2.
        
    Parameters:
    - infile_tree (str): The path to the input tree file.
    - df_ga3 (DataFrame): A pandas DataFrame containing the GA3 +/- data for the known species.
    """
    
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


def bootstrap_count_gains_losses(df_ga3, infile_bootstrap =  "Data/RAxML_bootstrap.fragilis_tree_markersv3_simplify_boot"):
    """
    Counts the number of gain and loss events in bootstrap phylogenetic trees based on ancestral state reconstruction.
    
    Parameters:
    - df_ga3 (DataFrame): Input data that contains GA3 +/- to be used in ancestral state reconstruction.
    - infile_bootstrap (str, optional): Path to the file containing bootstrap trees. 
    """
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
