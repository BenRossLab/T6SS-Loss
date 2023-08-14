
import re
from collections import Counter, defaultdict
from ete3 import Tree

class SpeciesToGenusTree():
    def __init__(self, tree_str, genera_oi, fraction_genus_cutoff=0.75, genus_purity=0.75):
        self.genera_oi = genera_oi
        self.load_tree(tree_str)
        self.record_genus_dict()
        nodes_use = self.determine_nodes_to_use(fraction_genus_cutoff, genus_purity)
        self.collapse_nodes_into_genera(nodes_use)
        self.remove_nongenera_nodes()
    def load_tree(self, tree_str):
        # Drop species not in the genera of interest
        self.t = Tree(tree_str)
        self.leafs_all = [x for x in self.t.get_leaf_names()]

        self.t.prune([x for x in self.leafs_all if x.split("|")[-2] in self.genera_oi])

        # Kill all "_unclassified" bugs
        self.leafs_all = [x for x in self.t.get_leaf_names()]
        self.t.prune([x for x in self.leafs_all if "_unclassified" not in x.split("|")[-2] ])
        self.leafs_all = [x for x in self.t.get_leaf_names()]
        
        self.name_nodes()

    def screen_purity(self, genera_descendants, genus_purity):
        top_genus = Counter(genera_descendants).most_common(1)[0]
        if top_genus[1] / len(genera_descendants)  < genus_purity:
            return True
        return False

    def name_nodes(self):
        # Name the internal nodes
        i = 1
        for node in self.t.traverse("postorder"):
            if node.is_leaf():
                continue
            node.name = "n{}".format(i)
            i += 1

    def record_genus_dict(self):
        # Record where the genus names of the nodes
        self.species_dict = defaultdict(list)
        self.genus_dict = {}
        for species in self.leafs_all:
            genera = species.split("|")[-2]
            self.species_dict[genera].append(species)
            self.genus_dict[species] = genera

    def determine_nodes_to_use(self, fraction_genus_cutoff=0.75, genus_purity=0.75):
        n_descendants_best = defaultdict(int)
        node_use = {}

        for node in self.t.traverse("postorder"):
            if node.is_leaf():
                genera_descendants = [self.genus_dict[node.name]]
            else:
                genera_descendants = [self.genus_dict[x.name] for x in node.get_leaves()]
            
            if self.screen_purity(genera_descendants, genus_purity):
                continue
           
            genus, _ = Counter(genera_descendants).most_common(1)[0]
            # Number of species in the genus
            n_species_full = len(self.species_dict[genus])
            # Number of species in the clade
            n_descendants = len(genera_descendants)
            if n_descendants / n_species_full >= fraction_genus_cutoff:
                if n_descendants_best[genus] < n_descendants:
                    n_descendants_best[genus] = n_descendants
                node_use[genus] = node.name
        print("Able to keep {} out of {} genera".format(len(node_use), len(self.genera_oi)))
        return node_use

    def collapse_nodes_into_genera(self, node_use):
        # Collapse all species into a genus node
        for key in node_use:
            nodes_search = self.t.search_nodes(name=node_use[key])
            if len(nodes_search) == 0:
                continue
            node = nodes_search[0]
            # Don't delete this if there are genus nodes in the children
            if any([bool(re.search("g__[A-Za-z_]+$",leaf.name)) for leaf in node.get_leaves()]):
                continue
            node.name = key
            if node.is_leaf():
                continue
            for child in node.get_children():
                node.remove_child(child)
         
        self.leafs_all = [x for x in self.t.get_leaf_names()]

    def remove_nongenera_nodes(self):
        # Drop extra species nodes that were not accounted for
        # For reasons I don't quite understand, we need to run this code repeatedly
        while True:
            n_deleted = 0
            for node in self.t.iter_leaves():
                if "s__" not in node.name:
                    continue

                # Go up the tree, find the last node you can knock off without removing a genera node
                while True:
                    parent = node.up
                    if any([bool(re.search("g__[A-Za-z]+$",x.name)) for x in parent.iter_leaves()]):
                        break
                    node = parent
                n_deleted += 1
                node.delete()
            if n_deleted == 0:
                break

    def write_tree(self):
        return self.t.write()

