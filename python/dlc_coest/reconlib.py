#
# Reconciliation library
# Adapted from dlcpar.reconlib
#

# python libs
import copy
import collections

# rasmus, compbio libs
from rasmus import treelib
from rasmus import util
from compbio import phylo

#=============================================================================
# tree factorization functions

def is_full_tree(tree, stree, recon, events):
    """Checks that the tree has all implied internal nodes AND no extra nodes

    Does NOT handle delay nodes
    """

    for node in tree:
        if events[node] == "gene":
            continue

        snode = recon[node]
        schildren = snode.children
        nschildren = len(schildren)
        if nschildren != 2 and nschildren != 0:
            raise Exception("Species tree must be binary")
        schildren2 = [recon[child] for child in node.children]

        if events[node] == "spec":
            if len(node.children) == 1:
                # speciation followed by loss
                if schildren[0] != schildren2[0] and \
                   schildren[1] != schildren2[0]:
                    print >>sys.stderr, "Reconciliation mismatch under speciation-loss node %s" % node.name
                    return False
            elif len(node.children) == 2:
                # speciation
                if set(schildren) != set(schildren2):
                    print >>sys.stderr, "Reconciliation mismatch under speciation node %s" % node.name
                    return False
            else:
                raise Exception("Cannot handle non-binary trees")

        elif events[node] == "dup":
            if len(node.children) == 1:
                # extra node
                print >>sys.stderr, "Single child under duplication node %s" % node.name
                return False
            elif len(node.children) == 2:
                # duplication
                if not (snode == schildren2[0] == schildren2[1]):
                    print >>sys.stderr, "Reconciliation mismatch under duplication node %s" % node.name
                    return False
            else:
                raise Exception("Cannot handle non-binary trees")

        else:
            raise Exception("Invalid event %s at node %s" % (events[node], node.name))

    return True


def add_spec_from_dup_nodes(node, tree, recon, events):
   """
   Relabel the current speciation node 'node' as a duplication.
   Insert new speciation nodes BELOW gene node 'node'.
   New nodes reconcile to same species node as 'node'.
   Modifies recon and events accordingly.
   """

   assert events[node] == "spec"
   snode = recon[node]
   events[node] = "dup"

   # insert new nodes into tree
   added = []
   for child in list(node.children):
       added.append(phylo.add_spec_node(child, snode, tree, recon, events))

   return added


def add_implied_spec_nodes(tree, stree, recon, events):
    """
    Add speciation nodes to tree that are implied but are not present because of gene losses.

    Extends phylo.add_implied_spec_nodes to handle non-MPR.
    Only guaranteed to work for binary trees.
    """

    added_spec = phylo.add_implied_spec_nodes(tree, stree, recon, events)

    added_dup = []
    for node in list(tree):
        schildren = [recon[child] for child in node.children]
        if len(schildren) > 1 and len(set(schildren)) == 1 and events[node] != "dup":
            added_dup.extend(add_spec_from_dup_nodes(node, tree, recon, events))

    assert is_full_tree(tree, stree, recon, events)

    return added_spec, added_dup


def add_delay_nodes(node, tree, recon, events):
    """
    Insert new delay nodes BELOW gene node 'node'.
    New nodes reconcile to same species node as 'node'.
    Modifies recon and events accordingly.

    TODO: same as add_spec_from_dup_nodes
    """

    assert events[node] == "spec"
    snode = recon[node]
    events[node] = "dup"    # move current node to be internal to species branch

    # insert new nodes into tree
    added_delay = []
    for child in list(node.children):
        added_delay.append(phylo.add_spec_node(child, snode, tree, recon, events))

    return added_delay


def add_implied_delay_nodes(tree, stree, recon, events):
    """
    Add nodes to tree after each speciation to allow
    for delay between coalescence and speciation
    """

    added = []
    for node in list(tree):
        if events[node] == "spec" and len(node.children) > 1:
            added.extend(add_delay_nodes(node, tree, recon, events))
    return added


def add_implied_nodes(tree, stree, recon, events, delay=True):
    """Wrapper to add implied speciation nodes and delay nodes

    If 'delay' is set, then add delay nodes as well.
    """

    added_spec, added_dup = add_implied_spec_nodes(tree, stree, recon, events)
    if delay:
        added_delay = add_implied_delay_nodes(tree, stree, recon, events)
    else:
        added_delay = []
    return added_spec, added_dup, added_delay


def remove_implied_nodes(tree, added, recon=None, events=None):
    """Wrapper to remove implied nodes from the tree"""
    for node in added:
        phylo.remove_spec_node(node, tree, recon, events)


def get_subtree(node, snode, recon, events, delay=True):
    """Returns the leaves of a duplication subtree

    Modified from dlcpar to handle coalescence only.
    """

    leaves = []

    def walk(node):
        if recon[node] != snode:
            return
        #if events[node] != "dup":
        #    leaves.append(node)
        if events[node] == "gene":
            leaves.append(node)
        elif events[node] == "spec":
            if not delay:
                # if delay nodes do not exist
                # go one past -- this should not be
                # used when determining coalescences
                leaves.append(node.children)
            else:
                leaves.append(node)
        else:
            for child in node.children:
                walk(child)
    walk(node)

    return leaves


def factor_tree(tree, stree, recon, events, delay = True):
    """Returns subtrees for each species branch

    Output is a dict with key = snode, val = (root(subtree), start, leaves(subtree)),
    where start is the node of the gene tree on which to recur
    (either the root of the subtree or its child).
    """

    # initialize
    subtrees = {}
    for snode in stree:
        subtrees[snode] = []

    # root is a dup
    if events[tree.root] != "spec":
        snode = recon[tree.root]
        subleaves = get_subtree(tree.root, snode, recon, events, delay = True)
        subtrees[snode].append((tree.root, tree.root, subleaves))

    # find subtrees
    for node in tree:
        if events[node] == "spec":
            snode = recon[node]
            for schild in snode.children:
                nodes2 = [x for x in node.children if recon[x] == schild]
                #print node.name, snode, nodes2
                if len(nodes2) > 0:
                    assert len(nodes2) == 1, (node, nodes2)
                    node2 = nodes2[0]
                    subleaves = get_subtree(node2, schild, recon, events)
                    subtrees[schild].append((node, node2, subleaves))
                else:
                    subtrees[schild].append((node, None, None))

    return subtrees

#=============================================================================
# labeled history functions

def enum_labeled_histories_subtree(tree, subtree):
    """
    This function enumerates all possible labeled histories for a subtree (a branch). The subtree is the value
    returned by factor_tree function.

    subtree is a list of subtrees returned by factor_tree
    """
    start = []
    leaves = []
    for (r, s, l) in subtree:
        if s is not None :
            if len(s.children) >= 2:
                start.append(s)
                leaves.extend(l)
    # print start
    # print leaves
    
    # recursive helper function to update current order (lst) and tochoose
    def helper(lst, tochoose):
        if len(tochoose) == 0:
            yield lst
        else:
            for child in tochoose:
                # print child
                next_lst = lst[:]
                next_lst.append(child)

                next_tochoose = tochoose[:]
                next_tochoose.remove(child)

                toadd = [node for node in child.children if not node in leaves]
                next_tochoose.extend(toadd)

                for local_order in helper(next_lst, next_tochoose):
                    yield local_order
     
    lst = []
    tochoose = start
    for order in helper(lst, tochoose):
        yield order


def enum_labeled_histories(tree):
    """
    This function enumerates all possible labeled histories for a tree topology.
    """

    # tree leaves
    leaves = tree.leaves()

    # recursive helper function to update current order (lst) and tochoose
    def helper(lst, tochoose):
        if len(tochoose) == 0:
            yield lst
        else:
            for child in tochoose:
                next_lst = lst[:]
                next_lst.append(child)

                next_tochoose = tochoose[:]
                next_tochoose.remove(child)

                toadd = [node for node in child.children if not node.is_leaf()]
                next_tochoose.extend(toadd)

                for local_order in helper(next_lst, next_tochoose):
                    yield local_order

    # starting case
    lst = [tree.root]
    tochoose = [child for child in tree.root.children if not child.is_leaf()]
    for order in helper(lst, tochoose):
        yield order
