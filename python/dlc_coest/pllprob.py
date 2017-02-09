"""
Wrapper around pllpy to compute tree likelihoods.
"""

# python libs
import pllpy
import os
import tempfile

# dlc_coest libs
from __init__ import ALIGNMENT_SEED

# rasmus, compbio libs
from rasmus import treelib
from compbio import phylip

#====================================
# optimization functions

def optimize_parameters(alnfile, partfile, coal_tree,
                        threads=1, seed=ALIGNMENT_SEED, eps=0.1):
    """
    The function takes in alignment file, partitions file, coal_tree,
    and return the rates, freqs, alphas after optimization. These 
    parameters are used when the alignment probability is calculated.
    """
    rates = []
    freqs = []
    alphas = []

    # convert coal_tree to filename
    tree_temp = tempfile.NamedTemporaryFile(delete=False)
    tree = treelib.unroot(coal_tree, newCopy=True) # unrooted tree required for ML
    tree.write(tree_temp, oneline=True)
    tree_temp.close()
    tree_filename = tree_temp.name

    # initialize pll instance
    pll = pllpy.pll(alnfile, partfile,
                    tree_filename, threads, seed)

    # set likelihood convergence
    pll.set_epsilon(eps)

    # optimize rates, freqs, alphas, and branches
    pll.optimise(True, True, True, True)

    # store optimal parameters
    for i in range(pll.get_number_of_partitions()):
        rates.append(pll.get_rates_vector(i))
        freqs.append(pll.get_frequencies_vector(i))
        alphas.append(pll.get_alpha(i))

    os.remove(tree_filename)
    return rates, freqs, alphas

#===============================================================================
# probability functions

def prob_alignment_nooptimize(alnfile, partfile, coal_tree,
                              rates, freqs, alphas,
                              threads=1, seed=ALIGNMENT_SEED, eps=0.1):

    """
    This function implements the pll function. It computes the log likelihood
    of alignment data given the coal_tree without optimizing parameters. 
    Mathematically, it computes: P(A | T^G, t^G).
        
    rates, freqs, alphas  -- parameters in pll computation
    alnfile               -- alignment file
    partfile              -- partition file
    coal_tree             -- coalescent tree
 
    """

    # convert coal_tree to filename
    tree_temp = tempfile.NamedTemporaryFile(delete=False)
    tree = treelib.unroot(coal_tree, newCopy=True) # unrooted tree required for ML
    tree.write(tree_temp, oneline=True)
    tree_temp.close()
    tree_filename = tree_temp.name

    # initialize pll instance
    pll = pllpy.pll(alnfile, partfile,
                    tree_filename, threads, seed)

    # initialize pll with previously optimized parameters
    for i in range(pll.get_number_of_partitions()):
        pll.set_alpha(alphas[i], i, True)
        pll.set_frequencies(freqs[i], i, True)
        if pll.is_dna(i):
            pll.set_rates(rates[i], i, True)

    # set likelihood convergence
    pll.set_epsilon(eps)

    # do not optimize any of the parameters
    # pll.optimise(False, False, False, False)

    # get (log) likelihood
    prob = pll.get_likelihood()

    os.remove(tree_filename)
    return prob


def prob_alignment(alnfile, partfile, coal_tree,
                   threads=1, seed=int("0xDEADBEEF", 16), eps=0.1,
                   opt_branches=False):
    """
    This function implements the pll function. It optimize the alpha, rates
    and frequencies, and use these parameters to compute the log likelihood
    of alignment data given the coal_tree. This function is not used because 
    it is computationally inefficient to optimize the parameters.
	 
    """


    # convert coal_tree to filename
    tree_temp = tempfile.NamedTemporaryFile(delete=False)
    tree = treelib.unroot(coal_tree, newCopy=True) # unrooted tree required for ML
    tree.write(tree_temp, oneline=True)
    tree_temp.close()
    tree_filename = tree_temp.name

    # initialize pll instance
    pll = pllpy.pll(alnfile, partfile,
                    tree_filename, threads, seed)

    # tell pll to optimize all model parameters
    for i in range(pll.get_number_of_partitions()):
        pll.set_optimisable_alpha(i, True)
        pll.set_optimisable_frequencies(i, True)
        if pll.is_dna(i):
            pll.set_optimisable_rates(i, True)

    # set likelihood convergence
    pll.set_epsilon(eps)

    # optimize the model
    pll.optimise(True, True, True, opt_branches) # rates, freqs, alphas, branches

    # get (log) likelihood
    prob = pll.get_likelihood()

    os.remove(tree_filename)
    return prob
