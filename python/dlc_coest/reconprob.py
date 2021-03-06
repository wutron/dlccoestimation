# python libs
from __future__ import division
from math import exp,log,sqrt

import random
import numpy

# dlc_coest libs
from __init__ import ALIGNMENT_SEED
import reconlib
import pllprob

# rasmus, compbio libs
from rasmus import treelib, util, stats
from compbio import birthdeath, phylo, phylip

# dlcoal libs
import dlcoal
from dlcoal import duploss, coal

# root finder
try:
    from scipy.optimize import brentq
except ImportError:
    def brentq(f, a, b, disp=False):
        return stats.bisect_root(f, a, b)

#=========================================================================================================
# coal times sampling

def sample_coal_times_one_branch(bottom_num, top_num, popsizes, T):
    """
    Return a list of times of coalescence of 'k' to 'k-1', 'k-1' to 'k-2', ..., 'l+1' to 'l', and pretime
    from k = bottom_num to l = top_num.
    """

    num_coal = bottom_num - top_num
    times = [0 for _ in xrange(num_coal + 1)]

    if num_coal == 0 : # no times to sample
        times[0] = T

    elif T == util.INF:
        t = coal.sample_coal(bottom_num, popsizes) # time is unbounded at root
        times[0] = t
        times[1] = 0.0
    else:
        k = bottom_num
        remain = T
        for i in xrange(num_coal):
            # sample
            #print k, top_num, remain, popsizes
            t = sample_coal_cond_counts(k, top_num, remain, popsizes) # sample bounded coalescent times
            times[i] = t

            # update
            k -= 1
            remain -= t

        times[num_coal] = remain
    return times

def sample_coal_times_topology(coal_tree, coal_recon, locus_tree, popsizes) :
    """
    Sample the coalescent times for a topology by doing in a two-step way.
    Sample the labeled history consistent with the topology and then sample branch length
    using the labeled history. This sampling process captures the conditional probability:
    P(t^G | T^G, R^G, T^L, t^L, N^L)  

    coal_tree    -- coalescent tree
    coal_recon   -- reconciliation of coalescent tree to locus tree
    locus_tree   -- locus tree
    popsizes     -- population sizes of the locus tree

    """
    times = {}

    coal_events = phylo.label_events(coal_tree, coal_recon)

    lineages = coal.count_lineages_per_branch(coal_tree, coal_recon, locus_tree)

    added_spec, added_dup, added_delay = reconlib.add_implied_nodes(coal_tree, locus_tree, coal_recon, coal_events, delay=True)
    added = added_spec + added_dup + added_delay

    subtrees = reconlib.factor_tree(coal_tree, locus_tree, coal_recon, coal_events)

    for coal_node in coal_tree.preorder():
        coal_node.dist = 0.0

    for lnode in locus_tree.preorder():
        subtree = subtrees[lnode]

        # enumerate all labeled histories
        all_labeled_histories_subtree = list(reconlib.enum_labeled_histories_subtree(coal_tree, subtree))

        # sample a labeled history
        sample_labeled_history = random.sample(all_labeled_histories_subtree, 1)[0]

        bottom_num, top_num = lineages[lnode]
        # print lineages[lnode]
        if lnode.parent:
            T = lnode.dist
        else:
            T = util.INF
        popsizes = popsizes
        # print T

        # sample coalescent times for a branch
        try:
            coal_times_subtree = sample_coal_times_one_branch(bottom_num, top_num, popsizes, T)
        except:
            # clean up coal tree nodes
            reconlib.remove_implied_nodes(coal_tree, added, coal_recon, coal_events)
            raise

        set_coal_tree_subtree(sample_labeled_history, coal_times_subtree, coal_tree, subtree)

    # remove the implied nodes
    reconlib.remove_implied_nodes(coal_tree, added, coal_recon, coal_events)


def set_coal_tree_subtree(labeled_history, coal_times, coal_tree, subtree):
    """
    Set the coal_tree distances according to the labeled history and the sampled coal times.
    """
    leaves = [] # exract all leavese in subtree
    for (r, s, l) in subtree:
        if s is not None :
            leaves.extend(l)

    n = len(labeled_history)
    m = len(coal_times)
    assert n == m - 1, "n = %d, m = %d, LH = %s, t = %s" % (n, m, str(labeled_history), str(coal_times))

    if labeled_history == []:
        for node in leaves:
            node.dist = coal_times[0]
    else:
        nodes = leaves
        for i in xrange(n + 1):
            t = coal_times[i]
            for node in nodes:
                node.dist += t

            if i == n:
                break

            else:
                for child in labeled_history[n - 1 - i].children:
                    nodes.remove(child)
                nodes.append(labeled_history[n - 1 - i])


#=======================================================================================================
# probability calculation

def prob_locus_gene_species_alignment_recon(alnfile, partfile, stree,
                                           popsizes, duprate, lossrate, subrate, beta,
                                            pretime, premean,
                                            coal_tree, coal_recon, nsamples_coal,
                                            locus_tree, locus_recon,  nsamples_locus, daughters,
                                            rates, freqs, alphas,
                                            threads=1, seed=ALIGNMENT_SEED, eps=0.1,
                                            info=None):
    """
    (Log) probability of the joint probability of locus_tree, locus_recon, coal_tree,
    coal_recon, daughters and alignment. Mathematically, it computes: 
    
    P(T^G, T^L, R^G, R^L, delta^L, A | S, theta) = P(delta^L | T^L, R^L, S) + P(T^L, R^L | S, theta^S) + 
    int int P(t^L | T^L, R^L, S, theta) * P(T^G, R^G, t^G | t^L, T^L, daughters, R^L, theta) * P(A | T^G, t^G) dt^L dt^G

    alnfile           -- alignment file
    partfile	      -- partition file
    stree	      -- species tree
    popsizes          -- population sizes in species tree
    duprate           -- duplication rate
    lossrate          -- loss rate
    subrate           -- substitution rate
    beta              -- regularization parameter
    pretime           -- starting time before species tree
    premean           -- mean starting time before species tree

    coal_tree         -- coalescent tree
    coal_recon        -- reconciliation of coalescent tree to locus tree
    nsamples_coal     -- number of times to sample coal times t^G
    locus_tree        -- locus tree (has dup-loss)
    locus_recon       -- reconciliation of locus tree to species tree
    nsamples_locus    -- number of times to sample the locus tree times t^L
    daughters         -- daughter nodes
    
    rates, freqs, alphas  -- optimization parameters  

    Note: Adapted from dlcoal.prob_dlcoal_recon_topology(...) [in __init.py]
    """

    # duploss proability: P(T^L, R^L | S, theta)
    locus_events = phylo.label_events(locus_tree, locus_recon)
    dl_prob = duploss.prob_dup_loss(locus_tree, stree, locus_recon, locus_events,
                                      duprate, lossrate)

    # daughters probability: P(daughters | T^L, R^L, S)
    dups = phylo.count_dup(locus_tree, locus_events)
    daughter_prob = dups*log(.5)

    # double integral
    double_integral = prob_gene_species_alignment_recon(alnfile, partfile, stree,
                                                        popsizes, duprate, lossrate, subrate, beta,
                                                        pretime, premean,
                                                        coal_tree, coal_recon, nsamples_coal,
                                                        locus_tree, locus_recon,  nsamples_locus, daughters,
                                                        rates, freqs, alphas,
                                                        threads=1, seed=ALIGNMENT_SEED, eps=0.1,
                                                        info=None)

    return dl_prob+daughter_prob+double_integral


def prob_gene_species_alignment_recon(alnfile, partfile, stree,
                                      popsizes, duprate, lossrate, subrate, beta,
                                      pretime, premean,
                                      coal_tree, coal_recon, nsamples_coal,
                                      locus_tree, locus_recon,  nsamples_locus, daughters,
                                      rates, freqs, alphas,
                                      threads=1, seed=ALIGNMENT_SEED, eps=0.1,
                                      info=None):
    """
    Evaluate terms that depend on T^G and R^G.

    That is, fix T^L, R^L, and daughters and evaluate the double integral:
    int int P(t^L | T^L, R^L, S, theta) * P(T^G, R^G, t^G | t^L, T^L, daughters, R^L, theta) * P(A | T^G, t^G) dt^L dt^G

    This is the probability we used in the searching process. 

    alnfile           -- alignment file
    partfile          -- partition file
    stree             -- species tree
    popsizes          -- population sizes in species tree
    duprate           -- duplication rate
    lossrate          -- loss rate
    subrate           -- substitution rate
    beta              -- regularization parameter
    pretime           -- starting time before species tree
    premean           -- mean starting time before species tree

    coal_tree         -- coalescent tree
    coal_recon        -- reconciliation of coalescent tree to locus tree
    nsamples_coal     -- number of times to sample coal times t^G
    locus_tree        -- locus tree (has dup-loss)
    locus_recon       -- reconciliation of locus tree to species tree
    nsamples_locus    -- number of times to sample the locus tree times t^L
    daughters         -- daughter nodes
    
    rates, freqs, alphas  -- optimization parameters  

    """

    locus_events = phylo.label_events(locus_tree, locus_recon)

    # optimize the parameters
    # util.tic("optimize parameter")
    # rates, freqs, alphas = pllprob.optimize_parameters(alnfile, partfile, coal_tree,
    #                                                   threads=threads, seed=seed, eps=eps)
    # util.toc()
    # double integral
    double_integral_list = []
    double_integral = 0.0
    util.tic("recon prob")
    for i in xrange(nsamples_locus):

        # sample t^L, the unit should be in myr
        #util.tic("topo prob")
        locus_times = duploss.sample_dup_times(locus_tree, stree, locus_recon,
                                               duprate, lossrate, pretime, premean, events=locus_events)
        treelib.set_dists_from_timestamps(locus_tree, locus_times)

        # calculate P(T^G, R^G | T^L, t^L, daughters, theta)
        topology_prob = prob_locus_coal_recon_topology(coal_tree, coal_recon,
                                                       locus_tree, popsizes, daughters)
        #util.toc()
        # for a fixed t^L, compute coal_prob
        # sample t^G for topology and compute the probabililty of observing the alignment using MonteCarlo integration
        coal_prob = 0.0
        alignment_prob_MonteCarlo = 0.0
        alignment_prob_list = []

        # check probability of lineage counts for this locus tree
        zero_lineage_prob = False

        #util.tic("set times")
        for lnode in locus_tree:
            lineages = coal.count_lineages_per_branch(coal_tree, coal_recon, locus_tree)
            bottom_num, top_num = lineages[lnode]
            if lnode.parent:
                T = lnode.dist
            else:
                T = util.INF

            popsizes = popsizes
            lineage_prob = prob_coal_counts(bottom_num, top_num, T, popsizes)

            # set zero_lineage_prob = TRUE if one lineage returns zero probability
            if (lineage_prob == 0.0):
                zero_lineage_prob = True

        #util.toc()
        # if lineage_prob is zero, coal_prob is zero
        if zero_lineage_prob:
            coal_prob = -float("inf")

        # otherwise, we calculate the coal_prob
        else:
            for j in xrange(nsamples_coal):

                # sample coal times and set the coal_tree accordingly
                # locus tree branch lengths are in myr
                # make sure the input popsizes are scaled to fit the time unit (typically myr)

                try:
                    sample_coal_times_topology(coal_tree, coal_recon, locus_tree, popsizes)
                except (ZeroDivisionError, ValueError):
                    # bad sample
                    util.log("bad sample")
                    alignment_prob = -util.INF
                    continue

                #===============================================================================
                # (log) probability of observing the alignment
                #util.tic("alignment probability")

                # convert branch lengths from myr to sub/site
                for node in coal_tree:
                    node.dist *= subrate

                #util.tic("alignment prob")
                # set a regularization parameter beta
                print beta
                alignment_prob = beta * prob_alignment(alnfile, partfile, coal_tree,
                                                       rates, freqs, alphas,
                                                       threads=threads, seed=seed, eps=eps)
                #util.toc()
                ### util.log("p = %.6f" % alignment_prob)
                #util.toc()

                #===============================================================================
                ### util.log("   log p = %.6g" % alignment_prob)
                ### util.log("   p = %.6g" % exp(alignment_prob))
                alignment_prob_list.append(alignment_prob)

            ### util.log("p = %f" % alignment_prob_MonteCarlo)

            # log_sum_exp function exponentiate the log probability of observing alignment,
            # add them up, and take log again
            if len(alignment_prob_list) == 0:
                # all bad samples
                alignment_prob_MonteCarlo = -util.INF
            else:
                alignment_prob_MonteCarlo = log_sum_exp(alignment_prob_list) - log(nsamples_coal)

            # P(T^G, R^G | T^L, t^L, daughters, theta) * $ P(t^G | ~) * P(A | T^G,t^G) dtG
            # coal_prob is a log probability
            coal_prob += topology_prob + alignment_prob_MonteCarlo

            # add coal probability to a list for further processing
        double_integral_list.append(coal_prob)

        # log_sum_exp function exponentiate the log probability of observing alignment,
        # add them up, and take log again
        double_integral = log_sum_exp(double_integral_list) - log(nsamples_locus)

        # logging info
        if info is not None:
            info["topology_prob"] = topology_prob              # one sample of t^L
            info["alignment_prob"] = alignment_prob_MonteCarlo # one sample of t^L, averaged over t^G
            info["coal_prob"] = double_integral
    util.toc()
    return double_integral



prob_locus_coal_recon_topology = dlcoal.prob_locus_coal_recon_topology



def prob_alignment(alnfile, partfile, coal_tree,
                   rates, freqs, alphas,
                   threads=1, seed=ALIGNMENT_SEED, eps=0.1):
    """
    Return the probability of seeing the alignment given the coal_tree with branch lengths,
    e.g. P(A | T^G, t^G).
    """

    prob = pllprob.prob_alignment_nooptimize(alnfile, partfile, coal_tree,
                                             rates, freqs, alphas,
                                             threads=threads, seed=seed, eps=eps)
    return prob

#===========================================================================================================
# Some helper functions

def log_sum_exp(num_list):
    """
    Given a list of numbers, take exponentials of them, sum these exponentials up,
    and take natural log of the sum.
    log(exp(x1)+exp(x2)+ ... + exp(xn)) = m + log(exp(x1 - m)+exp(x2 - m)+ ... + exp(xn - m)),
    where m is max(x1, x2, ..., xn)
    """
    m = max(num_list)
    summation = 0
    for i in xrange(len(num_list)):
        summation += exp(num_list[i]-m)
    return log(summation) + m

prob_coal_counts = coal.prob_coal_counts

sample_coal_cond_counts = coal.sample_coal_cond_counts





if __name__ == "__main__":
    #========================================
    # test cases for prob_locus_gene_species_alignment_recon
    # sim-flies, N = 1e6, g = 0.1, R = 1x, L = 100bp, mu = 5e-9

    import os
    import numpy
    from compbio import fasta
    import dlcoal
    import StringIO

    path = "/home/muddcs15/research/work/coestimation/"
    stree = treelib.read_tree(os.path.join(path, "simulation/config/flies.stree"))

    prob_raxml = []
    prob_treefix = []
    prob_dlca = []

    for i in range(100, 200):

        # read raxml recon
        coal_tree_raxml, extra_raxml = dlcoal.read_dlcoal_recon(os.path.join(
                         path, "simulation/data/1000/5e-9/1e6-1x/",str(i),str(i)+".raxml.dlcoal"), stree)
        locus_tree_raxml = extra_raxml["locus_tree"]
        locus_recon_raxml = extra_raxml["locus_recon"]
        coal_recon_raxml = extra_raxml["coal_recon"]
        daughters_raxml = extra_raxml["daughters"]

        # read treefix recon
        coal_tree_treefix, extra_treefix = dlcoal.read_dlcoal_recon(os.path.join(
                        path, "simulation/data/1000/5e-9/1e6-1x/",str(i),str(i)+".raxml.treefix.dlcoal"), stree)
        locus_tree_treefix = extra_treefix["locus_tree"]
        locus_recon_treefix = extra_treefix["locus_recon"]
        coal_recon_treefix = extra_treefix["coal_recon"]
        daughters_treefix = extra_treefix["daughters"]

        # read dlca recon
        coal_tree_dlca, extra_dlca = dlcoal.read_dlcoal_recon(os.path.join(path,
                                    "dlc_coestimation/data/1000/5e-9/1e6-1x/",str(i),str(i)+".dlca"), stree)
        locus_tree_dlca = extra_dlca["locus_tree"]
        locus_recon_dlca = extra_dlca["locus_recon"]
        coal_recon_dlca = extra_dlca["coal_recon"]
        daughters_dlca = extra_dlca["daughters"]


        # read alignment
        alnfile = os.path.join(path, "simulation/data/1000/5e-9/1e6-1x/",str(i),str(i)+".coal.fasta.phylip")
        partfile = os.path.join(path, "simulation/data/1000/5e-9/1e6-1x/0/0.coal.partitions")

        # modify coal_tree and coal_recon
        #subtree, newpos = phylo.propose_random_spr(coal_tree)
        #phylo.perform_spr(coal_tree, subtree, newpos)
        #coal_recon = phylo.reconcile(coal_tree, locus_tree)

        # species tree parameters
        duprate = 0.12  # events/lineages/myr
        lossrate = 0.12 # events/lineages/myr
        gentime = 0.1   # yr / gen
        popsizes = 2 * 1e6 * gentime / 1e6  # "normalized popsize" = 2 (diploid) * Ne * yr/gen * myr/yr (first 1e6 is the pop size)
        subrate = 5e-9 / gentime * 1e6 # sub/site/myr = sub/site/gen * gen/yr * yr/myr\


        rates, freqs, alphas = pllprob.optimize_parameters(alnfile, partfile, coal_tree_treefix,
                                                           threads=1, seed=ALIGNMENT_SEED, eps=1)

        nsamples_coal = 1
        nsamples_locus = 1

        times = treelib.get_tree_timestamps(stree)
        pretime = None
        premean = 0.5*times[stree.root]

        p_raxml = prob_locus_gene_species_alignment_recon(alnfile, partfile, stree,
                                                    popsizes, duprate, lossrate, subrate,
                                                    pretime, premean,
                                                    coal_tree_raxml, coal_recon_raxml, nsamples_coal,
                                                    locus_tree_raxml, locus_recon_raxml, nsamples_locus, daughters_raxml,
                                                    rates, freqs, alphas,
                                                    threads=1, seed=ALIGNMENT_SEED, eps=1,
                                                    info=None)
        prob_raxml.append(p_raxml)

        p_treefix = prob_locus_gene_species_alignment_recon(alnfile, partfile, stree,
                                                    popsizes, duprate, lossrate, subrate,
                                                    pretime, premean,
                                                    coal_tree_treefix, coal_recon_treefix, nsamples_coal,
                                                    locus_tree_treefix, locus_recon_treefix, nsamples_locus, daughters_treefix,
                                                    rates, freqs, alphas,
                                                    threads=1, seed=ALIGNMENT_SEED, eps=1,
                                                    info=None)
        prob_treefix.append(p_treefix)

        p_dlca = prob_locus_gene_species_alignment_recon(alnfile, partfile, stree,
                                                    popsizes, duprate, lossrate, subrate,
                                                    pretime, premean,
                                                    coal_tree_dlca, coal_recon_dlca, nsamples_coal,
                                                    locus_tree_dlca, locus_recon_dlca, nsamples_locus, daughters_dlca,
                                                    rates, freqs, alphas,
                                                    threads=1, seed=ALIGNMENT_SEED, eps=1,
                                                    info=None)

        prob_dlca.append(p_dlca)

    print numpy.mean(prob_raxml),numpy.std(prob_raxml)/sqrt(len(prob_raxml)),numpy.mean(prob_dlca),numpy.std(prob_dlca)/sqrt(len(prob_dlca)), numpy.mean(prob_treefix),numpy.std(prob_treefix)/sqrt(len(prob_treefix))

    # look at difference between treefix and dlca
    # diff = [a_i - b_i for a_i, b_i in zip(prob_treefix, prob_dlca)]
    # numpy.mean(diff)
    # 771.55121053161417
    # numpy.std(diff)/sqrt(len(diff))
    # 118.78254408292744
    # diff = np.array(diff)
    # len(diff[diff>=0])
    # 90

