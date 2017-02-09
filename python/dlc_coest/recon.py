# python libs
import sys
import copy

# dlc_coest libs
from __init__ import ALIGNMENT_SEED
import reconprob
import pllprob

# rasmus, compbio libs
from rasmus import util
from compbio import phylo, alignlib

# dlcoal libs
import dlcoal
from dlcoal import recon, coal

#=================================================================================================================
#=================================================================================================================
# Main Class

def dlc_coest_recon(alnfile, partfile, stree, gene2species,
                    popsizes, duprate, lossrate, subrate, beta,
                    pretime=None, premean=None,
                    niter=1,
                    nsearch_coal=1000, nsamples_coal=1000, nprescreen_coal=20, init_coal_tree=None,
                    nsearch_locus=1000, nsamples_locus=1000, nprescreen_locus=20, init_locus_tree=None,
                    threads=1, seed=ALIGNMENT_SEED, eps=2,
                    log=sys.stdout):

    reconer = DLC_Coest_Recon(alnfile, partfile, stree, gene2species,
                              popsizes, duprate, lossrate, subrate, beta,
                              pretime=pretime, premean=premean,
                              nsamples_coal=nsamples_coal, nprescreen_coal=nprescreen_coal,
                              nsamples_locus=nsamples_locus, nprescreen_locus=nprescreen_locus,
                              threads=threads, seed=seed, eps=eps,
                              init_coal_tree=init_coal_tree,
                              init_locus_tree=init_locus_tree,
                              log=log)
    maxrecon = reconer.recon(niter, nsearch_locus, nsearch_coal)

    return maxrecon.get_dict()


def get_init_coal_tree(alnfile):
    """
    Return an initial coal_tree from alignment by neighbor joining

    Works only with DNA sequences with alignment length that is multiple of 3.
    """
    aln = phylip.read_phylip_align(alnfile)
    distmat = alignlib.calc_four_fold_dist_matrix(aln)
    tree = phylo.neighborjoin(distmat, aln.keys())
    return tree

#============================================================================================================
#============================================================================================================
# Coestimation Recon Class

class DLC_Coest_Recon(object):

    def __init__(self, alnfile, partfile, stree, gene2species,
                 popsizes, duprate, lossrate, subrate, beta,
                 pretime=None, premean=None,
                 nsamples_coal=1000, nprescreen_coal=20, init_coal_tree=None,
                 nsamples_locus=1000, nprescreen_locus=20, init_locus_tree=None,
                 threads=1, seed=ALIGNMENT_SEED, eps=2,
                 name_internal="n", log=sys.stdout):

        self.alnfile = alnfile
        self.partfile = partfile
        self.stree = stree
        self.gene2species = gene2species
        self.popsizes = popsizes
        self.duprate = duprate
        self.lossrate = lossrate
        self.subrate = subrate
        self.beta = beta

        self.pretime = pretime
        self.premean = premean

        self.nsamples_coal = nsamples_coal
        self.nprescreen_coal = nprescreen_coal

        self.nsamples_locus = nsamples_locus
        self.nprescreen_locus = nprescreen_locus

        self.name_internal = name_internal
        self.log_stream = log

        self.threads = threads
        self.seed = seed
        self.eps = eps

        # init coal tree
        if init_coal_tree is None:
            init_coal_tree = get_init_coal_tree(self.alnfile)
        self.init_coal_tree = init_coal_tree.copy()
        self.coal_tree = init_coal_tree.copy()  # copy probably not necessarily

        # init locus tree
        if init_locus_tree is None:
            init_locus_tree = init_coal_tree.copy()
        self.init_locus_tree = init_locus_tree.copy()  # copy probably not necessarily
        self.locus_tree = init_locus_tree.copy()

        # note: the proposer should be externally set
        self.proposer = None


        # compute background parameters rates, frequencies, and alphas
        rates, freqs, alphas = pllprob.optimize_parameters(alnfile, partfile, init_coal_tree,
                                                           threads=threads, seed=seed, eps=eps)
        self.rates = rates
        self.freqs = freqs
        self.alphas = alphas


    def set_log(self, log):
        self.log_stream = log


    def recon(self, niter=1000, nsearch_locus=1000, nsearch_coal=1000):
        """Perform reconciliation"""

        util.tic("reconciling...")

        coal_tree = self.init_coal_tree.copy()
        locus_tree = self.init_coal_tree.copy()

        for i in xrange(niter):
            #========================================
            # fix coal_tree, estimate locus_tree
            util.tic("iter %d: estimating locus_tree" % i)

            self.stage = "locus_tree"
            locus_search = lambda ltree: DLCLocusTreeSearch(ltree, self.stree, self.gene2species,
                                                            self.duprate, self.lossrate,
                                                            nprescreen=self.nprescreen_locus)
            self.proposer = DLC_ReconProposer_LocusTree(
                coal_tree, self.stree, self.gene2species, search=locus_search)

            self.proposer.set_locus_tree(locus_tree.copy())
            maxrecon = self.recon_helper(nsearch_locus)

            util.toc()

            # update the reconciliation
            coal_tree = maxrecon.coal_tree.copy()
            locus_tree = maxrecon.locus_tree.copy()
            daughters = set([locus_tree.nodes[x.name] for x in maxrecon.daughters])

            #========================================
            # fix locus_tree, estimate coal_tree
            util.tic("iter %d: estimating coal_tree" % i)

            self.stage = "coal_tree"
            coal_search = lambda ctree: DLCCoalTreeSearch(ctree, locus_tree, daughters, self.popsizes,
                                                          nprescreen=self.nprescreen_coal)

            self.proposer = DLC_ReconProposer_CoalTree(
                self.stree, locus_tree, daughters, self.gene2species, search=coal_search)

            self.proposer.set_coal_tree(coal_tree.copy())
            maxrecon = self.recon_helper(nsearch_coal)

            # update the reconciliation
            coal_tree = maxrecon.coal_tree.copy()
            locus_tree = maxrecon.locus_tree.copy()
            daughters = set([locus_tree.nodes[x.name] for x in maxrecon.daughters])

            util.toc()

        util.toc()

        # return the best reconciliation
        return self.maxrecon


    def recon_helper(self, nsearch=1000):
        """Perform reconciliation"""

        self.maxp = -util.INF
        self.maxrecon = None
        proposal = self.proposer.init_proposal()
        init_proposal = proposal.copy()

        for i in xrange(nsearch):
            if i % 10 == 0:
                print "search", i

            # evaluate the probability of proposal
            util.tic("eval")
            p = self.eval_proposal(proposal)
            util.toc()

            # evaluate the search, then keep or discard the proposal
            util.tic("prop")
            self.eval_search(p, proposal)
            proposal = self.proposer.next_proposal() # set the next proposal
            util.toc()

        # all proposals bad, use initial proposal
        if not self.maxrecon:
            self.maxrecon = init_proposal

        # rename locus tree nodes
        dlcoal.rename_nodes(self.maxrecon.locus_tree, self.name_internal)  # how about coal_tree names?

        return self.maxrecon


    def next_proposal(self):
        """Returns next proposal"""
        self.proposal.next_proposal()


    def eval_proposal(self, proposal):
        """Compute probability of proposal"""

        ### util.tic("eval")
        # proposal.coal_tree.write()

        # compute recon probability
        info = {}

        # DEBUG
        counts = coal.count_lineages_per_branch(proposal.coal_tree,
                                                proposal.coal_recon,
                                                proposal.locus_tree)
        maxcount = max(x[0] for x in counts.values())

        # util.logger("max lineage count %d" % maxcount)
        # do not continue if there are too many deep coalescences
        if maxcount > 10:
            return -util.INF

        # estimate coal_tree
        if self.stage == "coal_tree":
             #util.tic("reconprob")
             p = reconprob.prob_gene_species_alignment_recon(self.alnfile, self.partfile, self.stree,
                                                             self.popsizes, self.duprate, self.lossrate, self.subrate, self.beta,
                                                             self.pretime, self.premean,
                                                             proposal.coal_tree, proposal.coal_recon, self.nsamples_coal,
                                                             proposal.locus_tree, proposal.locus_recon,  self.nsamples_locus,
                                                             proposal.daughters,
                                                             self.rates, self.freqs, self.alphas,
                                                             self.threads, self.seed, self.eps,
                                                             info=info)
             #util.toc()

             """
             p_joint = reconprob.prob_locus_gene_species_alignment_recon(self.alnfile, self.partfile, self.stree,
                                                                         self.popsizes, self.duprate, self.lossrate, self.subrate, self.beta,
                                                                         self.pretime, self.premean,
                                                                         proposal.coal_tree, proposal.coal_recon, self.nsamples_coal,
                                                                         proposal.locus_tree, proposal.locus_recon,  self.nsamples_locus,
                                                                         proposal.daughters,
                                                                         self.rates, self.freqs, self.alphas,
                                                                         self.threads, self.seed, self.eps,
                                                                         info=info)

             util.log("p = %f p_joint = %f [maxp = %f]" % (p, p_joint, self.maxp))
             print '(p,p_joint,maxp', (p, p_joint, self.maxp)
             """
        # estimate locus_tree
        else:
             p = dlcoal.prob_dlcoal_recon_topology(stree = self.stree, n = self.popsizes,
                                                   duprate = self.duprate, lossrate = self.lossrate,
                                                   pretime = self.pretime, premean = self.premean,
                                                   coal_tree = proposal.coal_tree, coal_recon = proposal.coal_recon,
                                                   locus_tree = proposal.locus_tree, locus_recon = proposal.locus_recon,
                                                   locus_events = proposal.locus_events,
                                                   daughters = proposal.daughters,
                                                   nsamples=self.nsamples_locus,
                                                   add_spec=False,
                                                   info=info)
             """
             p_joint = reconprob.prob_locus_gene_species_alignment_recon(self.alnfile, self.partfile, self.stree,
                                                                         self.popsizes, self.duprate, self.lossrate, self.subrate, self.beta,
                                                                         self.pretime, self.premean,
                                                                         proposal.coal_tree, proposal.coal_recon, self.nsamples_coal,
                                                                         proposal.locus_tree, proposal.locus_recon,  self.nsamples_locus,
                                                                         proposal.daughters,
                                                                         self.rates, self.freqs, self.alphas,
                                                                         self.threads, self.seed, self.eps,
                                                                         info=info)
             util.log("p = %f p_joint = %f [maxp = %f]" % (p, p_joint, self.maxp))
             print '(p,p_joint,maxp', (p, p_joint, self.maxp)
             """
        proposal.data = info
        #util.toc()

        return p


    def eval_search(self, p, proposal):
        """Evaluate a proposal for search"""

        self.log_proposal(proposal)

        # keep the proposal if it is better than maxp
        if p > self.maxp:
            #util.log("accept")

            self.maxp = p
            self.maxrecon = proposal.copy()

            # search with a new copy
            self.proposer.accept()
        else:
            #util.log("reject")
            self.proposer.reject()


    def log_proposal(self, proposal):
        self.log_stream.write(repr(proposal) + "\n")
        self.log_stream.flush()




#================================================================================================================
# ===============================================================================================================
# Proposer Class

class DLC_ReconProposer_LocusTree(dlcoal.recon.DLCoalReconProposer):

    """
    inherit the DLCoalReconPropose Class in dlcoal.recon;
    propose a locus tree when coal tree is fixed
    """

    def init_proposal(self):
        recon = super(DLC_ReconProposer_LocusTree, self).init_proposal()
        return Recon(self._coal_tree, recon.coal_recon,
                     recon.locus_tree, recon.locus_recon, recon.locus_events, recon.daughters)

    def next_proposal(self):
        recon = super(DLC_ReconProposer_LocusTree, self).next_proposal()
        return Recon(self._coal_tree, recon.coal_recon,
                     recon.locus_tree, recon.locus_recon, recon.locus_events, recon.daughters)

    def __recon_lca(self, locus_tree):
        recon = super(DLC_ReconProposer_LocusTree, self).__recon_lca(locus_tree)
        return Recon(self._coal_tree, recon.coal_recon,
                     recon.locus_tree, recon.locus_recon, recon.locus_events, recon.daughters)


class DLC_ReconProposer_CoalTree (object):

    def __init__(self, stree, locus_tree, daughters, gene2species,
                 search=phylo.TreeSearchNni,
                 num_coal_recons=1):
        self._stree = stree
        self._locus_tree = locus_tree
        self._daughters = daughters
        self._coal_search = search(None)

        # locus recon (static) -- propose LCA reconciliation
        self._locus_recon = phylo.reconcile(locus_tree, stree, gene2species)
        self._locus_events = phylo.label_events(locus_tree, self._locus_recon)

        # coal recon search
        self._num_coal_recons = num_coal_recons
        self._i_coal_recons = 1
        self._coal_recon_enum = None
        self._coal_recon_depth = 2
        self._accept_coal = False

        self._recon = None


    def set_coal_tree(self, coal_tree):
        self._coal_search.set_tree(coal_tree)

    def init_proposal(self):
        """Get first proposal"""

        if self._coal_search.get_tree() is None:
            self._coal_search.set_tree(self._locus_tree.copy())
        self._i_coal_recons = 1
        self._recon = self._recon_lca(self._coal_search.get_tree().copy())

        return self._recon


    def next_proposal(self):

        # if leaves <= 2, no need to propose another tree
        if len(self._coal_search.get_tree().leaves()) <= 2:
            return self._recon

        # always tree now because _num_coal_recon is set 0
        if self._i_coal_recons >= self._num_coal_recons:
            # propose new locus_tree

            # if locus_tree has not yet been accepted, then revert it
            if not self._accept_coal:
                self._coal_search.revert()

            # propose a new coal tree using _coal_search
            self._coal_search.propose()
            self._accept_coal = False
            self._i_coal_recons = 1

            # set coal_tree to be the proposed one
            coal_tree = self._coal_search.get_tree().copy()

            # make recon root optimal
            phylo.recon_root(coal_tree, self._locus_tree,
                             newCopy=False)
            dlcoal.rename_nodes(coal_tree)

            # propose remaining parts of dlcoal recon
            # reconciliation is given by lca
            self._recon = self._recon_lca(coal_tree)

        else:
            # modify coal_recon
            try:
                self._i_coal_recons += 1
                self._coal_recon_enum.next()
            except StopIteration:
                self._i_coal_recon = self._num_coal_recons
                return self.next_proposal()

        return self._recon


    def _recon_lca(self, coal_tree):
        # get coal tree, and LCA coal_recon
        coal_recon = phylo.reconcile(coal_tree,
                                     self._locus_tree, lambda x: x)

        # we do not explore the reconciliation space now
        self._coal_recon_enum = phylo.enum_recon(coal_tree,
                                                 self._locus_tree,
                                                 recon=coal_recon,
                                                 depth=self._coal_recon_depth)


        return Recon(coal_tree, coal_recon, self._locus_tree,
                     self._locus_recon, self._locus_events, self._daughters)


    def accept(self):
        self._accept_coal = True


    def reject(self):
        pass


#====================================================================================================
#====================================================================================================
# recon data Structure

class Recon (object):
    """
    The reconciliation datastructure for the DLCoal model
    (same as DLCoal.Recon except keep track of coal_tree)
    """

    def __init__(self, coal_tree, coal_recon, locus_tree, locus_recon, locus_events,
                 daughters, data=None):
        self.coal_tree = coal_tree
        self.coal_recon = coal_recon
        self.locus_tree = locus_tree
        self.locus_recon = locus_recon
        self.locus_events = locus_events
        self.daughters = daughters

        if data is None:
            self.data = {}
        else:
            self.data = data

    def copy(self):
        return Recon(self.coal_tree, self.coal_recon,
                     self.locus_tree, self.locus_recon, self.locus_events,
                     self.daughters, data=copy.deepcopy(self.data))

    def get_dict(self):
        return {"coal_tree": self.coal_tree,
                "coal_recon": self.coal_recon,
                "locus_tree": self.locus_tree,
                "locus_recon": self.locus_recon,
                "locus_events": self.locus_events,
                "daughters": self.daughters}

    def __repr__(self):
        return repr({"coal_tree": self.coal_tree.get_one_line_newick(
                         root_data=True),
                     "coal_recon": [(x.name, y.name) for x,y in
                                    self.coal_recon.iteritems()],
                     "locus_tree": self.locus_tree.get_one_line_newick(
                         root_data=True),
                     "locus_recon": [(x.name, y.name) for x,y in
                                     self.locus_recon.iteritems()],
                     "locus_events": [(x.name, y) for x,y in
                                      self.locus_events.iteritems()],
                     "daughters": [x.name for x in self.daughters],
                     "data": self.data})


#===================================================================================
#===================================================================================
# Tree Search

class DLCCoalTreeSearch (phylo.TreeSearch):
    """Adapted from dlcoal.DLCoalTreeSearch but to estimate coal_tree for fixed locus_tree"""

    def __init__(self, tree, locus_tree, daughters, popsizes,
                 tree_hash=None, nprescreen=20, weight=.2):

        phylo.TreeSearch.__init__(self, tree)

        self.locus_tree = locus_tree
        self.daughters = daughters
        self.popsizes = popsizes

        mix = phylo.TreeSearchMix(tree)
        mix.add_proposer(phylo.TreeSearchNni(tree), .4)
        mix.add_proposer(phylo.TreeSearchSpr(tree), .6)
        #self.search = phylo.TreeSearchUnique(tree, mix, tree_hash)

        prescreen = phylo.TreeSearchPrescreen(tree, mix,
                                              self.prescreen,
                                              poolsize=nprescreen)

        mix2 = phylo.TreeSearchMix(tree)
        mix2.add_proposer(prescreen, 1.0-weight)
        mix2.add_proposer(mix, weight)

        self.search = mix2

    def set_tree(self, tree):
        self.tree = tree
        self.search.set_tree(tree)

    def reset(self):
        self.search.reset()

    def propose(self):
        self.search.propose()
        return self.tree

    def revert(self):
        self.search.revert()
        return self.tree

    def prescreen(self, tree):
        # tree is coal tree, compute the associated lca reconciliation
        recon = phylo.reconcile(tree, self.locus_tree)

        # calculate the log probability of a reconciled coalescent tree (topology + reconciliation)
        # against the locus tree under the coalescent model
        return reconprob.prob_locus_coal_recon_topology(tree, recon,
                                                        self.locus_tree, self.popsizes,
                                                        self.daughters)


DLCLocusTreeSearch = dlcoal.recon.DLCoalTreeSearch
