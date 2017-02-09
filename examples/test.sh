#!/usr/bin/bash
#
# This is an example of how to use dlc_coestimation to reconcile gene families.
#
#
# Don't execute this script all at once.  Instead try copying and pasting
# each command to the command line one at a time in order to learn how the
# commands work.
#

#=============================================================================
# setup/install

# Make sure tools are compiled and installed before running the commands in
# this tutorial.  See INSTALL.txt for more information.

# Or you can run from the source directory:

cd ..
python setup.py install

cd examples
export PATH=$PATH:../bin
export PYTHONPATH=$PYTHONPATH:../python


#=============================================================================
# reconcile gene tree using DLCpar

# show help information
dlc_coest -h

# Usage: dlc_coest [options] <gene tree> ...
#
# Options:
#   Input/Output:
#     -A <alignment>, --fasta.phylip=<alignment>
#                       alignment file in phylip format 
#     -p <partition>, --partitions=<partition>
#                       partition file 
#     -s <species tree>, --stree=<species tree>
#                         species tree file in newick format
#     -S <species map>, --smap=<species map>
#                         gene to species map
#
#   Model Parameters
#     -n <population sizes>, --popsizes=<population sizes>
#                              population sizes
#     -D <duplication rate>, --duprate=<duplication rate>
#                              duplication rate in dups/gene/myr
#     -L <loss rate>,        --lossrate=<loss rate>
#                              loss rate in loss/gene/myr   
#     -g <generation time>,  --gentime=<generation time>
#                              generation time in years
#     -subrate <substitution rate>, --subrate=<substitution rate>
#                                   substitution rate in sub/site/myr           
#
#   Search Parameters
#     -i <number of iterations>, --iter=<number of iterations>
#                                 number of search iterations
#     -nsearch_coal <number of coal search>, --nsearch_coal=<number of coal search>
#                                          number of search iterations for coal tree 
#     -nsamples_coal <number of coal samples>, --nsamples_coal=<number of coal samples
#                                             number of samples when integrate over coal tree
#     -nprescreen_coal <number of coal prescreen>, --nprescreen_coal=<number of coal prescreen>
#                                          number of prescreening when estimating coal tree
#     -init_coal_tree <initial coal tree>, --init_coal_tree=<initial coal tree
#                                          initial coal tree to start search   
#     -nsearch_locus <number of locus search>, --nsearch_locus=<number of locus search>
#                              number of search iterations for the locus tree 
#     -nsamples_locus <number of locus samples>, --nsamples_locus=<number of locus samples
#                                             number of samples when integrate over locus tree
#     -nprescreen_locus <number of locus prescreen>, --nprescreen_locus=<number of locus prescreen>
#                                          number of prescreening when estimating locus tree
#     -init_locus_tree <initial locus tree>, --init_locus_tree=<initial locus tree>
#                                          initial locus tree to start search       
#   Miscellaneous
#     -x <random seed>, --seed=<random seed>
#                         random number seed
#     -l <debugging log>, -- log=<debugging log>
#                        if given, output debugging log
#     -threads <number of threads>, --threads=<number of threads>
#                       number of threads to run pll
#     -eps <epsilon>, --eps=<epsilon>
#                     criterion to determine convergence
#   File Extensions:
#     -I <input file extension>, --inputext=<input file extension>
#                         input file extension (default: "")
#     -O <output file extension>, --outputext=<output file extension>
#                         output file extension (default: ".dlca")

# by default, dlcpar outputs the reconciliation in LCT format
# this creates the files 0.dlcpar{.tree,.recon,.order}

# 5e-9 sub/site/gen * 10 gen/yr * 1e6 yr/myr = 5e-2 sub/site/myr
dlc_coestimate \
    -A sim-flies/0/0.coal.fasta.phylip -p sim-flies/0/0.coal.partitions \
    -s config/flies.stree -S config/flies.smap \
    -D 0.0012 -L 0.0012 -n 1e6 \
    -g 0.1 --subrate 5e-2 \
    -i 100 \
    --nsearch_coal 20 --nsamples_coal 1 --nprescreen_coal 20 \
    --nsearch_locus 100 --nsamples_locus 100 --nprescreen_locus 20 \
    --init_coal_tree sim-flies/0/0.raxml.tree \
    -x 1234 --threads 1 --eps 2 \
    --log \
    -I .coal.fasta.phylip -O .dlca
