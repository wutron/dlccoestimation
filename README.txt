DLC_COEST (duplications, losses, and coalescence coestimation)
http://XXX
Yi-Chieh Wu, Bo Zhang

=============================================================================
ABOUT

DLC_COEST is a software package containing the dlc_coest reconciliation program 
for working with gene trees.

dlc_coest is a reconciliation program that takes in a DNA sequence data and a species 
tree and determines how gene duplications and losses have occurred, under the presence 
of incomplete linegae sorting.

DLC_COEST citation: 


=============================================================================
DEPENDENCIES

DLC_COEST has the following requirements:

- GNU Scientific library (GSL) http://www.gnu.org/software/gsl/
- Python (2.5 or greater) http://python.org/


=============================================================================
INSTALL

NOTE: Makefile installation will work best on UNIX or CYGWIN (Windows).


To compile the DLC_COEST library use the Makefile.

    make

Once compiled, to install the DLC_COEST programs (default install in /usr) use:

    make install

To specify your own installation path use:
    
    make install prefix=/usr/local

DLC_COEST can also run directly from the source directory.  Simply add the
bin/ directory to your PATH or create symlinks to the scripts within bin/
to any directory on your PATH.


=============================================================================
USAGE

Running dlcoal_recon with no arguments will print out its command-line usage:


Usage: dlc_coestimate [options]

Options:
    -h, --help            show this help message and exit
    -A ALIGNMENT,    --alignment=ALIGNMENT
                        sequence data (phylip format))
    -p PARTITION_FILE  --partitions=PARTITION_FILE
                       partition file of DNA sequence data 
    -s SPECIES_TREE, --stree=SPECIES_TREE
                        species tree file in newick format (myr)
    -S GENE_TO_SPECIES_MAP, --smap=GENE_TO_SPECIES_MAP
                        gene to species map
  
  Model Parameters:
    -n POPULATION_SIZE, --popsize=POPULATION_SIZE
                        Effective population size
    -D DUPLICATION_RATE, --duprate=DUPLICATION_RATE
                        rate of a gene duplication (dups/gene/myr)
    -L LOSS_RATE, --lossrate=LOSS_RATE
                        rate of gene loss (losses/gene/myr)
    -g GENRATION_TIME, --gentime=GENRATION_TIME
                        generation time (years)
    --subrate SUBSTITUTION_RATE, --subrate=SUBSTITUTION_RATE
                        substitution rate in sub/site/myr
    --beta BETA,    --beta=BETA
                   log prob = beta * log likelihood + log prior
 
  Search Parameters:
    -i ITERATIONS, --iter=ITERATIONS
                        number of search iterations
    --nsearch_coal NUMBER_OF_COAL_SEARCH, --nsearch_coal=NUMBER_OF_COAL_SEARCH
          number of search iterations for the coalescent tree when fixing the locus tree
    --nsamples_coal NUMBER_OF_COAL_SAMPLES, --nsamples_coal=NUMBER_OF_COAL_SAMPLES
         number of samples taken when doing integration over coalescent tree branch lengths
    --nprescreen_coal NUMBER_OF_COAL_PRESCREENS, --nprescreen_coal=NUMBER_OF_COAL_PRESCREENS
         number of prescreening iterations when estimating the coalescent tree
    --init_coal_tree INITIAL_COAL_TREE, --init_coal_tree=INITIAL_COAL_TREE
                                         initial coal tree for search
    --nsearch_locus NUMBER_OF_LOCUS_SEARCH, --nsearch_locus=NUMBER_OF_LOCUS_SEARCH
          number of search iterations for the locus tree when fixing the coalescnet tree
    --nsamples_locus NUMBER_OF_LOCUS_SAMPLES, --nsamples_locus=NUMBER_OF_LOCUS_SAMPLES
         number of samples taken when doing integration over locus tree branch lengths
    --nprescreen_locus NUMBER_OF_LOCUS_PRESCREENS, --nprescreen_locus=NUMBER_OF_LOCUS_PRESCREENS
         number of prescreening iterations when estimating the locus tree
    --init_locus_tree INITIAL_LOCUS_TREE, --init_locus_tree=INITIAL_LOCUS_TREE
                                         initial locus tree for search

  File extensions:
    -I INPUT_EXT, --inext=INPUT_EXT
                        input file extension (default='')
    -O OUTPUT_EXT, --outext=OUTPUT_EXT
                        output file extension (default='.dlca')

  Miscellaneous:
    -x RANDOM_SEED, --seed=RANDOM_SEED
                        random number seed
    -l, --log           if given, output debugging log
    --threads NUMBER_OF_THREADS, --threads=NUMBER_OF_THREADS
                        number of threads to run pll
    --eps EPSILON, --eps=EPSILON
                        criterion to determine convergence
#=============================================================================
# Examples

See examples/make.sh for an example of how to use each program
in the DLCoal package.


#=============================================================================
# Documentation