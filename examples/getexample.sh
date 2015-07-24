path=~/research/work/coestimation/simulation

# get example from project (sim-flies, 1e6, 1x, famid 0)
mkdir config
cp -L $path/config/flies{.stree,.smap} config/

mkdir -p sim-flies/0
cp -L $path/data/1000/5e-9/1e6-1x/0/0{.coal.partitions,.coal.fasta.phylip} sim-flies/0/
cp -L $path/data/1000/5e-9/1e6-1x/0/0.raxml.dlcoal.coal.tree sim-flies/0/0.raxml.tree
