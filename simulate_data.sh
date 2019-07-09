#!/bin/bash

# simulate invariant-sites vcfs using msprime and a helper script
# pop gen parameters are hard-coded in script, but included in file name

rm -r data/msprime_sim
rm -r data/msprime_sim_invar

mkdir -p data/msprime_sim
mkdir -p data/msprime_sim_invar

# simulate data using the msprime wrapper script
python scripts/msprime_simulate_data.py --prefix data/msprime_sim/sim_dat_ --datasets 10

# inject invariant sites
sh scripts/inject_invariant_sites.sh data/msprime_sim data/msprime_sim_invar