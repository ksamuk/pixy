#!/usr/bin/env python
# coding: utf-8

#! conda install -y -c conda-forge msprime
import msprime
import argparse

# define arguments for the parser
parser = argparse.ArgumentParser(description='Simulate polymorphism data with msprime and write to VCF file(s)')
parser.add_argument('--prefix', type=str, nargs='?', help='file name prefix for output VCF(s)')
parser.add_argument('--datasets', type=int, nargs='?', help='number of datasets to simulate')

# pull out varizbles from the parsed args
# sim_args = parser.parse_args(['--prefix', 'data/msprime_sim/pi_sim_', '--datasets', '10'])
sim_args = parser.parse_args()

prefix = sim_args.prefix
datasets = sim_args.datasets

# define parameters + simulate some data
Ne = 1e6
Mu = 1e-8
sample_size = 100
n_sites = 10000
iterations = list(range(0, datasets, 1))

for iteration in iterations:
    # simulate the specified number of iterations/sites/samples
    tree_seq = msprime.simulate(sample_size=sample_size, Ne=Ne, 
                                length=n_sites, mutation_rate=Mu, random_seed = iteration + 1)
    # write to file
    filename = prefix + "Ne=" + str('%.1e' % Ne) + "_" + "mu=" + str(Mu) + "_" + "samples=" + str(sample_size) + "_" + "sites=" + str(n_sites) + "_" + str(iteration) + '.vcf'
    with open(filename, "w") as vcf_file:
        tree_seq.write_vcf(vcf_file, ploidy=2)



