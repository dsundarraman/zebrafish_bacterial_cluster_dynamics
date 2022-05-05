
import numpy as np
import random
from glob import glob
from pathlib import Path
from functions_gillespie import expulsion, aggregation, growth, fragmentation, expelled_singles

"K is the carrying capacity, the mean growth rate of individuals and clusters"
" (from experiments) 0.24 /hr and 0.66 /hr. Agg_rate is the value of alpha set to be 10**(-2.5) /hr. This simulation assumes" \
" size dependence of the rates for both aggregation and fragmentation with exponent 1/3 and 2/3. " \
"The fragmentation rates are generated in frag_rate_list at steps of df = 0.025 between 0.025 and 1 /hr. The expulsion" \
"rate of aggregates is fixed and set to 0.1 /hr and the expulsion rate of individuals is varied logarithmically" \
"in r_exp_list. One of the saved initial configurations (generated from the mono-association simulation) is  chosen. " \
"r_exp_list and frag_rate_list are then chosen and the simulation is run. The total time for the simulation" \
"T = 24hrs after which the final configuration is saved. All rates are in /hr"

K = 10000

growth_rate_inds = 0.24
growth_rate_clumps = 0.66

agg_rate = 10**(-2.5)
power_aggregation = 1/3

power_fragmentation = 2/3
frag_rate_list = [np.round(i * 0.025, 2) for i in range(1,41)]

r_exp_aggs = 0.1
exp_factor= np.logspace(np.log10(0.25), 1, num=10)
r_exp_list = [np.round(i * r_exp_aggs, 3) for i in exp_factor]

T_total = 24

direc_files = glob('/media/rplab/Deepika backup1/cluster_model_march/4_13_22/initial/*.npz')

time = 0
fish = 0
total_fish = 500


#### each function, expulsion, aggregation and fragmentation, calculates the rates for each process
#### and which clusters are undergoing the reaction


"""Simulation"""

for exp in range(len(r_exp_list)):
    ## choose the expulsion rate for the simulation and generate a path name with this expulsion rate for saving later ##

    r_exp_inds = r_exp_list[exp]
    save_direc = '/media/rplab/Deepika backup1/cluster_model_march/4_19_22/exp_' + str(r_exp_inds) + '/'
    Path(save_direc).mkdir(parents=True, exist_ok=True)

    for f in range(len(frag_rate_list)):

        ## choose the fragmentation rate for the simulation ##########

        frag_rate = frag_rate_list[f]
        print('starting f =' + str(frag_rate) + 'expulsion rate individuals = ' + str(r_exp_inds))

        for fish in range(total_fish):

            initial_cluster_file = np.random.choice(direc_files, 1)
            clusters = list(np.load(initial_cluster_file[0])['arr_0'])

            time = 0

            while time < T_total:

                ### calculate rates of each of the reactions in this time step and then calculate the total reaction ##
                ### rate.

                clumps, singles, cluster_exp_rates = expulsion(clusters, r_exp_aggs)
                if cluster_exp_rates == 'Extinct':

                    break

                cluster_agg_rates, pairs_list, non_aggregating_loner = aggregation(clusters, agg_rate, power_aggregation)
                cluster_frag_rates, clusters_greater_than_one, new_clusters = fragmentation(clusters, frag_rate, power_fragmentation)

                all_rates = [cluster_agg_rates, cluster_exp_rates, cluster_frag_rates]
                flattened_rates = [item for sublist in all_rates for item in sublist]

                total_reaction_rate = sum(flattened_rates)
                if total_reaction_rate == 0:
                    clusters = [0]
                    break

                ### calculate the time to the next reaction ###
                tau = np.random.exponential(1 / total_reaction_rate)
                u1 = np.random.random()
                time = tau + time

                probabilities = [i/ total_reaction_rate for i in flattened_rates]
                ### pick reaction based on probabilities ###
                i_reaction = np.random.choice(len(probabilities), p=probabilities)

                ### figure out which reaction is occuring based on the value of i_reaction and perform the reaction  ###

                if len(cluster_agg_rates) > i_reaction:
                    ## aggregation reaction, just take the sum of the specific pair aggregating and leave remainign as is ###
                    new_clusters = [np.array(np.sum(pairs_list[i])) if i == i_reaction else pairs_list[i] for \
                                   i in range(len(pairs_list))]
                    clusters = np.concatenate([x.ravel() for x in new_clusters])
                    if len(non_aggregating_loner) > 0:
                        clusters = np.concatenate([clusters, non_aggregating_loner])

                elif len(cluster_agg_rates) + len(cluster_exp_rates) > i_reaction:
                    ## expulsion reaction, remove cluster from list ###
                    i_exp = i_reaction - len(cluster_agg_rates)
                    clusters = list(np.delete(clumps, i_exp))
                    clusters.extend(singles)

                elif len(cluster_agg_rates) + len(cluster_exp_rates) + len(cluster_frag_rates) > i_reaction \
                        and len(cluster_agg_rates) + len(cluster_exp_rates) < i_reaction + 1:
                    ## fragmentation reaction, find out which cluster fragments and perform reaction based on rate ###
                    i_frag = i_reaction -  len(cluster_agg_rates) - len(cluster_exp_rates)
                    fragmented_clump = clusters_greater_than_one[i_frag]  # cluster undergoing fragmentation #
                    clusters = [i for i in clusters_greater_than_one if i != fragmented_clump]
                    cells_removed = fragmented_clump * flattened_rates[i_reaction] * tau
                    new_cluster_size = int(fragmented_clump - cells_removed)
                    cells_fragmented = int(cells_removed) * [1]
                    cells_fragmented.append(new_cluster_size)

                    if len(cells_fragmented) > 0:
                        clusters.extend(cells_fragmented)
                    else:
                        clusters.extend([fragmented_clump])

                    if len(new_clusters) > 0:
                        clusters.extend(new_clusters)
                n_tot = np.sum(clusters)

                ### growth of clusters ###
                clusters = growth(clusters, n_tot, tau, K, growth_rate_inds, growth_rate_clumps)
                ### expulsion individuals based on rate and time step ###
                clusters = expelled_singles(clusters, r_exp_inds, tau)

            exp_save_direc = save_direc + '/frag_' + str(frag_rate)
            Path(exp_save_direc).mkdir(parents=True, exist_ok=True)
            print(np.sum(clusters))

            #np.savez_compressed(exp_save_direc + '/fish_' + str(fish), clusters, np.sum(clusters))
            print('Completed fish ' + str(fish + 1) + ' of ' + str(total_fish) +  ' fish')