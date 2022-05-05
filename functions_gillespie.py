
import numpy as np
import random
from glob import glob
from pathlib import Path


def expulsion(all_clusters, r_exp_aggs):
    ## generates a list of aggregates (cluster size > or = 2) with cluster expulsion rates.
    if len(all_clusters) > 0:
        clusters_greater_than_one = [i for i in all_clusters if i > 2 or i == 2]
        individuals = [i for i in all_clusters if i < 2]
        cluster_rates = [r_exp_aggs] * len(clusters_greater_than_one)
    else:
        cluster_rates = 'Extinct'
        clusters_greater_than_one = []
        individuals = []
    return clusters_greater_than_one, individuals, cluster_rates

def aggregation(all_clusters, agg_rate, power_aggregation):
    ## if even number of clusters present, choose a random cluster that doesn't participate in this reaction- non_aggregating_loner##
    ## the remaining clusters are paired up and randomly and the rates of aggregation are calculated. Cluster rates is a list of rates  ###
    ## calculated for each pair of aggregating clusters. Pairs list keeps track of which clusters are paired up in order of cluster rates.##

    if len(all_clusters) % 2 != 0:
        not_aggregating = random.choice(list(enumerate(all_clusters)))[0]
        aggregating_clusters = [all_clusters[i] for i in range(len(all_clusters)) if i!= not_aggregating]
        non_aggregating_loner = [all_clusters[not_aggregating]]
    else:
        aggregating_clusters = all_clusters
        non_aggregating_loner = []

    if len(aggregating_clusters) > 0:
        pairs_list = np.random.choice(aggregating_clusters, size=(int(np.floor(len(aggregating_clusters)/2)), 1, 2), replace=False)
        cluster_rates = [agg_rate*(((pairs_list[k][0][0])*(pairs_list[k][0][1]))**(power_aggregation)) for k in range(len(pairs_list))]
    else:
        cluster_rates = [0]
        pairs_list = non_aggregating_loner

    return cluster_rates, pairs_list, non_aggregating_loner


def fragmentation(all_clusters, frag_rate, power_fragmentation):
    ### only clusters can fragment, so generate a separate list of fragmenting clusters and calculate their rates ###
    clusters_greater_than_one = [i for i in all_clusters if i > 2 or i == 2]
    individuals = [i for i in all_clusters if i < 2]
    cluster_rates = [(k**(power_fragmentation))*frag_rate for k in clusters_greater_than_one]

    return cluster_rates, clusters_greater_than_one, individuals

def growth(all_clusters, n_tot, dt, K, growth_rate_inds, growth_rate_clumps):
    ### treat growth of individuals and aggregates separately, based on their respective rates ####
    ### note that the concatenated growth_aggs consists of both the individuals and aggregates ####

    single_cells = [i for i in all_clusters if i < 2]
    clusters_greater_than_one = [i for i in all_clusters if i > 2 or i == 2]

    growth_inds = np.sum(single_cells) + growth_rate_inds * np.sum(single_cells) * (1 - (n_tot / K)) * dt
    growth_aggs = [np.round(cluster, 1) + growth_rate_clumps * np.round(cluster,1) * (1 - (n_tot / K)) * dt
                     for cluster in clusters_greater_than_one]
    growth_aggs.extend([1] * int(growth_inds))
    return growth_aggs

def expelled_singles(all_clusters, r_exp_inds, tau):

    ### expulsiom of indovodials simply involves that at ever time step tau, calculate the number of single cells expelled ####
    ### which is determined by the chosen value of r_exp_inds. ####

    single_cells = [i for i in all_clusters if i < 2]
    clusters_greater_than_one = [i for i in all_clusters if i > 2 or i == 2]

    number_inds_expelled = int(np.ceil(r_exp_inds * len(single_cells) * tau))
    single_cells = single_cells[number_inds_expelled::]
    clusters_greater_than_one.extend(single_cells)

    return clusters_greater_than_one

