# # Blastp Alignment Parser
#
# Script for analyzing files from the "launch_blast.sh" script. The best alignment of each sequence of one species against all the sequences of another is chosen according 
# to parameters defined by the user:
#     - [per_dentity] corresponds to the minimum percentage identity desired between the two aligned sequences,
#     - [per_coverage] to the minimum coverage percentage of the alignment compared to the smallest of the two sequences,
#     - [evalue] is the minimum value of the evalue desired for the alignment.
# During execution, the user has indications on the progress of the script. The script performs the same calculations with two different methods. The first method consists
# in recording the pairs of sequences meeting the user's constraints in a dictionary then constructing a graph from this dictionary to calculate the number of cliques. 
# The second method consists in directly constructing the graph and then calculating the number of cliques.
# At the end of the execution, 3 values for each method (6 values in total) are displayed to the user:
#     - the number of reciprocal best hits, corresponding to the sequences (number of vertices pointing to each other in a directed graph)
#     - the number of cliques found in the graphical representation
#     - the execution time of the method script
#


# Libraries Loading
import glob
import os
import networkx as nx
from datetime import datetime


# Function definitions


# Dictionnary Method

# Function that records the best hits meeting the constraints in a dictionary
def find_BestHits(files, per_identity_thr, per_coverage_thr, evalue_thr):
    BestHits_Dictionary = {}
    for file in files:
        basename = file.split(".")[0]
        x = basename.split("-vs-")
        genome1 = x[0]
        genome2 = x[1]
        gene_exist = False
        with open(file, "r") as f:
            for line in f:
                if line.startswith("# Query"):
                    best_hit_found = False
                    gene1 = line.split(" ")[2].split()[0]
                    gene_1 = genome1 + "_" + gene1
                    gene_exist = True
                elif gene_exist and line.startswith(gene1):
                    if gene_1 not in BestHits_Dictionary:
                        BestHits_Dictionary[gene_1] = []
                    infos = line.split("\t")
                    if not best_hit_found:
                        gene_2 = infos[1]
                        gene_2 = genome2 + "_" + gene_2
                        id_perc = float(infos[2])
                        couverture = float(infos[3])
                        evalue = float(infos[11])
                        query_length = float(infos[13])
                        couv_perc = couverture / query_length
                        if id_perc > per_identity_thr and couv_perc > per_coverage_thr and evalue < evalue_thr:
                            BestHits_Dictionary[gene_1].append(gene_2)
                            best_hit_found = True
    return BestHits_Dictionary


# Function that saves reciprocal hits in a list
def find_reciprocal_hits(BestHits_Dictionnary):
    Reciprocal_BestHits = []
    for gene1 in BestHits_Dictionnary:
        for gene2 in BestHits_Dictionnary[gene1]:
            if gene2 in BestHits_Dictionnary and gene1 in BestHits_Dictionnary[gene2]:
                Reciprocal_BestHits.append([gene1, gene2])
    return Reciprocal_BestHits

# Function that calculates the graph associated with the dictionary of best hits
def dictionnary_to_graph(Reciprocal_BestHits):
    Graph = nx.Graph()
    Graph.add_edges_from(Reciprocal_BestHits)
    return Graph


# Graph Methods

# Function that builds a graph from the best hits
def build_graph(files, perc_id_thr, perc_coverage_thr, evalue_thr):
    pairs = set()
    Graph = nx.Graph()
    for file in files:
        basename = file.split(".")[0]
        x = basename.split("-vs-")
        genome1 = x[0]
        genome2 = x[1]
        gene_exist = False
        with open(file, "r") as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith("# Query"):
                    best_hit_found = False
                    gene1 = line.split(" ")[2].split()[0]
                    gene_1 = genome1 + "_" + gene1
                    gene_exist = True
                elif gene_exist and line.startswith(gene1):
                    infos = line.split("\t")
                    if not best_hit_found:
                        gene_2 = infos[1]
                        gene_2 = genome2 + "_" + gene_2
                        id_perc = float(infos[2])
                        couverture = float(infos[3])
                        evalue = float(infos[11])
                        query_length = float(infos[13])
                        perc_couverture = couverture / query_length
                        if id_perc > perc_id_thr and perc_couverture > perc_coverage_thr and evalue < evalue_thr:
                            pairs.add((gene_1, gene_2))
                            best_hit_found = True
                            if (gene_2, gene_1) in pairs:
                                if gene_1 not in Graph:
                                    Graph.add_node(gene_1)
                                if gene_2 not in Graph:
                                    Graph.add_node(gene_2)
                                Graph.add_edge(gene_1, gene_2)
    return Graph


# All Methods

# Function that search cliques in the graph
def find_sized_cliques(Graph, size):
    CliquesList = list(nx.find_cliques(Graph))
    h=0
    for i in range(len(CliquesList)):
        if len(CliquesList[h]) != size:
            del CliquesList[h]
        else :
                h+=1
    return CliquesList


# Main function
def main(per_identity, per_coverage, evalue):
    os.chdir("Blast_output")
    files = sorted(glob.glob('*.bl'))

    # Dictionary Method
    print("DICTIONARY METHOD")
    start_time_dic = datetime.now()
    print("   FILES PARSING")
    best_hits_dic = find_BestHits(files, per_identity_thr = per_identity,per_coverage_thr = per_coverage, evalue_thr = evalue)
    print("   RECIPROCAL HITS SEARCHING")
    reciprocal_hits_dic = find_reciprocal_hits(best_hits_dic)
    print("   BUILDING GRAPH")
    graph_dic = dictionnary_to_graph(reciprocal_hits_dic)
    print("   CLIQUES SEARCHING ")
    cliques_dic = find_sized_cliques(graph_dic, 21)
    end_time_dic = datetime.now()
    delta_time_dic = format(end_time_dic-start_time_dic)


    # Graph Method
    print("GRAPH METHOD")
    start_time_graph = datetime.now()
    print("   GRAPH BUILDING")
    graph_graph = build_graph(files, perc_id_thr=per_identity, perc_coverage_thr=per_coverage, evalue_thr=evalue)
    print("   CLIQUES SEARCHING")
    cliques_graph = find_sized_cliques(graph_graph, 21)
    end_time_graph = datetime.now()
    delta_time_graph = format(end_time_graph-start_time_graph)


    # Printing results
    print(" ")
    print(" ")
    print("=============================================")
    print("======             RESULTS             ======")
    print("=============================================")
    print("=             Dictionary method             =")
    print("=   Number of reciprocal hits : " + str(len(reciprocal_hits_dic)) + "     =")
    print("=   Number of cliques : "+ str(len(cliques_dic)) +"                =")
    print("=   Duration : "+ delta_time_dic +"               =")
    print("=============================================")
    print("=                Graph method               =")
    print("=   Number of edges : " + str(graph_graph.number_of_edges()) + "                =")
    print("=   Number of cliques : "+ str(len(cliques_graph)) +"                =")
    print("=   Duration : "+ delta_time_graph +"               =")
    print("=============================================")
    print(" ")
    print(" ")
    
# Running the script
main(per_identity=70, per_coverage=0, evalue=float("1e-10"))
