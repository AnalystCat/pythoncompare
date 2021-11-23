import glob
import time

# Début du programme, on calcule le temps d'exécution
start_time = time.time()

# Stockage des fichiers .bl
files = []

for fi in glob.glob("*.bl"):
    files.append(fi)


def get_best_hit(fichier):
    """
    Donne le best hit pour chaque gène du génome 1
    :param fichier: genome1-vs-genome2.bl
    :return: blastp dico où chaque clé correspond à un gène de génome 1 et sa valeur est son best hit dans génome 2
    """

    f = open(fichier)
    lines = f.readlines()

    # Initialisation des variables
    blastp = {}
    hit = 0 # vaut 1 si la ligne est un hit, 0 sinon
    signif = 0
    # signif vaut 1 si le hit est significatif (identité > 70%, e-value < 0.0000000001)
    # 2 s'il est meilleur que le 1er hit

    # On parcourt le fichier
    for line in lines:

        if line.startswith("#"):
            hit = 0 # la ligne n'est pas un hit
            signif = 0

        else:
            if hit == 0: # la ligne d'avant n'était pas un hit
                gene = line.split("\t")[0] # on stocke le nom du gène
            hit = 1 # la ligne est un hit

            # Si le hit est significatif
            if float(line.split("\t")[2]) > 70 and float(line.split("\t")[11]) < 0.0000000001:

                if signif == 1:
                    # Le hit est meilleur que celui d'avant
                    if identite >= float(line.split("\t")[2]) and evalue <= float(line.split("\t")[11]):
                        blastp[gene] = line.split("\t")[1]
                        signif = 2

                else:
                    signif = 1
                    blastp[gene] = line.split("\t")[1]
                    identite = float(line.split("\t")[2])
                    evalue = float(line.split("\t")[11])                    

    f.close()

    return blastp


def get_gene_dict(fichier):
    """
    Renvoie un dictionnaire où les clés sont les noms des gènes du génome et les valeurs associées sont "None"
    :param fichier: fichier .fa
    :return: genes_dict
    """

    f = open(fichier)
    lines = f.readlines()

    genes = []

    for line in lines:
        if line.startswith(">"): # la ligne est un identifiant de séquence
            genes.append(line[1:-1])

    f.close()

    genes_dict = dict.fromkeys(genes)

    return genes_dict


# Stockage des fichiers .fa
fa_files = []

for fi in glob.glob("*.fa"):
    fa_files.append(fi)

# Dictionnaire où chaque clé est le nom d'un génome
# Pour chaque génome, un autre dictionnaire où les clés sont les gènes du génome
# Pour chaque gène, la valeur associée est la liste des best hits
dico_genome = {}

# Dictionnaire qui associe à chaque génome le préfixe commun à chacun de ses gènes (ex : Eco1 pour Escherichia_coli_536)
prefixe_genes = {}

# Remplissage de dico_genome avec les noms de génomes et les gènes associés
for fi in fa_files:
    # Stockage du préfixe
    f = open(fi)
    firstline = f.readline()
    f.close()
    prefixe_genes[firstline.split("_")[0][1:]] = fi[:-3]

    genome_name = fi[:-3]
    dico_genome[genome_name] = {}
    dico_genome[genome_name] = get_gene_dict(fi)
    for gene in dico_genome[genome_name]:
        dico_genome[genome_name][gene] = []

# Remplissage des best hits pour chaque gène de chaque génome
for fi in files:
    best_hits = get_best_hit(fi)
    for gene in dico_genome[fi.split("-vs-")[0]]:
        for cle in best_hits:
            if cle == gene:
                dico_genome[fi.split("-vs-")[0]][gene].append(best_hits[cle])

# Liste des gènes du core génome
core_genome = []

# Remplissage du core génome
for genome_name in dico_genome:
    for gene in dico_genome[genome_name]:
        i = 0
        if len(dico_genome[genome_name][gene]) == len(fa_files): # le gène a un best hit pour chaque génome
            for best_hit in dico_genome[genome_name][gene]:
                if best_hit.split("_")[0] != gene.split("_")[0]:
                    genome = prefixe_genes[best_hit.split("_")[0]]
                    if gene in dico_genome[genome][best_hit]:
                        # Le best hit a lui-même pour best hit "gene" (le best hit est réciproque)
                        i = i + 1
        if i == len(fa_files)-1: # tous les best hits de gene ont "gene" pour best hit (réciprocité)
            core_genome.append(gene)

# On stocke les gènes du core génome dans un fichier
with open('core_genome.txt', 'w') as file:
    for gene in core_genome:
        file.write(gene + "\n")

print("Temps d'exécution : %s secondes" % (time.time() - start_time))
