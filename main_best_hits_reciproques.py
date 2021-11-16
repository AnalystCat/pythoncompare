import glob

files = []

for fi in glob.glob("*.bl"):
    files.append(fi)

"""def get_best_hit(fichier):
    f = open(fichier)
    lines = f.readlines()
    blastp = {}

    hit = 0
    signif = 0

    for line in lines:
        if line.startswith("#"):
            hit = 0 # la ligne n'est pas un hit
            signif = 0
        else:
            if hit == 0: # la ligne d'avant n'était pas un hit
                gene = line.split("\t")[0]
            hit = 1 # la ligne est un hit

            if float(line.split("\t")[2]) > 70 and float(line.split("\t")[11]) < 0.0000000001:
                if signif == 1:
                    if identite < float(line.split("\t")[2]) or evalue > float(line.split("\t")[11]):
                        blastp[gene] = line.split("\t")[1]
                        signif = 2
                else:
                    signif = 1
                    blastp[gene] = line.split("\t")[1]
                    identite = float(line.split("\t")[2])
                    evalue = float(line.split("\t")[11])                    

    f.close()

    return blastp

def hits_reciproques(dico1, dico2):
    blastp = {}

    for cle, valeur in dico1.items():
        if valeur in dico2:
            if cle == dico2.get(valeur):
                blastp[cle] = valeur

    return blastp

best_hit = {}

for fi in files:
    best_hit[fi] = get_best_hit(fi)

hits = {}
done = []

for fi in files:
    if fi not in done:
        dico1 = best_hit[fi]
        reciproque = fi.split("-vs-")[1][:-3] + "-vs-" + fi.split("-vs-")[0] + ".bl"
        dico2 = best_hit[reciproque]
        hits[fi] = hits_reciproques(dico1, dico2)
        done.append(reciproque)"""

def get_gene_list(fichier):
    f = open(fichier)
    lines = f.readlines()
    genes = []

    for line in lines:
        if line.startswith(">"):
            genes.append(line[1:-1])

    f.close()

    return genes

fa_files = []

for fi in glob.glob("*.fa"):
    fa_files.append(fi)

dico_genome = {}

for fi in fa_files:
    genome_name = fi[:-3]
    dico_genome[fi] = {}
    dico_genome[fi] = get_gene_list(fi)

"""with open('hits_reciproques.txt', 'w') as file:
    for paire in hits:
        file.write(paire + "\n")
        for cle, valeur in hits[paire].items():
            file.write(cle + " " + valeur + "\n")
        file.write("\n")"""
    