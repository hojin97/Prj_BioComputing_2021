import os
import sys
import matplotlib.pyplot as plt
import datetime
import numpy as np
from collections import defaultdict
import math

BP = "BP"
MF = "MF"

def OBO_fileopen(filename):
    file = open(os.path.join(os.getcwd(), filename), 'r')
    OBO_dict = dict()
    OBO_dict[BP] = dict()
    OBO_dict[MF] = dict()
    first_line = True
    temp_id = ''
    temp_namespace = ''
    temp_relationship_is_a = []
    temp_relationship_part_of = []
    temp_obsolete = False

    for line in file.readlines():
        if line == '[Typedef]\n': break
        if line == '[Term]\n':
            if temp_obsolete or temp_namespace == 'cellular_component':
                pass
            else:
                if first_line:
                    first_line = False
                    continue
                temp_relationship = temp_relationship_is_a
                temp_relationship.extend(temp_relationship_part_of)
                temp_relationship = set(temp_relationship)

                if temp_namespace == "biological_process":
                    OBO_dict[BP][temp_id] = set(temp_relationship).copy()
                elif temp_namespace == "molecular_function":
                    OBO_dict[MF][temp_id] = set(temp_relationship).copy()

            temp_id = ''
            temp_namespace = ''
            temp_relationship_is_a = []
            temp_relationship_part_of = []
            temp_obsolete = False
        else:
            if 'alt_id: ' in line:
                pass
            elif 'id: ' in line:
                temp_id = line.replace('id: ', '')
                temp_id = temp_id.replace('\n', '')

            if 'namespace: ' in line:
                temp_namespace = line.replace('namespace: ', '')
                temp_namespace = temp_namespace.replace('\n', '')

            if 'is_a: ' in line:
                is_a = line.replace('is_a: ', '')
                is_a = is_a.split(' !')
                temp_relationship_is_a.append(is_a[0])

            if 'relationship: part_of ' in line:
                part_of = line.replace('relationship: part_of ', '')
                part_of = part_of.split(' !')
                temp_relationship_part_of.append(part_of[0])

            if 'is_obsolete: ' in line:
                temp_obsolete = True
    if temp_obsolete or temp_namespace == 'cellular_component':
        pass
    else:
        if temp_namespace == "biological_process":
            OBO_dict[BP][temp_id] = set(temp_relationship).copy()
        elif temp_namespace == "molecular_function":
            OBO_dict[MF][temp_id] = set(temp_relationship).copy()

    return OBO_dict

def FIND_root(OBO_dict):
    # term_type = BP, MF
    for term in OBO_dict[BP].keys():
        if OBO_dict[BP][term] == set():
            bp_root = term
            break

    for term in OBO_dict[MF].keys():
        if OBO_dict[MF][term] == set():
            mf_root = term
            break
    return bp_root, mf_root

def FIND_error_relationship(OBO_dict):
    for term in OBO_dict[BP].keys():
        OBO_dict[BP][term] = OBO_dict[BP][term] & OBO_dict[BP].keys()

    for term in OBO_dict[MF].keys():
        OBO_dict[MF][term] = OBO_dict[MF][term] & OBO_dict[MF].keys()

def Anotation_PPI_fileopen(filename1, filename2):
    # PPI read
    file = open(os.path.join(os.getcwd(), filename2), 'r')
    protein_dict = {}

    for protein_interaction in file.readlines():
        for p in protein_interaction.split():
            protein_dict[p] = []

    file.close()

    # Annotation read
    file = open(os.path.join(os.getcwd(), filename1), 'r')

    for annotation_line in file.readlines():
        anno_info = annotation_line.split('\t')
        # print(anno_info)

        db_symbol = anno_info[2]
        db_qualifier = anno_info[3]
        db_term = anno_info[4]
        db_evidence = anno_info[6]
        db_aspect = anno_info[8]

        if 'NOT' in db_qualifier or db_evidence == "IEA" or db_aspect == 'C':
            continue

        try:
            protein_dict[db_symbol] = protein_dict[db_symbol].append(db_term).copy()

        # Process for none exist protein.
        except Exception:
            continue
    file.close()

    key_list = list(protein_dict.keys())
    for protein in key_list:
        protein_dict[protein] = set(protein_dict[protein])
        if protein_dict[protein] == set():
            protein_dict.pop(protein)

    return protein_dict

def FIND_OBO_ancestor(OBO_dict):
    ancestor_OBO_dict = dict()
    ancestor_OBO_dict[BP] = dict()
    ancestor_OBO_dict[MF] = dict()

    # BP
    for term in OBO_dict[BP].keys():
        parents_node = set(OBO_dict[BP][term])
        remain_path = list(OBO_dict[BP][term])

        while len(remain_path) > 0:
            tmp = remain_path.pop(0)
            remain_path.extend(OBO_dict[BP][tmp] - parents_node)
            parents_node.update(OBO_dict[BP][tmp])

        ancestor_OBO_dict[BP][term] = parents_node

    # MF
    for term in OBO_dict[MF].keys():
        parents_node = set(OBO_dict[MF][term])
        remain_path = list(OBO_dict[MF][term])

        while len(remain_path) > 0:
            tmp = remain_path.pop(0)
            remain_path.extend(OBO_dict[MF][tmp] - parents_node)
            parents_node.update(OBO_dict[MF][tmp])

        ancestor_OBO_dict[MF][term] = parents_node

    return ancestor_OBO_dict

# Dict[BP/MF][term] = protein_list
def FIND_OBO_term_protein(protein_dict, ancestor_OBO_dict):
    OBO_term_protein_dict = dict()
    OBO_term_protein_dict[BP] = dict()
    OBO_term_protein_dict[MF] = dict()

    for protein in protein_dict.keys():
        bp_protein = set(protein_dict[protein]) & ancestor_OBO_dict[BP].keys()
        mf_protein = set(protein_dict[protein]) - bp_protein

        for b_term in bp_protein:
            try:
                OBO_term_protein_dict[BP][b_term].update({protein})
            except Exception:
                OBO_term_protein_dict[BP][b_term] = set({protein})

            for ancestor_term in ancestor_OBO_dict[BP][b_term]:
                try:
                    OBO_term_protein_dict[BP][ancestor_term].update({protein})
                except Exception:
                    OBO_term_protein_dict[BP][ancestor_term] = set({protein})

        for m_term in mf_protein:
            try:
                OBO_term_protein_dict[MF][m_term].update({protein})
            except Exception:
                OBO_term_protein_dict[MF][m_term] = set({protein})

            for ancestor_term in ancestor_OBO_dict[MF][m_term]:
                try:
                    OBO_term_protein_dict[MF][ancestor_term].update({protein})
                except Exception:
                    OBO_term_protein_dict[MF][ancestor_term] = set({protein})

    return OBO_term_protein_dict

def Two_term_similarity(termA, termB, term_Ancestor, root_term):
    PC1 = len(termA)/len(root_term)
    PC2 = len(termB)/len(root_term)
    PC0 = len(term_Ancestor)/len(root_term)
    return (2 * math.log2(PC0)) / (math.log2(PC1) + math.log2(PC2))

def Compare_term(T1, T2, Term_Type, ancestor_OBO_dict, OBO_term_protein_dict):
    if Term_Type == BP:
        root_term = OBO_term_protein_dict[Term_Type]['GO:0008150']
    else:
        root_term = OBO_term_protein_dict[Term_Type]['GO:0003674']

    similarity_list = np.zeros((len(T1), len(T2)))

    T1 = list(T1)
    T2 = list(T2)

    for i in range(len(T1)):
        Ti = T1[i]
        termA = ancestor_OBO_dict[Term_Type][Ti]

        # TA : Term common ancestor
        for j in range(len(T2)):
            Tj = T2[j]
            termB = ancestor_OBO_dict[Term_Type][Tj]

            common_ancestor_term = termA & termB
            common_ancestor = 1000000
            TA = None
            for tc in common_ancestor_term:
                min_v = min(common_ancestor, len(OBO_term_protein_dict[Term_Type][tc]))
                if common_ancestor > min_v:
                    TA = tc
                    common_ancestor = min_v

            term_A = OBO_term_protein_dict[Term_Type][Ti]
            term_B = OBO_term_protein_dict[Term_Type][Tj]
            if TA != None:
                term_Ancestor = OBO_term_protein_dict[Term_Type][TA]
                similarity_list[i][j] = Two_term_similarity(term_A, term_B, term_Ancestor, root_term)
            else:
                similarity_list[i][j] = 0

    sl1 = np.max(similarity_list, axis=0)
    sl2 = np.max(similarity_list, axis=1)

    return (np.sum(sl1) + np.sum(sl2)) / (len(T1) + len(T2))

def Annotation_Based_Method(filename, protein_dict, ancestor_OBO_dict, OBO_term_protein_dict):
    similarity_distribution = defaultdict(int)

    file = open(os.path.join(os.getcwd(), filename), 'r')

    for line in file.readlines():
        line = line.split()
        P1 = line[0]
        P2 = line[1]
        try:
            P1_term = protein_dict[P1]
            P2_term = protein_dict[P2]

            P1_BP = set(ancestor_OBO_dict[BP].keys()) & P1_term
            P1_MF = set(ancestor_OBO_dict[MF].keys()) & P1_term

            P2_BP = set(ancestor_OBO_dict[BP].keys()) & P2_term
            P2_MF = set(ancestor_OBO_dict[MF].keys()) & P2_term
        except Exception:
            continue

        bp_simlarity_score = 100
        mf_simlarity_score = 100

        if P1_BP == set() or P2_BP == set():
            pass
        else:
            bp_simlarity_score = Compare_term(P1_BP, P2_BP, BP, ancestor_OBO_dict, OBO_term_protein_dict)

        if P1_MF == set() or P2_MF == set():
            pass
        else:
            mf_simlarity_score = Compare_term(P1_MF, P2_MF, MF, ancestor_OBO_dict, OBO_term_protein_dict)

        if bp_simlarity_score == 100 and mf_simlarity_score == 100:
            pass
        else:
            if bp_simlarity_score == 100:
                bp_simlarity_score = 0
            if mf_simlarity_score == 100:
                mf_simlarity_score = 0

            similarity_distribution[int(max(bp_simlarity_score, mf_simlarity_score) * 10)] += 1



    return similarity_distribution


if '__main__' == __name__:
    try:
        if os.stat(sys.argv[1]).st_size == 0 or os.stat(sys.argv[2]).st_size == 0 or os.stat(sys.argv[3]).st_size == 0:
            print('No string found')

        file_name = [sys.argv[1], sys.argv[2], sys.argv[2]]

        # Ontology
        OBO_dict = OBO_fileopen(file_name[0])
        bp_root, mf_root = FIND_root(OBO_dict)
        FIND_error_relationship(OBO_dict)
        print('Complete Ontology')

        # Annotation / Dict[protein] = term
        protein_dict = Anotation_PPI_fileopen(file_name[1], file_name[2])
        print('Complete Annotation')

        # Dict[BP/MF][term] = ancestor set
        ancestor_OBO_dict = FIND_OBO_ancestor(OBO_dict)

        # Dict[BP/MF][term] = protein
        OBO_term_protein_dict = FIND_OBO_term_protein(protein_dict, ancestor_OBO_dict)
        print('Complete Protein by term')

        start = datetime.datetime.now()
        similarity_distribution = Annotation_Based_Method('PPI.txt', protein_dict, ancestor_OBO_dict, OBO_term_protein_dict)
        print((datetime.datetime.now() - start).microseconds)

        plt.bar(similarity_distribution.keys(), similarity_distribution.values())
        plt.xticks(np.arange(0, 11, 1), labels=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        plt.show()

    except Exception:
        print('INPUT FILE ERROR')