import os
import sys
import matplotlib.pyplot as plt
import datetime
import numpy as np
from collections import defaultdict
import math
from sklearn.manifold import TSNE
import pandas as pd
import copy

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

def FIND_depth(OBO_dict, OBO_depth_dict, root, type):
    for term in OBO_dict[type].keys():
        remain = list(OBO_dict[type][term])
        have_been_path = set()
        move_cnt = 1
        if {term} == {root}:
            OBO_depth_dict[type][term] = 0
        else:
            while True:
                if {root} & set(remain) == {root}:
                    break
                tmp = []
                for r_term in remain:
                    tmp.extend(OBO_dict[type][r_term] - have_been_path)
                    have_been_path.update(OBO_dict[type][r_term])
                remain = tmp.copy()
                move_cnt += 1

            OBO_depth_dict[type][term] = move_cnt
    return OBO_depth_dict

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

def Annotation_similarity(termA, termB, term_Ancestor, root_term):
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
                similarity_list[i][j] = Annotation_similarity(term_A, term_B, term_Ancestor, root_term)
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

def Edge_Based_Method(filename, protein_dict, ancestor_OBO_dict, OBO_depth_dict, OBO_dict):
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

        # Compare BP
        if P1_BP == set() or P2_BP == set():
            pass
        elif P1_BP == {'GO:0008150'} or P2_BP == {'GO:0008150'}:
            bp_simlarity_score = 0
        else:
            bp_simlarity_score = Edge_Compare_term(P1_BP, P2_BP, BP, ancestor_OBO_dict, OBO_depth_dict, OBO_dict)

        # Compare MF
        if P1_MF == set() or P2_MF == set():
            pass
        elif P1_MF == {'GO:0003674'} or P2_MF == {'GO:0003674'}:
            mf_simlarity_score = 0
        else:
            mf_simlarity_score = Edge_Compare_term(P1_MF, P2_MF, MF, ancestor_OBO_dict, OBO_depth_dict, OBO_dict)

        # Accumulate similarity score
        if bp_simlarity_score == 100 and mf_simlarity_score == 100:
            pass
        else:
            if bp_simlarity_score == 100:
                bp_simlarity_score = 0
            if mf_simlarity_score == 100:
                mf_simlarity_score = 0

            similarity_distribution[int(max(bp_simlarity_score, mf_simlarity_score) * 10)] += 1

    file.close()
    return similarity_distribution

# T1, T2 = set( term )
def Edge_Compare_term(T1, T2, TYPE, ancestor_OBO_dict, OBO_depth_dict, OBO_dict):
    if TYPE == BP:
        ROOT = 'GO:0008150'
    elif TYPE == MF:
        ROOT = 'GO:0003674'

    similarity_list = np.zeros((len(T1), len(T2)))

    for i, C1 in enumerate(T1):
        for j, C2 in enumerate(T2):
            if C1 == ROOT or C2 == ROOT:
                similarity_list[i][j] = 0
            else:
                T0 = ancestor_OBO_dict[TYPE][C1] & ancestor_OBO_dict[TYPE][C2]
                C0_Root_dist = 0
                for c0 in T0:
                    if C0_Root_dist <= OBO_depth_dict[TYPE][c0]:
                        C0_Root_dist = OBO_depth_dict[TYPE][c0]
                        C0 = c0

                if C0 == ROOT:
                    similarity_list[i][j] = 0
                else:
                    C0_C1_dist = Between_C0_term(C0, C1, TYPE, OBO_dict)
                    C0_C2_dist = Between_C0_term(C0, C2, TYPE, OBO_dict)

                    # Edge Similarity
                    sim = (2 * C0_Root_dist) / (C0_C1_dist + C0_C2_dist + 2 * C0_Root_dist)
                    similarity_list[i][j] = sim

    sl1 = np.max(similarity_list, axis=0)
    sl2 = np.max(similarity_list, axis=1)

    return (np.sum(sl1) + np.sum(sl2)) / (len(T1) + len(T2))

def Between_C0_term(C0, C1, TYPE, OBO_dict):
    remain = list(OBO_dict[TYPE][C1])
    have_been_path = set()
    move_cnt = 1
    while True:
        if {C0} & set(remain) == {C0}:
            break
        tmp = []
        for r_term in remain:
            tmp.extend(OBO_dict[TYPE][r_term] - have_been_path)
            have_been_path.update(OBO_dict[TYPE][r_term])
        remain = tmp.copy()
        move_cnt += 1

    return move_cnt

def PPI_Similar_CSVopen(filename):
    PPI_Similarity_dict = dict()

    file = open(os.path.join(os.getcwd(), filename), 'r')
    first_line = False
    for line in file.readlines():
        if not first_line:
            first_line = True
            continue
        line = line.split('\t')
        ppi_set = line[0]
        ppi_set = ppi_set.replace('"', '')
        ppi_set = ppi_set.replace('\'', '')
        ppi_set = ppi_set.replace('\\', '')
        ppi_set = ppi_set.replace('(', '')
        ppi_set = ppi_set.replace(')', '')
        ppi_set = ppi_set.replace(',', '')

        node_score = float(line[1])
        anno_score = float(line[2])
        edge_score = float(line[3].replace('\n', ''))

        PPI_Similarity_dict[ppi_set] = [node_score, anno_score, edge_score]

    return PPI_Similarity_dict

def Term_vectorization(OBO_depth_dict, ancestor_OBO_dict, OBO_term_protein_dict):
    OBO_term_vector = OBO_depth_dict.copy()

    for term in list(OBO_term_vector[BP].keys()):
        try:
            depth = OBO_term_vector[BP][term]
            num_parents = ancestor_OBO_dict[BP][term]
            num_annotations = len(OBO_term_protein_dict[BP][term])
            infomation_content = len(OBO_term_protein_dict[BP][term]) / len(OBO_term_protein_dict[BP]['GO:0008150'])
            OBO_term_vector[BP][term] = [depth, num_parents, num_annotations, -1 * math.log2(infomation_content)]
        except Exception:
            OBO_term_vector[BP].pop(term)

    for term in list(OBO_term_vector[MF].keys()):
        try:
            depth = OBO_term_vector[MF][term]
            num_parents = ancestor_OBO_dict[MF][term]
            num_annotations = len(OBO_term_protein_dict[MF][term])
            infomation_content = len(OBO_term_protein_dict[MF][term]) / len(OBO_term_protein_dict[MF]['GO:0003674'])
            OBO_term_vector[MF][term] = [depth, num_parents, num_annotations, -1 * math.log2(infomation_content)]
        except Exception:
            OBO_term_vector[MF].pop(term)

    return OBO_term_vector

def Embbed_term(TERM, TYPE, OBO_term_vector, TERM_deepest_len):
    if len(TERM) == 1:
        vector = OBO_term_vector[TYPE][str(TERM.pop())].copy()
        vector[1] = len(vector[1])
        return vector
    else:
        tmp_vector = []
        for term in TERM:
            tmp_vector.append(OBO_term_vector[TYPE][term])

        tmp_vector = np.array(tmp_vector)

        DEPTH = np.max(tmp_vector[:, 0])

        ancestor_set = list(tmp_vector[:, 1])

        ancestor_union = ancestor_set.pop(0)
        ancestor_intersection = ancestor_union.copy()

        for a in ancestor_set:
            ancestor_union = ancestor_union | a
            ancestor_intersection = ancestor_intersection & a

        ANCESTOR = len(ancestor_union) - len(ancestor_intersection)

        ANNOTATION = np.min(tmp_vector[:, 2])

        INFORM_CONTENT = 0
        for depth, inform_content in zip(tmp_vector[:, 0], tmp_vector[:, 3]):
            INFORM_CONTENT = INFORM_CONTENT + inform_content * depth / TERM_deepest_len

        return [DEPTH, ANCESTOR, ANNOTATION, INFORM_CONTENT]

def Protein_vectorization(filename, protein_dict, OBO_term_vector, BP_largest_len, MF_largest_len):
    Protein_vector_dict = {
        BP:dict(),
        MF:dict()
    }
    file = open(os.path.join(os.getcwd(), filename), 'r')

    protein_list_dict = {}

    for protein_interaction in file.readlines():
        for p in protein_interaction.split():
            protein_list_dict[p] = []

    file.close()

    for protein in protein_list_dict.keys():
        try:
            P1_term = protein_dict[protein]
        except Exception:
            continue

        P1_BP = set(ancestor_OBO_dict[BP].keys()) & P1_term
        P1_MF = set(ancestor_OBO_dict[MF].keys()) & P1_term

        if P1_BP == set():
            BP_vector = [0,0,0,0]
        else:
            BP_vector = Embbed_term(P1_BP, BP, OBO_term_vector, BP_largest_len)

        if P1_MF == set():
            MF_vector  = [0,0,0,0]
        else:
            MF_vector = Embbed_term(P1_MF, MF, OBO_term_vector, MF_largest_len)

        Protein_vector_dict[BP][protein] = BP_vector
        Protein_vector_dict[MF][protein] = MF_vector

    return Protein_vector_dict

def Protein_Visualization(Protein_vector_dict):
    pvd = copy.deepcopy(Protein_vector_dict)

    for term in list(pvd[BP].keys()):
        if pvd[BP][term] == [0, 0, 0, 0]:
            pvd[BP].pop(term)

    for term in list(pvd[MF].keys()):
        if pvd[MF][term] == [0, 0, 0, 0]:
            pvd[MF].pop(term)

    bp_target = list(pvd[BP].keys()).copy()
    mf_target = list(pvd[MF].keys()).copy()

    # target = list(Protein_vector_dict[BP].keys())
    tsne = TSNE(random_state=30, perplexity=10)

    # plt.figure(0)
    # plt.title('BP')
    tsne_result_bp = tsne.fit_transform(list(pvd[BP].values()))
    tsne_result_bp = pd.DataFrame(tsne_result_bp, columns=['tsne1', 'tsne2'])
    # plt.scatter(tsne_result_bp['tsne1'], tsne_result_bp['tsne2'])
    # for x, y, t in zip(tsne_result_bp['tsne1'], tsne_result_bp['tsne2'], target):
    #     plt.text(x, y, str(t))

    # plt.figure(1)
    # plt.title('MF')
    tsne_result_mf = tsne.fit_transform(list(pvd[MF].values()))
    tsne_result_mf = pd.DataFrame(tsne_result_mf, columns=['tsne1', 'tsne2'])
    # plt.scatter(tsne_result_mf['tsne1'], tsne_result_mf['tsne2'])
    # for x, y, t in zip(tsne_result_mf['tsne1'], tsne_result_mf['tsne2'], target):
    #     plt.text(x, y, str(t))

    # plt.show()

    return tsne_result_bp, tsne_result_mf, bp_target, mf_target

def Distance_based_method(filename, tsne_result_bp, tsne_result_mf, bp_target, mf_target):
    similarity_distribution = defaultdict(int)
    proetin_similarity_dict = dict()

    tsne_result_bp = np.array(tsne_result_bp)
    tsne_result_mf = np.array(tsne_result_mf)

    protein_bp_min = np.min(tsne_result_bp, axis=0)
    protein_bp_max = np.max(tsne_result_bp, axis=0)

    protein_mf_min = np.min(tsne_result_mf, axis=0)
    protein_mf_max = np.max(tsne_result_mf, axis=0)

    file = open(os.path.join(os.getcwd(), filename), 'r')

    for line in file.readlines():
        line = line.split()
        P1 = line[0]
        P2 = line[1]

        try:
            P1_bp = tsne_result_bp[bp_target.index(P1)]
        except Exception:
            P1_bp = None

        try:
            P2_bp = tsne_result_bp[bp_target.index(P2)]
        except Exception:
            P2_bp = None

        try:
            P1_mf = tsne_result_mf[mf_target.index(P1)]
        except Exception:
            P1_mf = None

        try:
            P2_mf = tsne_result_mf[mf_target.index(P2)]
        except Exception:
            P2_mf = None

        if type(P1_bp) == type(None) or type(P2_bp) == type(None):
            bp_score = -10
        else:
            bp_x_dist = abs((P1_bp[0] - protein_bp_min[0])/(protein_bp_max[0] - protein_bp_min[0]) \
                        - (P2_bp[0] - protein_bp_min[0])/(protein_bp_max[0] - protein_bp_min[0]))

            bp_y_dist = abs((P1_bp[1] - protein_bp_min[1])/(protein_bp_max[1] - protein_bp_min[1]) \
                        - (P2_bp[1] - protein_bp_min[1])/(protein_bp_max[1] - protein_bp_min[1]))

            bp_score = math.sqrt(bp_x_dist**2 + bp_y_dist**2)

        if type(P1_mf) == type(None) or type(P2_mf) == type(None):
            mf_score = -10
        else:
            mf_x_dist = abs((P1_mf[0] - protein_mf_min[0]) / (protein_mf_max[0] - protein_mf_min[0]) \
                         - (P2_mf[0] - protein_mf_min[0]) / (protein_mf_max[0] - protein_mf_min[0]))

            mf_y_dist = abs((P1_mf[1] - protein_mf_min[1]) / (protein_mf_max[1] - protein_mf_min[1]) \
                         - (P2_mf[1] - protein_mf_min[1]) / (protein_mf_max[1] - protein_mf_min[1]))

            mf_score = math.sqrt(mf_x_dist ** 2 + mf_y_dist ** 2)

        if bp_score == -10 and mf_score == -10:
            continue

        else:

            if bp_score == -10:
                bp_score = 0

            if mf_score == -10:
                mf_score = 0

            proetin_similarity_dict[(P1, P2)]= max(bp_score, mf_score)
            similarity_distribution[int(max(bp_score, mf_score) * 10)] += 1

    return similarity_distribution, proetin_similarity_dict

def SAVE_file(protein_similarity_dict):
    file = open('distance_method.txt', 'w')

    for k, v in zip(protein_similarity_dict.keys(), protein_similarity_dict.values()):
        file.write(str(k) + '\t' + str(v) + '\n')

    file.close()

if '__main__' == __name__:
    file_name = ['ontology.obo', 'goa_human.gaf', 'PPI.txt', 'PPI_Similarity.txt']

    # Similarity CSV
    # PPI_Similarity_dict = PPI_Similar_CSVopen(file_name[3])

    # Ontology
    OBO_dict = OBO_fileopen(file_name[0])
    bp_root, mf_root = FIND_root(OBO_dict)
    FIND_error_relationship(OBO_dict)
    print('done, (ontology)')

    OBO_depth_dict = {
        'BP': dict(),
        'MF': dict()
    }

    FIND_depth(OBO_dict, OBO_depth_dict, bp_root, BP)
    FIND_depth(OBO_dict, OBO_depth_dict, mf_root, MF)

    BP_largest_len = max(OBO_depth_dict[BP].values())
    MF_largest_len = max(OBO_depth_dict[MF].values())

    # Depth Distribution
    # plt.figure(0)
    # plt.hist(OBO_depth_dict[BP].values())
    # plt.title("BP")
    #
    # plt.figure(1)
    # plt.hist(OBO_depth_dict[MF].values())
    # plt.title("MF")
    #
    # plt.show()

    # Annotation / Dict[protein] = term
    protein_dict = Anotation_PPI_fileopen(file_name[1], file_name[2])
    print('done, (annotation)')

    # Dict[BP/MF][term] = ancestor set
    ancestor_OBO_dict = FIND_OBO_ancestor(OBO_dict)

    # Dict[BP/MF][term] = protein
    OBO_term_protein_dict = FIND_OBO_term_protein(protein_dict, ancestor_OBO_dict)
    print('done, (protein by term)')

    OBO_term_vector = Term_vectorization(OBO_depth_dict, ancestor_OBO_dict, OBO_term_protein_dict)
    print('done, (term vector)')

    Protein_vector_dict = Protein_vectorization(file_name[2], protein_dict, OBO_term_vector, BP_largest_len, MF_largest_len)
    print('done, (protein vector)')

    tsne_result_bp, tsne_result_mf, bp_target, mf_target = Protein_Visualization(Protein_vector_dict)
    print('done, (tsne)')

    similarity_distribution, protein_similarity_dict = Distance_based_method(file_name[2], tsne_result_bp, tsne_result_mf, bp_target, mf_target)
    print('done, (similarity)')

    SAVE_file(protein_similarity_dict)
    print('done,(save)')

    plt.bar(similarity_distribution.keys(), similarity_distribution.values())
    plt.xticks(np.arange(0, 11, 1), labels=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    plt.savefig('Project.png')