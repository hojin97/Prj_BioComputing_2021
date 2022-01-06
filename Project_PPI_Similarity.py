import datetime
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from itertools import product

def introduce(Topic, Num, sub=0):
    print("\n----------------------------------------------------------------")
    if (sub != 0):
        print(f"Bio Computing Assignment #{Num}-{sub}\nTopic:[ {Topic} ]")
    else:
        print(f"Bio Computing Assignment #{Num}\nTopic:[ {Topic} ]")

    print("Departure: Dept. of Computer TeleCommunication Engineerings")
    print("Student ID: 2016253020")
    print("Name: Ji-Hyeon Yoo")
    print("Student ID: 2016253004")
    print("Name: Ho-Jin Jung")
    print("----------------------------------------------------------------")

def readOntologyFile(oboFilename):
    # oboFilename="test.txt"
    structflag = False
    writeflag = False
    splitline = []

    dictlabel = "X"
    dictBP = dict()
    dictMF = dict()
    dictCC = dict()

    is_a = []
    part_of = []
    val = []

    duplicatedID = []

    obofile = open(oboFilename, 'r')
    while True:
        oboline = obofile.readline()

        if not oboline:
            break
        if oboline == "[Term]\n":
            structflag = True
            while structflag:
                line = obofile.readline()

                if ": " in line:
                    splitline = line.split(": ")
                    label = splitline[0]

                    # processing case by case
                    if label == "id":
                        # id = splitline[1].split(":")[1] # have to delete 'Go:'?
                        id = splitline[1]
                        id = id.rstrip()
                    elif label == "namespace":

                        if splitline[1] == "biological_process\n":
                            dictlabel = "BP"
                            writeflag = True
                        elif splitline[1] == "molecular_function\n":
                            dictlabel = "MF"
                            writeflag = True
                        elif splitline[1] == "cellular_component\n":
                            dictlabel = "CC"
                            writeflag = False

                    elif label == "is_obsolete":
                        if splitline[1] == "true\n":
                            is_a.clear()
                            part_of.clear()
                            break

                    elif label == "is_a":
                        # temp = splitline[1].split(":")[1]
                        temp = splitline[1]
                        is_a.append(temp.split(" ! ")[0])  # have to delete 'Go:'?

                    elif label == "relationship":
                        if "part_of " in splitline[1]:
                            # temp = splitline[1].split(":")[1]
                            # part_of.append(temp.split(" ! ")[0])
                            temp = splitline[1].split(" ")[1]
                            part_of.append(temp)
                    else:
                        continue

                else:  # read for data before starting line("[Term]")
                    structflag = False  # line == '\n'

                    tempset = set()
                    for parent in is_a:
                        val.append(parent)
                        tempset.add(parent)

                    for parent1 in part_of:
                        paraset = set()
                        paraset.add(parent1)
                        if writeflag:
                            if tempset.intersection(
                                    paraset) == set():  # pass 'part_of' where is it duplicated with 'is_a'
                                val.append(parent1)
                            else:
                                duplicatedID.append(id)

                    if dictlabel == "BP":
                        dictBP[id] = val.copy()
                    elif dictlabel == "MF":
                        dictMF[id] = val.copy()
                    else:  # CC
                        dictCC[id] = val.copy()

                    is_a.clear()
                    part_of.clear()
                    val.clear()
                    tempset.clear()

        else:
            continue

    return dictBP, dictMF, duplicatedID

def findError(dictBP, dictMF):
    errorCnt = 0

    MFkeylist = list(dictMF.keys())
    BPkeylist = list(dictBP.keys())

    for key, val in dictBP.items():  # dictionary's value : list type
        templist = val.copy()
        for BPid in val:  # for list 0 to the end

            if BPid in MFkeylist:
                # dictBP[key].remove(BPid)
                templist.remove(BPid)
                errorCnt = errorCnt + 1
        dictBP[key] = templist

    for key2, val2 in dictMF.items():
        templist = val2.copy()
        for MFid in val2:
            if MFid in BPkeylist:
                # dictMF[key2].remove(MFid)
                templist.remove(MFid)
                errorCnt = errorCnt + 1
        dictMF[key2] = templist

    return errorCnt

def findRoot(dictBP,dictMF):
    BProotList=[]
    MFrootList=[]

    for bpkey, bpval in dictBP.items():
        if len(bpval)==0:
            BProotList.append(bpkey)

    for mfkey,mfval in dictMF.items():
        if len(mfval)==0:
            MFrootList.append(mfkey)

    if len(BProotList)==1:
        if len(MFrootList)==1:
            return BProotList[0], MFrootList[0]

    print("Not the only one root")

def findDepth(dictBP,BProot,dictMF,MFroot):
    MFkeylist = list(dictMF.keys())
    BPkeylist = list(dictBP.keys())

    BPdepthDict=dict()
    MFdepthDict=dict()


    for id in BPkeylist:
        path_length = 0
        parents= dictBP[id]
        while True:
            if id == BProot:
                BPdepthDict[id]=0
                break
            elif BProot in parents:
                path_length=path_length+1
                BPdepthDict[id]=path_length
                break

            else:
                nextParents=set()
                path_length=path_length+1
                for parent in parents:
                    nextParents.update(dictBP[parent])

                parents=list(nextParents)

    for id in MFkeylist:
        path_length = 0
        parents= dictMF[id]
        while True:
            if id == MFroot:
                MFdepthDict[id]=0
                break
            elif MFroot in parents:
                path_length=path_length+1
                MFdepthDict[id]=path_length
                break

            else:
                nextParents=set()
                path_length=path_length+1
                for parent in parents:
                    nextParents.update(dictMF[parent])

                parents=list(nextParents)

    return BPdepthDict, MFdepthDict

def Distance_Edge(Ontologydict,ChildTerm,ParentTerm):

    path_length = 0
    parents = Ontologydict[ChildTerm]

    while True:
        if ChildTerm == ParentTerm:
            break

        if ParentTerm in parents:
            path_length = path_length + 1
            break

        else:
            nextParents = set()
            path_length = path_length + 1
            for parent in parents:
                nextParents.update(Ontologydict[parent])

            parents = list(nextParents.copy())


    return path_length

def readAnnotationFile(gafFilename):
    BPAnnotationDict = dict()
    MFAnnotationDict = dict()

    # initialization

    gafFile = open(gafFilename, 'r')
    while True:
        gafline = gafFile.readline()

        if not gafline:
            break
        # Comment -> Continue
        if gafline == "\n":
            continue
        if gafline[0] == '!':
            continue
        else:
            gaflist = gafline.strip().split('\t')

            # 4:qualified(except 'Not' ), 5: gene symbol 7: Evidence code(except IEA 9: Namespace //-1 for indexing
            # print(f"Qualified: {gaflist[3]}\tgene symbole: {gaflist[4]}\tEvidence: {gaflist[6]}\tNamespace: {gaflist[8]}")

            if "NOT" in gaflist[3] or gaflist[6] == "IEA" or gaflist[8] == 'C':
                continue
            else:
                if gaflist[8] == "P":  # BP 2번: 프로틴, 4번 프로틴 구성 텀
                    if gaflist[2] in BPAnnotationDict:
                        if gaflist[4] not in BPAnnotationDict[gaflist[2]]:  # eliminate duplicated term
                            BPAnnotationDict.get(gaflist[2]).append(gaflist[4])
                    else:
                        BPAnnotationDict[gaflist[2]] = [gaflist[4]]

                elif gaflist[8] == "F":  # MF
                    if gaflist[2] in MFAnnotationDict:
                        if gaflist[4] not in MFAnnotationDict[gaflist[2]]:  # eliminate duplicated term
                            MFAnnotationDict.get(gaflist[2]).append(gaflist[4])
                    else:
                        MFAnnotationDict[gaflist[2]] = [gaflist[4]]

    # print(BPAnnotationDict)
    # print(MFAnnotationDict)
    return BPAnnotationDict, MFAnnotationDict

def readPPI(ppiFilename):
    ppiFile = open(ppiFilename, 'r')
    returnppi = []  # {protein : interactive protein}
    #proteinSet = set()
    while True:
        ppiline = ppiFile.readline()

        if not ppiline: break
        ppilist = ppiline.split()
        #proteinSet.update(ppilist)
        returnppi.append(ppilist)

    return returnppi #, proteinSet

def exceptPPI(ppi, BPAnnotationDict, MFAnnotationDict):
    BPKeyList = list(BPAnnotationDict.keys())
    MFKeyList = list(MFAnnotationDict.keys())

    i = 0
    for interaction in ppi:
        #print(interaction)
        ## A-B interaction : (A is not in BP and MF) && (B is not in BP and MF)
        if (interaction[0] not in BPKeyList and interaction[0] not in MFKeyList) or \
           (interaction[1] not in BPKeyList and interaction[1] not in MFKeyList):
            ppi[i] = []
        # ## A-B interaction : only one side is element of BP or MF
        # elif (interaction[0] not in BPKeyList and interaction[0] not in MFKeyList) and \
        #      (interaction[1] in BPKeyList or interaction[1] in MFKeyList):
        #       ppi[i] = [interaction[1]]
        #
        # elif (interaction[0] in BPKeyList or interaction[0] in MFKeyList) and \
        #      (interaction[1] not in BPKeyList and interaction[1] not in MFKeyList):
        #       ppi[i] = [interaction[0]]

        ## A-B interaction : BP - MF or MF - BP
        if interaction[0] in BPKeyList and interaction[0] not in MFKeyList and\
           interaction[1] not in BPKeyList and interaction[1] in MFKeyList:
            ppi[i] = []

        if interaction[0] not in BPKeyList and interaction[0] in MFKeyList and\
           interaction[1] in BPKeyList and interaction[1] not in MFKeyList:
            ppi[i] = []

        i = i + 1

    returnppi = list(filter(([]).__ne__, ppi))

    return returnppi

def TraceAncester_Term(TermDict):
    termKeyList = list(TermDict.keys())
    parentsDict = dict()

    for term in termKeyList:
        temp = TermDict[term].copy()
        temp=set(temp)

        parents = TermDict[term].copy()
        parents=set(parents)

        while temp:
            check = temp.pop()
            temp.update(TermDict[check])
            parents.update(TermDict[check])

        parentsDict[term] = list(parents.copy())

    return parentsDict

def TraceAncester_Protein(AnnotationDict,AncesterDict):
    proteinAncesterDict=dict()
    for key, terms in AnnotationDict.items():
        tempset=set()
        for term in terms:
            tempset.update(AncesterDict[term])

        proteinAncesterDict[key]=list(tempset.copy())

    return proteinAncesterDict

def ProteinCode(protein,BPKeyList,MFKeyList):
    # BPCode='X'
    # MFCode='X'
    if protein not in BPKeyList:
        BPCode='0'
    else:
        BPCode='1'

    if protein not in MFKeyList:
        MFCode='0'
    else:
        MFCode='1'

    Code=MFCode+BPCode

    return Code

def ComputeJaccard(ppi,AnnotationDict,ProteinAncesterDict):
    AnnotationKeyList=list(AnnotationDict.keys())

    result=dict()

    for interaction in ppi:
        Gene1=interaction[0]
        Gene2=interaction[1]

        if Gene1 not in AnnotationKeyList or Gene2 not in AnnotationKeyList:
            result[(Gene1, Gene2)] = None
            #print("없어야되는데")
            continue
        else:
            pass

        tempA = ProteinAncesterDict[interaction[0]].copy()
        tempA = set(tempA)
        tempB = ProteinAncesterDict[interaction[1]].copy()
        tempB = set(tempB)

        Numer = len(tempA.intersection(tempB))
        Denom = len(tempA.union(tempB))

        if Denom == 0:
            Score = 0
        else:
            Score = Numer / Denom

        # result.append(Score)
        result[(Gene1, Gene2)] = Score

    return result

def TraceAnnotation(AnnotationDict,AncesterDict):
    ContentDict=dict()

    for Term, Ancesters in AncesterDict.items():
        GeneList=[k for k, v in AnnotationDict.items() if Term in v]
        tempSet=set(GeneList.copy())
        try:
            ContentDict[Term].update(tempSet.copy()) # recording all genes as value which contains a Term{Term: [gene1,gene2,...]}
            for Ancester in Ancesters:  # recording all genes with a term's all parents
                ContentDict[Ancester].update(tempSet.copy())
        except KeyError:
            ContentDict[Term] = tempSet.copy()
            for Ancester in Ancesters:
                try:
                    ContentDict[Ancester].update(tempSet.copy())
                except KeyError:
                    ContentDict[Ancester] = tempSet.copy()

        GeneList.clear()
        tempSet.clear()

    return ContentDict

def ComputeInformationContent(root,ContentDict,TermDict):
    InformationContentDict=dict()

    Denom = len(ContentDict[root])  # ''has contained with root's value
    if '' in ContentDict[root]:
        Denom = Denom-1

    for Term in TermDict.keys():
        Num = len(ContentDict[Term])

        if '' in ContentDict[Term]:
            Num = Num - 1

        InformationContent = Num/Denom
        if InformationContent != 0:
            InformationContentDict[Term] = -np.log2(InformationContent)    # commercial log?  natural log?
        else:
            InformationContentDict[Term] = -1

    #print(InformationContentDict)

    return InformationContentDict   # value: 0.0 ==> ContentDict[Term]=set{} data

def SpecificAncester_Edge(Term1,Term2,AncesterDict,DepthDict):
    if Term1=='GO:0008150' or Term2== 'GO:0008150': # Root ID of BP
        return 'Root'
    elif Term1=='GO:0003674' or Term2=='GO:0003674':  # Root ID of MF
        return 'Root'

    Ancester1=AncesterDict[Term1]
    Ancester2=AncesterDict[Term2]

    AncesterSet1=set(Ancester1)
    AncesterSet2=set(Ancester2)

    tempdict=dict()
    Candidate_Ancester=AncesterSet1.intersection(AncesterSet2)

    for Ancester in Candidate_Ancester:
        tempdict[Ancester]=DepthDict[Ancester]

    if len(tempdict)==0: Most_Specific_Common_Ancester='Root'
    else: Most_Specific_Common_Ancester=max(tempdict,key=tempdict.get)

    return Most_Specific_Common_Ancester

def SpecificAncester_Annotation(Term1,Term2,AncesterDict,InformationContentDict):
    # if term1 or term 2 == root > score == 0
    if Term1=='GO:0008150' or Term2== 'GO:0008150': # Root ID of BP
        return 'Root'
    elif Term1=='GO:0003674' or Term2=='GO:0003674':  # Root ID of MF
        return 'Root'

    Ancester1=AncesterDict[Term1]
    Ancester2=AncesterDict[Term2]

    AncesterSet1=set(Ancester1)
    AncesterSet2=set(Ancester2)

    tempdict=dict()
    Candidate_Ancester=AncesterSet1.intersection(AncesterSet2)

    for Ancester in Candidate_Ancester:
        tempdict[Ancester]=InformationContentDict[Ancester]

    if len(tempdict) == 0: Most_Specific_Common_Ancester='Root'
    else: Most_Specific_Common_Ancester=max(tempdict,key=tempdict.get)

    return Most_Specific_Common_Ancester

def SemanticSimilarity_Annotation(Term1,Term2,InformationContentDict,AncesterDict):
    Most_Specific_Common_Ancestor=SpecificAncester_Annotation(Term1,Term2,AncesterDict,InformationContentDict)

    if Most_Specific_Common_Ancestor=='Root':
        return 0
    if InformationContentDict[Term1] == -1 or InformationContentDict[Term2] == -1:
        return -1

    # if InformationContentDict[Most_Specific_Common_Ancestor] == -0: InformationContentDict[Most_Specific_Common_Ancestor] = 0
    # if int(InformationContentDict[Term1]) == -0:  InformationContentDict[Term1] = 0
    # if int(InformationContentDict[Term2]) == -0:  InformationContentDict[Term2] = 0

    Num= 2 * InformationContentDict[Most_Specific_Common_Ancestor]
    Denom= InformationContentDict[Term1] + (InformationContentDict[Term2])

    if Denom == 0: Similarity = -1
    else: Similarity = Num/Denom

    return Similarity

def SemanticSimilarity_Edge(Term1,Term2,OntologyDict,AncesterDict,DepthDict):
    Most_Specific_Common_Ancestor=SpecificAncester_Edge(Term1,Term2,AncesterDict,DepthDict)

    if Most_Specific_Common_Ancestor=='Root':
        return 0
    if DepthDict[Term1] == 0 or DepthDict[Term2] == 0:
        return 0

    # if InformationContentDict[Most_Specific_Common_Ancestor] == -0: InformationContentDict[Most_Specific_Common_Ancestor] = 0
    # if int(InformationContentDict[Term1]) == -0:  InformationContentDict[Term1] = 0
    # if int(InformationContentDict[Term2]) == -0:  InformationContentDict[Term2] = 0

    Num= 2*DepthDict[Most_Specific_Common_Ancestor]
    len1=Distance_Edge(OntologyDict,Term1,Most_Specific_Common_Ancestor)
    len2=Distance_Edge(OntologyDict,Term2,Most_Specific_Common_Ancestor)
    Denom= len1+len2+Num

    if Denom == 0: Similarity = 0
    else: Similarity = Num/Denom

    return Similarity

def Best_Match_Averaging(ppi,AnnotationDict,InformationContentDict,AncesterDict):
    AnnotationKeyList=AnnotationDict.keys()
    MatchDict=dict()

    for interaction in ppi:
        Gene1=interaction[0]
        Gene2=interaction[1]

        if Gene1 not in AnnotationKeyList or Gene2 not in AnnotationKeyList:
            MatchDict[(Gene1, Gene2)] = None
            #print("없어야되는데")
            continue
        else:
            pass

        Gene1Terms=AnnotationDict[Gene1]
        Gene2Terms=AnnotationDict[Gene2]

        maxlist=[]
        for Gene1Term in Gene1Terms:
            tmp=[Gene1Term]
            templist=[tmp,Gene2Terms]
            Combinations = list(product(*templist))

            score=[]
            for Combination in Combinations:
                sim=SemanticSimilarity_Annotation(Combination[0],Combination[1],InformationContentDict,AncesterDict)
                if sim != -1:
                    score.append(sim)

            if len(score) != 0:
                maxlist.append(max(score))

        for Gene2Term in Gene2Terms:
            tmp=[Gene2Term]
            templist=[Gene1Terms,tmp]
            Combinations=list(product(*templist))

            score = []
            for Combination in Combinations:
                sim = SemanticSimilarity_Annotation(Combination[0], Combination[1], InformationContentDict,
                                                    AncesterDict)
                if sim!=-1:
                    score.append(sim)

            if len(score)!=0:
                maxlist.append(max(score))

        Denom = len(maxlist)

        if Denom == 0: MatchDict[(Gene1, Gene2)] = 0
        else:
            summation = sum(maxlist)
            MatchDict[(Gene1, Gene2)] = summation/Denom

        maxlist.clear()

    # print(MatchDict)
    return MatchDict

def Best_Match_Averaging_Edge(ppi,OntologyDict,AnnotationDict,AncesterDict,DepthDict):
    AnnotationKeyList=AnnotationDict.keys()
    MatchDict=dict()

    for interaction in ppi:
        Gene1=interaction[0]
        Gene2=interaction[1]

        if Gene1 not in AnnotationKeyList or Gene2 not in AnnotationKeyList:
            MatchDict[(Gene1, Gene2)] = None
            continue
        else:
            pass

        Gene1Terms=AnnotationDict[Gene1]
        Gene2Terms=AnnotationDict[Gene2]

        maxlist=[]
        for Gene1Term in Gene1Terms:
            tmp=[Gene1Term]
            templist=[tmp,Gene2Terms]
            Combinations = list(product(*templist))

            score=[]
            for Combination in Combinations:
                sim=SemanticSimilarity_Edge(Combination[0],Combination[1],OntologyDict,AncesterDict,DepthDict)
                if sim != -1:
                    score.append(sim)

            if len(score) != 0:
                maxlist.append(max(score))

        for Gene2Term in Gene2Terms:
            tmp=[Gene2Term]
            templist=[Gene1Terms,tmp]
            Combinations=list(product(*templist))

            score = []
            for Combination in Combinations:
                sim = SemanticSimilarity_Edge(Combination[0], Combination[1],OntologyDict, AncesterDict, DepthDict)
                if sim!=-1:
                    score.append(sim)

            if len(score)!=0:
                maxlist.append(max(score))

        Denom = len(maxlist)

        if Denom == 0: MatchDict[(Gene1, Gene2)] = 0
        else:
            summation = sum(maxlist)
            MatchDict[(Gene1, Gene2)] = summation/Denom

        maxlist.clear()

    # print(MatchDict)
    return MatchDict

def CompareBP_MF_Similarity(BPSimilarityDict,MFSimilarityDict):
    scoreCandidate=[]
    result=dict()

    for ppi, BPval in BPSimilarityDict.items():
        MFval=MFSimilarityDict[ppi]
        if BPval != None and MFval != None:
            scoreCandidate.append(BPval)
            scoreCandidate.append(MFval)
            result[ppi] = max(scoreCandidate)
        elif BPval != None and MFval == None:
            result[ppi]=BPval
        elif BPval ==None and MFval !=None:
            result[ppi]= MFval

        scoreCandidate.clear()
    #print(result)
    return result

def GatheringDictionaries(ppi,NodeDict,AnnotationDict,EdgeDict):
    result=dict()

    for interaction in NodeDict.keys():

        Gene1 = interaction[0]
        Gene2 = interaction[1]

        myKey = (Gene1, Gene2)

        # if Gene1 not in AnnotationKeyList or Gene2 not in AnnotationKeyList:
        #     result[myKey]=None
        #     continue
        # else:
        #     pass
        tempdict = dict()

        tempdict['Node']=NodeDict[myKey]
        tempdict['Annotation']=AnnotationDict[myKey]
        tempdict['Edge']=EdgeDict[myKey]

        result[myKey]=tempdict.copy()

    #print(result)
    return result

def plotHistogram(val):
    #plt.figure(figsize=(8, 4))
    figBP=plt.figure()
    plt.title("Semantic Similarity - Edge Based, Best-Match Averaging")
    plt.hist(val,range=(0,1),bins=10)
    plt.show()

def exportExcel(exportDict):
    N_A_E_list=list(exportDict.values())
    nodeval=[]
    annoval=[]
    edgeval=[]
    for informations in N_A_E_list:
        nodeval.append(informations['Node'])
        annoval.append(informations['Annotation'])
        edgeval.append(informations['Edge'])

    excelFormat=pd.DataFrame({
        'PPI': list(exportDict.keys()),
        'Node_Based': nodeval,
        'Annotation_Based': annoval,
        'Edge_Based': edgeval
    })
    excelFormat.to_excel('PPI_Similarity.xlsx',index=False)

if __name__ == '__main__':  # ontology , ppi , annotation
    # try:
    #     if os.stat(sys.argv[1]).st_size == 0 or os.stat(sys.argv[2]).st_size == 0 or os.stat(sys.argv[3]).st_size == 0:
    #         print('No string found')
    #
    #     oboFilename = sys.argv[1]
    #     gafFilename = sys.argv[2]
    #     ppiFilename = sys.argv[3]
    # except FileNotFoundError:
    #     print('No input file')

    oboFilename = "ontology.obo"
    gafFilename = "goa_human.gaf"
    ppiFilename = "ppi.txt"

    introduce("Project", '*')

    # Ontology Parsing
    dictBP, dictMF, duplicatedID = readOntologyFile(oboFilename)
    errorCnt = findError(dictBP, dictMF)
    print("Ontology Parsing Over")

    # Annotation Parsing
    BPAnnotationDict, MFAnnotationDict = readAnnotationFile(gafFilename)
    print("GAF Parsing Over")

    # PPI Parsing
    ppi = readPPI(ppiFilename)  # [[protein1, protein2],[protein1,protein3]...], {protein1, protein2,....}
    #ppi = exceptPPI(ppi, BPAnnotationDict, MFAnnotationDict)
    print("PPI Parsing Over")
    print("Parsing End")

    # Trace ancester of Term
    BPAncesterDict = TraceAncester_Term(dictBP)
    MFAncesterDict = TraceAncester_Term(dictMF)
    print("Ancester Trace Over (Term)")

    # Trace ancester of Protein
    BPProteinAncesterDict=TraceAncester_Protein(BPAnnotationDict,BPAncesterDict)
    MFProteinAncesterDict=TraceAncester_Protein(MFAnnotationDict,MFAncesterDict)
    print("Ancester Trace Over (Protein)")

    # Trace all Term - Annotation Relation {Term: [Annotation,...]} All Parents contains child's annotation
    BPContentDict=TraceAnnotation(BPAnnotationDict,BPAncesterDict)
    MFContentDict=TraceAnnotation(MFAnnotationDict,MFAncesterDict)
    #print(MFContentDict)
    print("Annotation Trace Over")

    #depth
    rootBP, rootMF = findRoot(dictBP, dictMF)
    BPdepthDict, MFdepthDict = findDepth(dictBP, rootBP, dictMF, rootMF)
    BPICDict = ComputeInformationContent(rootBP, BPContentDict, dictBP)
    MFICDict = ComputeInformationContent(rootMF, MFContentDict, dictMF)
    print("PreProcessing has finished")

    # Edge Based Method
    print("Wait...")
    start = datetime.datetime.now()

    # Node Based Result
    BPMatchDict_Node = ComputeJaccard(ppi, BPAnnotationDict, BPProteinAncesterDict)
    MFMatchDict_Node = ComputeJaccard(ppi, MFAnnotationDict, MFProteinAncesterDict)
    result_Node= CompareBP_MF_Similarity(BPMatchDict_Node,MFMatchDict_Node)

    # Annotation Based Result
    BPMatchDict_Annotation = Best_Match_Averaging(ppi, BPAnnotationDict, BPICDict, BPAncesterDict)  # {(gene1,gene2):similarity}
    MFMatchDict_Annotation = Best_Match_Averaging(ppi, MFAnnotationDict, MFICDict, MFAncesterDict)
    result_Annotation = CompareBP_MF_Similarity(BPMatchDict_Annotation, MFMatchDict_Annotation)

    # Edge Based Result
    BPMatchDict_Edge=Best_Match_Averaging_Edge(ppi,dictBP,BPAnnotationDict,BPAncesterDict,BPdepthDict)  # {(gene1,gene2):similarity}
    MFMatchDict_Edge=Best_Match_Averaging_Edge(ppi,dictMF,MFAnnotationDict,MFAncesterDict,MFdepthDict)
    result_Edge=CompareBP_MF_Similarity(BPMatchDict_Edge,MFMatchDict_Edge)

    end = datetime.datetime.now()
    # print(result)
    print(f"\nElapsed time: {(end - start)} micro-sec")


    # plot
    # plotHistogram(list(result_Node.values()))
    # plotHistogram(list(result_Annotation.values()))
    # plotHistogram(list(result_Edge.values()))

    # export to excel
    Similarity_Result_Dict=GatheringDictionaries(ppi,result_Node,result_Annotation,result_Edge)
    exportExcel(Similarity_Result_Dict)