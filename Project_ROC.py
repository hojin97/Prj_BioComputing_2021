import datetime
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from openpyxl import load_workbook

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

def Parsing(ExcelFileName):
    workbook=load_workbook(ExcelFileName,data_only=True)
    myexcel=workbook['Sheet1']

    NodeSimDict=dict()
    AnnoSimDict=dict()
    EdgeSimDict=dict()

    cnt=0

    for row in myexcel.rows:
        if cnt == 0:
            cnt+=1
            continue

        aa= row[0].value.split(',')
        aa[0]=aa[0].replace('(','')
        Gene1=aa[0].replace('\'','')
        aa[1]=aa[1].replace(')','')
        aa[1]=aa[1].replace(' ','')
        Gene2=aa[1].replace('\'', '')

        myKey=(Gene1, Gene2)

        NodeSimDict[myKey] = row[1].value
        AnnoSimDict[myKey] = row[2].value
        EdgeSimDict[myKey] = row[3].value


    return NodeSimDict, AnnoSimDict, EdgeSimDict

def Parsing_txt(txtName):
    textFile=open(txtName, 'r')
    returnDict=dict()

    while True:
        ppi_sim=textFile.readline()

        if not ppi_sim: break

        ppi_simlist= ppi_sim.split()

        ppi_simlist[0]=ppi_simlist[0].replace('(','')
        ppi_simlist[0]=ppi_simlist[0].replace(',','')
        Gene1=ppi_simlist[0].replace('\'','')

        ppi_simlist[1]=ppi_simlist[1].replace(')','')
        Gene2=ppi_simlist[1].replace('\'','')

        similarity=ppi_simlist[2]

        returnDict[(Gene1,Gene2)]=similarity

    return returnDict

def SortingSim(SimDict):

    templist=sorted(SimDict.items(), key=lambda x: x[1], reverse=True)
    ReturnDict=dict()

    for info in templist:
        ReturnDict[info[0]]=info[1]

    return ReturnDict

def MarkingGroundTruth(SimDict,Complexity):
    PPIKeyList=SimDict.keys()
    tempDict=dict()
    returnDict=dict()

    for PPI in PPIKeyList:
        tempDict['Score'] = SimDict[PPI]
        tempDict['GT'] = 'F'

        returnDict[PPI] = tempDict.copy()
        tempDict.clear()

        for Complex in Complexity:

            if PPI[0] in Complex and PPI[1] in Complex:

                tempDict['Score'] = SimDict[PPI]
                tempDict['GT'] = 'T'

                returnDict[PPI] = tempDict.copy()
                tempDict.clear()
                break

    return returnDict

def ReadComplexity(ComplexProtein):
    ComplexProteinFile = open(ComplexProtein, 'r')
    templist = []  # {protein : interactive protein}

    while True:
        ProteinComplexity = ComplexProteinFile.readline()

        if not ProteinComplexity:
            break

        ppilist = ProteinComplexity.split()
        templist.append(ppilist)

    return templist

def plotROC(MarkedSimDict,title):
    total=len(MarkedSimDict)
    TP=0
    FP=0
    TPlist=[]
    FPlist=[]

    TPlist.append(TP)
    FPlist.append(FP)

    AUC=0

    groundPKey = [k for k, v in MarkedSimDict.items() if v['GT'] == 'T']
    groundP=len(groundPKey)
    groundN= total-groundP

    for ppiInfo in MarkedSimDict.values():
        if ppiInfo['GT'] == 'T':
            TP+=1/groundP
            TPlist.append(TP)
            FPlist.append(FP)

        if ppiInfo['GT'] == 'F':
            FP+=1/groundN
            TPlist.append(TP)
            FPlist.append(FP)

            AUC=AUC+ TP*1/groundN

    y=np.array(TPlist)
    x=np.array(FPlist)

    plt.figure()
    plt.title(title)
    plt.plot(x,y,marker='o')
    plt.show()
    return AUC

def plotHistogram(SimDict,title):
    val=list(SimDict.values())
    # plt.figure(figsize=(8, 4))
    plt.figure()
    plt.title("Semantic Similarity - "+title)
    plt.hist(val, range=(0, 1), bins=10)
    plt.show()


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

    excelName = "PPI_Similarity.xlsx"
    ComplexProteinText="Complex_Protein.txt"
    txtName="distance_method_10vec.txt"
    introduce("Project", '*')

    Complexity=ReadComplexity(ComplexProteinText)
    NodeSimDict, AnnoSimDict, EdgeSimDict = Parsing(excelName)
    DistSimDict = Parsing_txt(txtName)
    print("Parsing Over")

    # plot Histogram
    plotHistogram(NodeSimDict, "Node-Based")
    plotHistogram(AnnoSimDict, "Annotation-Based (Best Match)")
    plotHistogram(EdgeSimDict, "Edge-Based (Best Match)")
    plotHistogram(DistSimDict, "Distance-Based")

    # Sort
    SortedNodeSimDict = SortingSim(NodeSimDict)
    SortedAnnoSimDict = SortingSim(AnnoSimDict)
    SortedEdgeSimDict = SortingSim(EdgeSimDict)
    SortedDistSimDict = SortingSim(DistSimDict)
    print("Sort Over")

    # Marking Ground Truth
    Marked_NodeSimDict = MarkingGroundTruth(SortedNodeSimDict,Complexity)
    Marked_AnnoSimDict = MarkingGroundTruth(SortedAnnoSimDict,Complexity)
    Marked_EdgeSimDict = MarkingGroundTruth(SortedEdgeSimDict,Complexity)
    Marked_DistSimDict = MarkingGroundTruth(SortedDistSimDict,Complexity)
    print("Marking Finished")

    # Plot ROC Curve and Compute Area Under the ROC
    NodeAUC=plotROC(Marked_NodeSimDict,'Node-Based')
    AnnoAUC=plotROC(Marked_AnnoSimDict,'Annotation-Based')
    EdgeAUC=plotROC(Marked_EdgeSimDict,'Edge-Based')
    DistAUC=plotROC(Marked_DistSimDict,'Distance-Based')
    print(f"Node-Based:{NodeAUC}\nAnnotation-Based:{AnnoAUC}\n\
          Edge-Based:{EdgeAUC}\nDistance-Based:{DistAUC}")
    # print(DistAUC)
