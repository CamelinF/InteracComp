from pulp import *
import os,sys,math,time,random,bisect,gc
from skmine.datasets.fimi import fetch_file
def dataTransToDataItem(data):
    dataI={}
    for t in range(len(data)):
        for i in data[t]:
            if i in dataI:
                dataI[i].add(t)
            else:
                dataI[i]={t}
    return dataI
# OBTENTION DU DECODAGE de dataFile 
# FALSE => original to  Krimp,  True => Krimp to original
def getDecodage(dataFile,sens):
    f=open(f"datasets/{dataFile}.db")
    s=f.readline()
    while s[0]!="a" and s[1] !="b":
        s=f.readline()
    alphabet=[]
    splitS=s.split()
    for i in range(1,len(splitS)):
        alphabet.append(int(splitS[i]))

    while s[0]!="i" and s[1] !="t":
        s=f.readline()
    items=[]
    splitS=s.split()
    for i in range(1,len(splitS)):
        items.append(int(splitS[i]))
    f.close()
    dico={}
    for i in range(len(items)):
        if sens :
            dico[alphabet[i]]=items[i]
        else:
            dico[items[i]]=alphabet[i]
    return dico

def getCT(fileNamD,fileNamCT):
    dicoK= getDecodage(fileNamD,True)
    f=open(f"codetables/{fileNamCT}.ct")
    s=f.readline()
    s=f.readline()
    s=f.readline()
    splitS=s.split()
    CTD=[]
    CTK=[]
    while s!="":
        itemset=[]
        motif=[]
        for indIt in range(len(splitS)-1):
            item=splitS[indIt]
            itemset.append(int(item))
            motif.append(dicoK[int(item)])
        itemset.sort()
        motif.sort
        tmp=splitS[-1].split(",")[0]
        tmpp=int(tmp.split("(")[1])
        tmp=splitS[-1].split(",")[1]
        freq=int(tmp.split(")")[0])
        tmp=splitS[-1].split(",")[0]
        nbUtil=int(tmp.split("(")[1])
        if len(motif)>0:
            CTK.append([freq,len(itemset),itemset,nbUtil])
            CTD.append([freq,len(itemset),motif,nbUtil])
        s=f.readline()
        splitS=s.split()
    f.close()
    return CTK,CTD

def getCand(fileNamD,fileNamCT):
    dicoK= getDecodage(fileNamD,True)
    f=open(f"candidates/{fileNamCT}.isc")
    s=f.readline()
    s=f.readline()
    s=f.readline()
    splitS=s.split()
    CTD=[]
    CTK=[]
    while s!="":
        itemset=[]
        motif=[]
        for indIt in range(1,len(splitS)-1):
            item=splitS[indIt]
            itemset.append(int(item))
            motif.append(dicoK[int(item)])
        itemset.sort()
        motif.sort()
        tmp=splitS[-1].split("(")[1]
        freq=int(tmp.split(")")[0])
        if len(motif)>1:
            CTK.append([freq,len(itemset),itemset,freq])
            CTD.append([freq,len(itemset),motif,freq])
        s=f.readline()
        splitS=s.split()
    f.close()
    return CTK,CTD

def getFreq(d,m):
    couv=d[m[0]]
    for i in m:
        couv=couv.intersection(d[i])
        #print(f"{m} à uen couv de {len(couv)}")
    return couv,len(couv)




fileNameData="hepatitis"
data=fetch_file(f'datasets/{fileNameData}.dat',int_values=True)
dataI=dataTransToDataItem(data)
sumItemFreq=0
for i in dataI.keys():
    sumItemFreq+=len(dataI[i])

fileNameCT=f"{fileNameData}-c10f20k25"
CTK,CTD=getCT(fileNameData,fileNameCT)
# fileNameCand="pumsb-c10f75"
# CTK,CTD=getCand(fileNameData,fileNameCand)

dicoCouv=dict()
compteurFreq=0
comptFreqK=0
comptSinglK=0
ajoutSingle=False
nbMotifs=0
for i in range(len(CTD)):
    if CTD[i][1]>1:
        comptFreqK+=CTD[i][3]
        nbMotifs+=1
    else:
        comptSinglK+=CTD[i][3]
        ajoutSingle=True
    couv,freqM=getFreq(dataI,CTD[i][2])
    for j in couv:
        compteurFreq+=1
        if  j in dicoCouv:
            dicoCouv[j].append(i)
        else:
            dicoCouv[j]=[i]

print(f"nb motifs:{nbMotifs} et nb item:{len(CTD)-nbMotifs}")
print(f"nb util singleton Krimp={comptSinglK}")
print(f"nb util motifs Krimp={comptFreqK}")
print(f"nb utilisations Krimp ={comptFreqK+comptSinglK}")

print(f"nb utilisations max ={compteurFreq}")

#print(f"clé dicoCouv : {sorted(dicoCouv.keys())}")
cptCouv=0
for i in range(len(data)):
    if not(i in dicoCouv.keys()):
        cptCouv+=1
        # print(f"trans {i} n'est couverte par aucun motif : {data[i]} ")
print(f" il y a {cptCouv} sur {len(data)} de transactions non couvertes par un motif")


# Maximize score
prob=LpProblem("Problème_de_couverture_motifs", LpMaximize)
#Minimize nb utils
# prob=LpProblem("Problème_de_couverture_motifs", LpMinimize)
# couvTM = [[LpVariable(f"couvTM_{i}_{j}", cat='Binary') for j in range(len(CTD))] for i in range(len(data))]

# long=0
# for i in couvTM:
#     long+=len(i)
# print(f" long 1: {long}")
couvTM = {
    line: { motif: LpVariable(f"Var_{line}_{motif}", cat=LpBinary) for motif in motifs}
    for line, motifs in dicoCouv.items()
    
}
useM=[LpVariable(f"use_{i}" ,cat=LpBinary) for i in range(len(CTD)) ]
long=len(useM)
for i in couvTM:
    long+=len(couvTM[i])
# print(f" long 1: {long}")
print(f"{long} variables crées")

#Maxizime score
prob += (
    lpSum(couvTM[i][j] for i in dicoCouv for j in dicoCouv[i] if len(CTD[j][2])>1),
    "Somme_couverture",
)
#Minizime nombre util
# prob += (
#     lpSum(couvTM)
# )
# Nomnre de motifs différrents utilisés
# prob+=(
#     lpSum(useM[i] for i in range(len(CTD)) if len(CTD[j][2])==1)
# )
for i in range(len(data)):
    if i in dicoCouv:
        for j in data[i]:
            couvPossibles=[]
            for k in dicoCouv[i]:
                if j in CTD[k][2]:
                    couvPossibles.append((i,k))
            prob+=(lpSum([couvTM[x][y] for (x,y) in couvPossibles]) == 1,
            f"couverture_item_{j}_dans_transac_{i}",)
for i in range(len(CTD)):
    listSum=[]
    for j in dicoCouv.keys():
        if i in dicoCouv[j]:
            listSum.append((j,i))
    prob+=(lpSum([couvTM[x][y] for (x,y) in listSum])/CTD[i][0] <= useM[i],f"motif {i} utilise au moins une fois",)     

# for i in dicoCouv.keys():
#     prob+=(lpSum([couvTM[i][j] for j in dicoCouv[i]]) >= 1,f" transation {i} couverte par au moins 1 motif")

prob+=(lpSum(useM[i] for i in range(len(CTD))if len(CTD[i][2])>1) <=10,"limite du nombre de motifs dfferrents utilisé")
optionsG = [
    ("TimeLimit" , 60.0),
    ("WLSACCESSID", "146955fd-fb4c-409b-bcd0-7eb0472c8d4b"),
    ("WLSSECRET", "ae52ec3b-abb7-4028-9d66-285e40efdc7d"),
    ("LICENSEID", 2610563)
]
solver = GUROBI()
status = prob.solve(solver)
print(LpStatus[status])

print("analyse solution")
comptSingl=0
comptMotif=0
tailleCompLp=0
tailleEncodeur=0
nbUtilTot=value(prob.objective)
for i in range(len(CTD)):
    compt=0
    for j in dicoCouv.keys():
        if i in dicoCouv[j]:
            if value(couvTM[j][i])>0:
                compt+=1
    if compt>0:
        tailleCode=-1*math.log2(compt/nbUtilTot)
        tailleEncode=0
        for item in CTD[i][2]:
            tailleEncode+=-1*math.log2(len(dataI[item])/sumItemFreq)
        tailleEncodeur+=tailleEncode+tailleCode
        tailleCompLp+=compt*tailleCode
    if compt!=CTD[i][3]:
        if compt>CTD[i][3]:
            print(f" motif{i} utilisé plus de fois {compt} > {CTD[i][3]}")
        else:
            print(f" motif{i} utilisé moins de fois {compt} < {CTD[i][3]}")
    if CTD[i][1]>1:
        comptMotif+=compt 
    else:
        comptSingl+=compt 
print(f"nb util Singl={comptSingl} contre {comptSinglK}")
print(f"nb util motifs={comptMotif} contre {comptFreqK}")
print(f"nb util tot= {comptSingl+comptMotif} contre {comptSinglK+comptFreqK}")
print(f"taille du code compressé avec LP {tailleCompLp} ")
print(f"taille de l'encodeur avec LP {tailleEncodeur} ")
print(f"taille de la compression avec LP {tailleCompLp+tailleEncodeur} ")
print(f"valeur objectif = {value(prob.objective)}")