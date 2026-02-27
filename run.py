import os
import numpy as np
from skmine.datasets.fimi import fetch_file
from skmine.datasets.utils import describe
import pandas as pd
import matplotlib.pyplot as plt

from skmine.itemsets import LCM
from random import sample,random,choices,randint
import time 
couleur=['red','blue','green',"cyan","magenta","brown","purple","yellow","orange","pink","grey"]
dataset=["hepatitis","chess","mushroom","retail","splice1","eisen","pumsb","pumsb_star","connect","weatherAUS","twitter","dota"]
dataset=["hepatitis","chess","mushroom","retail","splice1","pumsb","pumsb_star","connect","weatherAUS","twitter"]
dataset=["hepatitis","chess","mushroom","splice1","connect"]
# dataset=["hepatitis","chess","mushroom","retail","splice1","pumsb","pumsb_star","connect","weatherAUS","twitter"]

# dataset=["hepatitis","chess","mushroom","pumsb_star","twitter"]

# dataset=["hepatitis","chess","mushroom","retail","splice1","pumsb","pumsb_star","connect","weatherAUS","twitter"]
# dataset=["hepatitis","chess","mushroom","retail","splice1","connect","weatherAUS","twitter"]
# dataset=["splice1","eisen","pumsb","pumsb_star","connect","weatherAUS","twitter"]
# dataset=["connect","weatherAUS","twitter"]
# dataset=["hepatitis","chess",'mushroom']
# dataset=["pumsb","pumsb_star"]
# dataset=["hepatitis"]


# treshold for 145k Frequent itemset
threshold={
    "hepatitis":48,
    "chess":2014,
    "mushroom":913,
    "retail":11,
    "splice1":92,
    "eisen":28,
    "pumsb":39200,
    "pumsb_star":16750,
    "connect":57370,
    "weatherAUS":540,
    "twitter":23050
}
# treshold for 300k Closed itemset
threshold={
    "hepatitis":30,
    "chess":1650,
    "mushroom":10,
    "retail":11,
    "splice1":75,
    "eisen":28,
    "pumsb":39200,
    "pumsb_star":16750,
    "connect":24000,
    "weatherAUS":540,
    "twitter":23050
}
nbIter=20
nbRun=10


# for i in dataset:
   
#     for j in range(nbRun):
#         print(i)
#         # nombre de transactions C=10,  nombre d'itération n=25, Fréquence dans les transactions échantillonné f= 0.2 
#         # nombre de motifs à présenter k=5 et nombre de motifs à garder d'une requête à l'autre l=1
#         os.system(f"python3 LCS-IH.py {i} 30 {nbIter} 0.5 5 1")
#         os.system(f"mkdir Res/{i}/Candidates/{j}")
#         os.system(f"mkdir Res/{i}/Codetables/{j}")
#         os.system(f"mkdir Res/{i}/GraphPrefs/{j}")
#         os.system(f"mkdir Res/{i}/Interacts")
#         os.system(f"mv Res/{i}/patternShown.txt Res/{i}/Interacts/patternShown{j}.txt")
#         os.system(f"mv Res/{i}/feedbacks.txt Res/{i}/Interacts/feedbacks{j}.txt")
#         os.system(f"mv Res/{i}/weightLearned.txt Res/{i}/Interacts/weightLearned{j}.txt")
#         for k in range(nbIter):
#             os.system(f"mv Res/{i}/Candidates/run{k}.isc  Res/{i}/Candidates/{j}/run{k}.isc")
#             os.system(f"mv Res/{i}/Codetables/run{k}.ct  Res/{i}/Codetables/{j}/run{k}.ct")
#             os.system(f"mv Res/{i}/GraphPrefs/run{k}.eps  Res/{i}/GraphPrefs/{j}/run{k}.eps")



def getSupp(d,m):
    # print(m)
    couv=d[int(m[0])]
    for i in m:
        couv=couv.intersection(d[int(i)])
        #print(f"{m} à uen couv de {len(couv)}")
    return couv,len(couv)

def get_Jaccard(covX,covY):
    temp = covX
    covX = covX.union(covY)
    covY = covY.intersection(temp)
    jaccard = 0
    if(len(covX) > 0):
        jaccard = len(covY)/len(covX)
    return jaccard

def isClosed(d,m,couv):
    couv=list(couv)
    closure=set(d[int(couv[0])])
    for i in couv:
        closure=closure.intersection(set(d[int(i)]))
    return len(closure)==len(m)
def getVecteurFromLine(line):
    splitS=line.split(",")
    vecteur=[float(splitS[0].split('[')[1])]
    for i in range(1,len(splitS)-1):
        vecteur.append(float(splitS[i]))
    vecteur.append(float(splitS[-1].split(']')[0]))
    return vecteur
# OBTENTION DU DECODAGE
# sens ==FALSE => original vers Krimp, sens True => Krimp vers original
def getDecodage(dataFile,sens):
    f=open(f"data/{dataFile}.db")
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

def getFreq(d,m):
    couv=d[m[0]]
    for i in m:
        couv=couv.intersection(d[i])
        #print(f"{m} à uen couv de {len(couv)}")
    return couv,len(couv)

def getPatternShown(d,nbItems,nbTransac,indRun,keysIndex):
    f=open(f"Res/{d}/Interacts/patternShown{indRun}.txt")
    motifs=[]
    s=f.readline()
    while s!="":
        splitS=s.split("[")[2].split(",")
        itemset=[]
        cpt=0
        for i in range(len(splitS)-2):
            if int(splitS[i]) == 1:
                cpt+=1

        # print(len(splitS))
        for i in range(nbItems-1):
            if int(splitS[i+nbTransac])==1:
                itemset.append(keysIndex[i])
        if int(splitS[nbItems+nbTransac-1][1])==1:
            itemset.append(keysIndex[nbItems-1])
        s=f.readline()
        tmpX=splitS[len(splitS)-1].split('\n')[0]
        # print(f"itemset ={itemset}, lenFeatures={len(splitS)} , cpt ={cpt}, -1 ={tmpX}, -2 ={splitS[len(splitS)-2]}")

        motifs.append(itemset)
    f.close()
    return motifs
tempsDetails=["Krimp","Learn","Req","Sample","Weight","Total"]




fileCumulRgt=open(f"Res/cumulRegret.txt","w")
fileCumulRgt.write(f"dataset;min;mavg;max\n")
plt.figure(figsize=(15, 10))
plt.xlabel("iterations")
plt.ylabel("regret")
plt.title(f"Regret Avg and Max ")
plt.tight_layout()

for d in dataset:
    print(f"Pour {d}")
    theta=threshold[f'{d}']
    dicoK=getDecodage(d,True)
    data=fetch_file(f'data/{d}.dat',int_values=True)
    dataI={}
    maxIndex=0
    for t in range(len(data)):
        maxIndex=max(maxIndex,t)
        for i in data[t]:
            if i in dataI:
                dataI[i].add(t)
            else:
                dataI[i]={t}
    infoD=describe(data)
    print(infoD)
    print(f" {len(dataI)},{len(data)}")
    freqItems=np.zeros(maxIndex)
    for i in range(infoD['n_items']):
        couv,supp=getSupp(dataI,[f"{dicoK[i]}"])
        freqItems[dicoK[i]-1]=supp/len(data)
    
    #Analyse temps
    # dictTemps=dict()
    # for t in tempsDetails:
    #     dictTemps[t]=[]
    # dictTemps["trueTotal"]=[]
    #     # Get temps
    # for t in tempsDetails:
    #     f=open(f"Res/{d}/temps{t}.txt")
    #     s=f.readline()
    #     while s!="":
    #         dictTemps[t].append(getVecteurFromLine(s))
    #         s=f.readline()
    #     f.close()
    # for t in dictTemps:
    #     print(f"{t} : {len(dictTemps[t])}" )
    #     # somme trueTotal
    # for indRun in range(len(dictTemps["Krimp"])):
    #     tempsTrue=[]
    #     for iterTemps in range(nbIter):
    #         tempsIter=0
    #         for indT in range(5):    
    #             t=tempsDetails[indT]
    #             tempsIter+=dictTemps[t][indRun][iterTemps]
    #         tempsTrue.append(tempsIter)
    #     dictTemps["trueTotal"].append(tempsTrue)
    # # calcul proportion true
    # proportionTempsIter=dict()            
    # for indT in range(5):
    #     t=tempsDetails[indT]
    #     proportionTempsIter[t]=[]
    # for iterTemps in range(nbIter):
    #     for indT in range(5):    
    #         t=tempsDetails[indT]
    #         tempsIterations=[]
    #         for indRun in range(len(dictTemps["Krimp"])):
    #             tempsIter=dictTemps[t][indRun][iterTemps]
    #             tempsIterations.append(tempsIter/dictTemps["trueTotal"][indRun][iterTemps])
    #         proportionTempsIter[t].append(np.mean(tempsIterations))

     ### get patternShownAll
    patternShownAll=[]
    patternShownDictSurp=dict()
    requetes=dict()
    jaccALLRun=[]
    for indRun in range(nbRun):
        requetes[indRun]=[]
        f=open(f"Res/{d}/Interacts/feedbacks{indRun}.txt")

        # Requete from list of feedback
        s=f.readline()
        while s!="":
            # print(s)
            if s!="[]\n":
                req=[]
                splitS=""
                for ch in s:
                    if (ch!='[') and ch!="]" and ch!="\\" and ch!="n":
                        splitS+=ch
                splitS=splitS.split("\n")[0].split(",")
            
                for ch in splitS:
                    tmpInt=int(ch)
                    if not(tmpInt in req):
                        req.append(tmpInt)
                
                reqPrevious=req
                # print(reqPrevious)
                requetes[indRun].append(req)  
            else:
                # print(reqPrevious)
                requetes[indRun].append(reqPrevious)

       
            s=f.readline()
        f.close()
        ###### Requete from list of pairs 
        # # # s=f.readline()
        # # # splitS=s.split(",")
        # # # first=int(splitS[0].split("[")[1])
        # # # second=int(splitS[1].split("]")[0])
        # # # req=[first,second]
        # # # s=f.readline()
        # # # while s !="":
            
        # # #     splitS=s.split(",")
        # # #     splitS=s.split(",")
        # # #     prevSecond=second
        # # #     first=int(splitS[0].split("[")[1])
        # # #     second=int(splitS[1].split("]")[0])
        # # #     if prevSecond == first or prevSecond ==second :
        # # #         req.append(second)
        # # #     else:
        # # #         requetes[indRun].append(req)
        # # #         req=[first,second] 
        # # #     s=f.readline()
        # # # requetes[indRun].append(req)
        # # # print(len( requetes[indRun]))
        # # # # print(requetes[0])
       ###########

        patternShownDictSurp[indRun]=dict()
        # print(f"indRun ={indRun}")
        patternShown=getPatternShown(d,len(dataI),len(data),indRun,sorted(list(dataI.keys())))


    ############## Surprise alll
        for pS in patternShown:
            # print(f" itemset = {pS}, nbItems ={len(pS)}, supp ={supp}")
            couv,supp=getSupp(dataI,pS)
           
            surprise=supp/len(data)
            surpriseItems=1
            for item in pS:
                surpriseItems*=freqItems[int(item)-1]
            patternShownDictSurp[indRun][len(patternShownDictSurp[indRun])]=[surprise-surpriseItems,0,0]
            if not(pS in patternShownAll):
                patternShownAll.append(pS)

    ######################################## JACCARD of shown pattern

        # jaccDict=np.zeros(101)
        # for iP in range(len(patternShown)):
        #     couvI,suppI=getSupp(dataI,patternShown[iP])
            
        #     for jP in range(iP+1,len(patternShown)):
        #         couvJ,suppJ=getSupp(dataI,patternShown[jP])
        #         jacc=get_Jaccard(couvI,couvJ)
        #         for jaccValue in range(100):
        #             if jaccValue/100>=jacc:
        #                 jaccDict[jaccValue]+=1
        #         jaccDict[100]+=1

        # for i in range(101):
        #     jaccDict[i]=jaccDict[i]/jaccDict[100]
        # jaccALLRun.append(jaccDict)


                


    # valueToPrintJacc=""
    # for i in [0,5,10,20,30,50]:
    #     valueTmp=0
    #     for j in jaccALLRun:
    #         valueTmp+=j[i]
    #     valueToPrintJacc+=f"{i}:{valueTmp/len(jaccALLRun)}"
    #     valueToPrintJacc+="; "
    # print(f"Jaccard pour {d} = {valueToPrintJacc}")
    
    print("Gettings Test set")



#     ############################################################################# rang percentile 

#      # get 145 k motifs and surprise  
    # LCM C
    os.system(f"./lcm53 C data/{d}.dat {theta} resLCM.txt >log.txt")
    os.system(f"wc -l resLCM.txt")
    print(" extract fini")
    f=open("resLCM.txt")
    # nbp=np.round(len(motifs)*(len(motifs)+1)/2)
    # upd=5
    s=f.readline()
    allItemsets=[]
    compteurMotif=0
    while s!="":
        compteurMotif=compteurMotif+1
        if compteurMotif % 35000 ==0:
            print(f"{compteurMotif} fait")
        i=1
        splitS=s.split()
        itemset=[]
        while splitS[i]!="]":
            itemset.append(splitS[i])
            i=i+1
        nbI=len(itemset)
        if nbI>1:
            couv,supp=getSupp(dataI,itemset)
            surprise=supp/len(data)
            surpriseItems=1
            for item in itemset:
                surpriseItems*=freqItems[int(item)-1]
            for indRunDict in patternShownDictSurp:
                for indPattern in patternShownDictSurp[indRunDict]:
                    if patternShownDictSurp[indRunDict][indPattern][0] > surprise-surpriseItems :
                        patternShownDictSurp[indRunDict][indPattern][1]+=1
                    else:
                        if patternShownDictSurp[indRunDict][indPattern][0] == surprise-surpriseItems :
                            patternShownDictSurp[indRunDict][indPattern][2]+=1
        s=f.readline()
    f.close()
   

   
    
    Y=[]
    Z=[]
    X=[]
    maxReq=0
    for indRunDict in patternShownDictSurp:
        y=[]
        z=[]
        x=[]
        # maxReq=max(maxReq,len(requetes[indRunDict]))
        for req in requetes[indRunDict]:
            rankPercentile=[]
            if requetes[indRunDict].index(req)!=0:
                newPrevious=motifFromPrevious
                for indPattern in req:
                    if indPattern!=motifFromPrevious:
                        rankPercentile.append((patternShownDictSurp[indRunDict][indPattern][1]+ (patternShownDictSurp[indRunDict][indPattern][2]/2))/compteurMotif)
                    if patternShownDictSurp[indRunDict][indPattern][0]>patternShownDictSurp[indRunDict][motifFromPrevious][0]:
                               newPrevious=indPattern
                motifFromPrevious=newPrevious
            else:
                motifFromPrevious=req[0]
                for indPattern in req:
                    rankPercentile.append((patternShownDictSurp[indRunDict][indPattern][1]+ (patternShownDictSurp[indRunDict][indPattern][2]/2))/compteurMotif)
                    if patternShownDictSurp[indRunDict][indPattern][0]>patternShownDictSurp[indRunDict][motifFromPrevious][0]:
                            motifFromPrevious=indPattern
            
            
            # if len(rankPercentile)>5:
            #     print(f" iter {indRunDict}  req {requetes[indRunDict].index(req)} rankPercentile = {rankPercentile}")
            # if rankPercentile==[]:
            #     print(f" iter {indRunDict}  req {requetes[indRunDict].index(req)} rankPercentile = {rankPercentile}")
            y.append(1-np.mean(rankPercentile))
            z.append(1-max(rankPercentile))
            x.append(1-min(rankPercentile))
        # print(f"x={x},\n y={y},\n z={z}")
        Y.append(y)
        Z.append(z)
        X.append(x)

        # plt.scatter(x,y,label=f'run {indRunDict} Avg') 
        # plt.scatter(x,z,label=f'run {indRunDict} Max') 
        # plt.plot(y,linewidth=0.5,color=couleur[indRunDict],linestyle="dashed",label=f'run {indRunDict} Avg') 
        # plt.plot(z,linewidth=0.5,color=couleur[indRunDict],label=f'run {indRunDict} Max') 

    # print(f"X={X},\n Y={Y},\n Z={Z}")
    meanY=[]
    maxZ=[]
    minZ=[]
    resRegret=open(f"Res/{d}/regret.txt","w")
    # X=[i for i in range(maxReq)]
    print(maxReq)
    for i in range(nbIter) :
        meantmp=[]
        maxtmp=[]
        mintmp=[]
        for j in range(nbRun):
            meantmp.append(np.round(Y[j][i],5))
            maxtmp.append(np.round(Z[j][i],5))
            mintmp.append(np.round(X[j][i],5))
        if meanY==[]:

            meanY.append(np.mean(meantmp))
            maxZ.append(np.mean(maxtmp))
            minZ.append(np.mean(mintmp))
        else:
            if meantmp == []:
                meanY.append(meanY[-1])
                maxZ.append(maxZ[-1])
                minZ.append(minZ[-1])
            else:
                meanY.append(meanY[-1]+np.mean(meantmp))
                maxZ.append(maxZ[-1]+np.mean(maxtmp))
                minZ.append(minZ[-1]+np.mean(mintmp))
    
    print(f"minZ={minZ},\nmeanY={meanY},\nmaxZ={maxZ}")
    plt.plot(meanY,linewidth=1.5,color=couleur[dataset.index(d)],linestyle="dashed",label=f'{d} Avg') 
    plt.plot(minZ,linewidth=1.5,color=couleur[dataset.index(d)],linestyle="dotted",label=f'{d} Min')
    plt.plot(maxZ,linewidth=1.5,color=couleur[dataset.index(d)],label=f'{d} Max') 
    resRegret.write(f"mean:{meanY}\n")
    resRegret.write(f"1-max:{maxZ}\n")
    resRegret.write(f"1-min:{minZ}")
    fileCumulRgt.write(f"{d};{np.round(minZ[-1],2)};{np.round(meanY[-1],2)};{np.round(maxZ[-1],2)}\n")
    resRegret.close()
    
    ###########################################################
   
    # Analyse pct rank

    # ######   Analyse distribution


     ### get Candidates or code tables
#     motifs=[]
#     allItemsets=[]
#     indCandidatsRun=[]
   
#     # patternShownAll=[]
#     for indRun in range(nbRun):
#         # print(f"indRun ={indRun}")
#         # patternShown=getPatternShown(d,len(dataI),len(data),indRun)
#         # for pS in patternShown:
#         #     if not(pS in patternShownAll):
#         #         patternShownAll.append(pS)
#         ### GET ensemble des candidats 
#         if os.path.isfile(f"Res/{d}/Candidates/{indRun}/run0.isc"):
#             for indIter in [0,3,5,8,10,13,15,17,19]:# range(nbIter): #
#                 # print(f"indIter ={indIter}")
#                 f=open(f"Res/{d}/C,andidates/{indRun}/run{indIter}.isc")
#                 s=f.readline()
#                 s=f.readline()
#                 s=f.readline()
#                 splitS=s.split(":")[1].split()
#                 while len(splitS)>2:
#                     itemset=[]
#                     for ind in range(len(splitS)-1):
#                         itemset.append(dicoK[int(splitS[ind])])
#                     itemset=sorted(itemset)
#                     if not(itemset in allItemsets): 
#                         surprise=int(splitS[-1].split('(')[1].split(')')[0])/len(data)
#                         estimation=1
#                         for i in itemset:
#                             estimation*=len(dataI[i])/len(data)
#                         surprise-=estimation
#                         if surprise <0:
#                             surprise=0
#                         couv,freq=getFreq(dataI,itemset)
#                         # itemsetFeatures=[]
#                         # for i in range(len(data)):
#                         #     if i in couv:
#                         #         itemsetFeatures.append(1)
#                         #     else:
#                         #         itemsetFeatures.append(0)    
#                         # for i in range(len(dataI)):
#                         #     if i in itemset:
#                         #         itemsetFeatures.append(1)
#                         #     else:
#                         #         itemsetFeatures.append(0)  
#                         # itemsetFeatures.append(freq)
#                         # itemsetFeatures.append(len(itemset))
#                         # motifs.append([surprise,itemset,itemsetFeatures])
#                         motifs.append([surprise,itemset,freq,len(itemset)])
#                         allItemsets.append(itemset)

#                     s=f.readline()
#                     if s!="":
#                         splitS=s.split(":")[1].split()
#                     else:
#                         splitS=[""]
#                 f.close()
#                 indCandidatsRun.append([len(allItemsets),indIter])
#         ### Get ensemble codetable 
#         # if os.path.isfile(f"Res/{d}/Codetables/{indRun}/run0.ct"):
#         #     for indIter in range(nbIter): # [0,3,5,8,10,12,15,17,19]:
#         #         # print(f"indIter ={indIter}")
#         #         f=open(f"Res/{d}/Codetables/{indRun}/run{indIter}.ct")
#         #         s=f.readline()
#         #         s=f.readline()
#         #         s=f.readline()
#         #         splitS=s.split("(")[0].split()
#         #         while len(splitS)>2:
#         #             itemset=[]
#         #             for ind in range(len(splitS)-1):
#         #                 itemset.append(dicoK[int(splitS[ind])])
#         #             itemset=sorted(itemset)
#         #             if not(itemset in allItemsets): 
#         #                 couv,freq=getFreq(dataI,itemset)
#         #                 surprise=freq/len(data)
#         #                 estimation=1
#         #                 for i in itemset:
#         #                     estimation*=len(dataI[i])/len(data)
#         #                 surprise-=estimation
                        
#         #                 # itemsetFeatures=[]
#         #                 # for i in range(len(data)):
#         #                 #     if i in couv:
#         #                 #         itemsetFeatures.append(1)
#         #                 #     else:
#         #                 #         itemsetFeatures.append(0)    
#         #                 # for i in range(len(dataI)):
#         #                 #     if i in itemset:
#         #                 #         itemsetFeatures.append(1)
#         #                 #     else:
#         #                 #         itemsetFeatures.append(0)  
#         #                 # itemsetFeatures.append(freq)
#         #                 # itemsetFeatures.append(len(itemset))
#         #                 # motifs.append([surprise,itemset,itemsetFeatures])
#         #                 motifs.append([surprise,itemset,freq,len(itemset)])
#         #                 allItemsets.append(itemset)

#         #             s=f.readline()
#         #             if s!="":
#         #                 splitS=s.split("(")[0].split()
#         #             else:
#         #                 splitS=[""]
#         #         f.close()
#         #         indCandidatsRun.append([len(allItemsets),indIter])

#     # print(f"        indice Candidats={indCandidatsRun}")
#     # print(len(motifs))
#     # print(f"        Motifs total ={len(allItemsets)}")
#     # print(f"        motifs Shown ={len(patternShownAll)}")




#     # distribFreq=[]
#     # distribLen=[]
#     distribSurpInv=[]
#     distribSurp=[]
#     for m in motifs:
#         distribSurp.append(m[0])
#         distribSurpInv.append(-m[0])
#         pS=0
#     #     distribFreq.append(m[-2])
#     #     distribLen.append(m[-1])
#     # counts, bin_edges = np.histogram(distribFreq, bins=15)  # tu peux choisir le nombre de bins
#     # np.set_printoptions(suppress=True) 
#     # print(f"         Bins (bords) : {bin_edges}")
#     # print("         Fréquences :", counts)
#     # counts, bin_edges = np.histogram(distribLen, bins=15)  # tu peux choisir le nombre de bins
#     # print(f"         Bins (bords) : {bin_edges}")
#     # print("         Longueur :", counts)
#     # counts, bin_edges = np.histogram(distribSurp, bins=20)  # tu peux choisir le nombre de bins
#     # print(f"        Bins (bords) : {bin_edges}")
#     # print("         Surprise :", counts)
#     listTrue=np.argsort(distribSurpInv)
#     # aTrouver=listTrue[:30]
#     # # print(f"        listTrue :{aTrouver}") 

#     ##### Analyse sur le classement du top 
 
    
#     # # indRunTrouver=[]
#     # # HasBeenShown=[]
#     # # for i in aTrouver:
#     # #     j=0
#     # #     if allItemsets[i] in patternShownAll:
#     # #         HasBeenShown.append(1)
#     # #     else:
#     # #         HasBeenShown.append(0)
#     # #     while j<len(indCandidatsRun):
#     # #         if i > indCandidatsRun[j][0]:
#     # #             j+=1
#     # #         else:
#     # #             indRunTrouver.append(indCandidatsRun[j][1])
#     # #             j=len(indCandidatsRun)+1
#     # cptTop30=0
#     # for i in listTrue:
#     #     if allItemsets[i] in patternShownAll:
#     #         j=0
#     #         while j<len(indCandidatsRun):
#     #             if i > indCandidatsRun[j][0]:
#     #                 j+=1
#     #             else:
#     #                 print(f"    {motifs[i][1]} avec surp={motifs[i][0]} a iter {indCandidatsRun[j][1]}")
#     #                 j=len(indCandidatsRun)+1
#     #                 cptTop30+=1
#     #     if cptTop30>29:
#     #         break
#     # print(f"        Surprise:{[ motifs[i][0] for i in aTrouver]}")
#     # print(f"        IndRunTrouver:{indRunTrouver}")
#     # print(f"        HasBeenShown:{HasBeenShown}")
#     # dictSurpriseByIter=dict()
#     # for i in range(nbIter):
#     #     dictSurpriseByIter[i]=[]
#     # for i in listTrue:
#     #     j=0
#     #     while j<len(indCandidatsRun):
#     #         if i > indCandidatsRun[j][0]:
#     #             j+=1
#     #         else:
#     #             dictSurpriseByIter[indCandidatsRun[j][1]].append(motifs[i][0])
#     #             j=len(indCandidatsRun)+1

#     # for i in dictSurpriseByIter:
#     #     print(f"        mean surprise at iter {i} : {np.mean(dictSurpriseByIter[i])}")
#     ##### Analyse sur le rank globale
#     # Get learn weight
#     print(" Analyse rank globale")
#     weightLearned=[]
#     valuesKendall=[]
#     for indRun in range(nbRun):
#         weightL=[]
#         # f=open(f"Res/{d}/Interacts/weightLearned{indRun}.txt")
#         # f=open(f"Res/dispale/weigth{d}{indRun}.txt")  
#         f=open(f"Res/letsip/weigth{d}{indRun}.txt")       
#         s=f.readline()
#         cptLinge=0

#         while s!="":
            
#             ###### LCS
#             # splitS=s.split('],')
#             # bias=float(splitS[1])
#             # weight=getVecteurFromLine(splitS[0]+"]")
#             ###### Letsip dispale
#             weight=getVecteurFromLine(s)
#             weight.pop(0)
#             weight.pop(len(dataI.keys()))
#             bias=0
#             #####
#             if cptLinge==19:
#                 weightL.append([weight,bias])
#                 scoreML=[]
#                 scoreMLbis=[]
#                 for m in motifs:
#                     itemsetFeatures=[]
#                     itemset=m[1]
#                     couv,freq=getFreq(dataI,itemset)

#                     ## ####LCS
#                     for i in range(len(data)):
#                         if i in couv:
#                             itemsetFeatures.append(1)
#                         else:
#                             itemsetFeatures.append(0)    
                    
#                     for i in sorted(list(dataI.keys())):
#                         if i in itemset:
#                             itemsetFeatures.append(1)
#                         else:
#                             itemsetFeatures.append(0) 
#                     # itemsetFeatures.append(freq)
#                     # itemsetFeatures.append(len(itemset))
#                     ##############Letsip
#                     # for i in sorted(list(dataI.keys())):
#                     #     if i in itemset:
#                     #         itemsetFeatures.append(1)
#                     #     else:
#                     #         itemsetFeatures.append(0) 
#                     # for i in range(len(data)):
#                     #     if i in couv:
#                     #         itemsetFeatures.append(1)
#                     #     else:
#                     #         itemsetFeatures.append(0)    
#                     # itemsetFeatures.append(freq)
#                     # itemsetFeatures.append(len(itemset))
                    
#                     ##################""
#                     x=-(np.array(itemsetFeatures) @ np.array(weight)+bias)
#                     # x=-(np.array(itemsetFeatures[:len(data)]) @ np.array(weight[:len(data)]))
                    
#                     # x-=(np.array(itemsetFeatures[:-len(dataI.keys())]) @ np.array(weight[:-len(dataI.keys())]))
#                     # xbis=-(np.array(itemsetFeatures[:len(data)]) @ np.array(weight[:len(data)]))-(np.array(itemsetFeatures[:-len(dataI.keys())]) @ np.array(weight[:-len(dataI.keys())]))
#                     scoreML.append(x)
#                     # scoreMLbis.append(xbis)
#                 listML=np.argsort(scoreML)
#                 # listMLbis=np.argsort(scoreMLbis)
#                 # print(f"Nombre de motifs :{len(listML)}")
#                 # print(f"listTrue :{listTrue[:10]}")
#                 # print(f"listML:{listML[:10]}")
#                 # print("    recall@30 ML: ",len(set(listML[:30])& set(listTrue[:30]))/30)
#                 # print(f"listMLbis:{listMLbis[:10]}")
#                 # print("    recall@30 MLbis: ",len(set(listMLbis[:30])& set(listTrue[:30]))/30)
#                 # Testing with paires by sublist
#                 nbFauxML=0
#                 nbVraiML=0
#                 listTrue=list(listTrue)
#                 listML=list(listML)
#                 listDeajVu=set()
#                 x=[]
#                 y=[]
#                 for i in range(len(listTrue)):
#                     # if (i+1)%(len(listTrue)/5)==0:
#                         # print(" un cinquièeme de fait")
#                     motifCheck=listTrue[i]
#                     motifsPlusGrand=set(listTrue[0:i])
#                     motifsPlusPetit=set(listTrue[i+1:len(listTrue)])

#                     indML=listML.index(motifCheck)
#                     MPGML=set(listML[0:indML])
#                     MPPML=set(listML[indML+1:len(listML)])
#                     vraiML=(motifsPlusGrand & MPGML).union(motifsPlusPetit & MPPML)-listDeajVu
#                     fauxML=(motifsPlusGrand ^ MPGML).union(motifsPlusPetit ^ MPPML)-listDeajVu

#                     nbVraiML+=len(vraiML)
#                     nbFauxML+=len(fauxML)
#                     listDeajVu.add(motifCheck)
#                     y.append(i)
#                     x.append(listML.index(listTrue[i]))
                
#                 # print("     nbVraiML: ",nbVraiML)
#                 # print("     nbFauxML: ",nbFauxML)
#                 # print(f"    Kendal: {(nbVraiML-nbFauxML)/(nbVraiML+nbFauxML)}")
#                 valuesKendall.append((nbVraiML-nbFauxML)/(nbVraiML+nbFauxML))
#             # plt.figure(figsize=(20, 20))
#             # plt.xlabel("Learned rank")
#             # plt.ylabel("True rank")

#             # plt.title(f"{d} run{indRun}-iter{cptLinge} ")
#             # plt.tight_layout()
#             # plt.scatter(x,y,linewidth=0,s=2)
#             # plt.scatter(y,y,linewidth=0,s=2,c='black')
#             # plt.savefig(f"Res/{d}/Interacts/run{indRun}-iter{cptLinge}.eps", format='eps', dpi=300, bbox_inches='tight')
#             # plt.clf()
#             s=f.readline()
#             cptLinge+=1
#         weightLearned.append(weightL)
#         f.close()

#     print(f"Pour letsip et {d}")
#     print(f"valuesKendall:{valuesKendall}")
#     print(f"mean valuesKendall:{np.mean(valuesKendall)}")
plt.plot([1+tmp for tmp in range(20)],linewidth=1.5,color='black',label=f'Maximum') 
plt.legend()
plt.show()
# plt.savefig(f"Res/{d}/Interacts/run{indRun}-iter{cptLinge}.eps", format='eps', dpi=300, bbox_inches='tight')
plt.clf()  
fileCumulRgt.close()