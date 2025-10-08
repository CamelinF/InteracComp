
import os
import sys
import time 
import math
import random
import skmine
import time
import bisect
import numpy as np 

from skmine.datasets.utils import describe
from skmine.datasets.fimi import fetch_file
from skmine.itemsets import LCM
from random import sample,random,choices,randint
import pulp
from pulp import GUROBI

import networkx as nx
import matplotlib.pyplot as plt

dataTypeKrimp = {
    "connect" : "bai32",
    "hepatitis" :"bm128",
    "mushroom" : "bm128",
    "splice1" : "bai32",
    "pumsb" : "uint16",
    "pumsb_star" : "uint16",
    "eisen" : "uint16",
    "retail" :"uint16",
    "chess" :"bm128",
}

def getFreq(d,m):
    couv=d[m[0]]
    for i in m:
        couv=couv.intersection(d[i])
        #print(f"{m} à uen couv de {len(couv)}")
    return couv,len(couv)

def isRectangular(data):
    ref=len(data[0])
    rectangular=True
    for i in data:
        if len(i)!=ref:
            rectangular=False
    return rectangular


def motifDecode(m,decode):
    motifD=[]
    for i in m:
        motifD.append(decode[i])
    return motifD

def motifStrToInt(m):
    motifInt=[]
    tmp=m.split(',')
    for i in tmp:
        motifInt.append(int(i))
    return motifInt
def motifIntToStr(m):
    motifStr=""
    for i in range(len(m)-2):
        motifStr+=f"{m[i]},"
    motifStr+=f"{m[-2]}"
    return motifStr

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



def extract(dataFile,dataT,dataI,ensTrans,ensMotif,typeM,minSupp):
    lcm=LCM(min_supp=minSupp)
    #print("extraction des motifs")
    # temps=time.time()
    newCandidats=lcm.fit_transform(ensTrans)
    dicoK= getDecodage(dataFile,False)
    nbS=0
    longEns=len(ensMotif)
    for i in newCandidats.index:
        itemset=newCandidats.iat[i,0]
        motifK=[]
        for i in itemset:
            motifK.append(dicoK[int(i)])
        motifK.sort()
        nbI=len(itemset)
        if nbI>1:
            present=False 
            for i in range(len(ensMotif)):
                if motifK ==ensMotif[i][2]:
                    present=True
                    break
            if not(present):
                couv,freq=getFreq(dataI,itemset)
                #print(f" on a fait le motif {itemset} avec un freq de {freq}")
                ensMotif.append([freq/len(dataT),nbI,motifK,couv])
        else:
            nbS=nbS+1
    r=open(f"Res/{dataFile}/trouve.txt","a")
    r.write(str(len(newCandidats)-nbS))
    r.write(",")
    r.close()
    #print(" on a extraits "+str(len(newCandidats)-nbS))
    r=open(f"Res/{dataFile}/ajout.txt","a")
    r.write(str(len(ensMotif)-longEns))
    r.write(",")
    r.close()
    #print(f" Il y a {len(ensMotif)-longEns} nouveaux motifs")
    ensMotif.sort(key=lambda x: (-x[0],-x[1],x[2]))
    #print(" on encode pour Krimp")
    r=open(f"data/candidates/{dataFile}-{typeM}-"+str(minSupp)+"d.isc","w")
    r.write("ficfis-1.3\n")
    r.write("mi: numSets="+str(len(ensMotif)))
    maxL=[1]
    for i in ensMotif:
        maxL.append(i[1])
    r.write(" minSup="+str(minSupp)+" maxLen="+str(max(maxL)))
    r.write(f" sepRows=0 iscOrder=d patType={typeM} dbName={dataFile}\n")
    for i in range(len(ensMotif)):
        freqM=len(ensMotif[i][3])
        longM=ensMotif[i][1]
        motif=ensMotif[i][2]
        r.write(str(longM)+":")
        for j in range(longM):
            r.write(" "+str(motif[j]))
        r.write(" (" + str(freqM)+")\n")
    r.close()
    #print(f" ça a pris {time.time()-temps} secondes")

    return ensMotif





def extractCT(dataFile,dataI,indRun,ensMotifA,dim):
    ensMotif=[]
    ensCT=[]
    ensCTFC=[]
    sumTot=0
    poids=[]
    dicoK= getDecodage(dataFile,True)
    os.system(f"rm experiments/{dataFile}{indRun}/compress/{dataFile}*/*-0.ct")
    os.system(f"mv experiments/{dataFile}{indRun}/compress/{dataFile}*/*.ct resKrimp.ct")
    f=open(f"resKrimp.ct")
    s=f.readline()
    s=f.readline()
    s=f.readline()
    splitS=s.split()
    while s!="":
        itemset=[]
        motif=[]
        for indIt in range(len(splitS)-1):
            item=splitS[indIt]
            itemset.append(int(item))
            motif.append(dicoK[int(item)])
        itemset.sort()
        tmp=splitS[-1].split(",")[0]
        nbUsage=int(tmp.split("(")[1])
        tmp=splitS[-1].split(",")[1]
        freq=int(tmp.split(")")[0])
        if len(motif)>1:
            index=0
            while (ensMotifA[index][2]!=itemset) and (index<len(ensMotifA)):
                index+=1
            true_score=userFunction(dataI,[freq,motif],dim)
            
            motifFeatures=[]
            for i in range(dim[0]):
                if i in ensMotifA[index][3]:
                    motifFeatures.append(1)
                else:
                    motifFeatures.append(0)    
            for i in range(dim[1]):
                if i in motif:
                    motifFeatures.append(1)
                else:
                    motifFeatures.append(0)    
            ensMotif.append([freq,len(itemset),itemset,ensMotifA[index][3],true_score,motifFeatures])
        sumTot+=nbUsage
        poids.append(nbUsage)
        ensCT.append([motif,freq,nbUsage])
        s=f.readline()
        splitS=s.split()
    
    tailleComp=[]
    for i in range(len(poids)):
        x=poids[i]/sumTot
        tailleComp.append(x)
    os.system(f"mv resKrimp.ct Res/{dataFile}/Codetables/run{indRun}.ct")
    

    return ensMotif,ensCT,tailleComp

def newSampler(data,poidsD,c):
    ensTrans=[]
    #print(sum(poidsD))
    ensIndex=np.random.choice(len(data),size=c,replace=False,p=poidsD)
    
    # check de np.random.choice
    # for i in range(c):
    #     for j in range(i+1,c):
    #         if ensIndex[i]==ensIndex[j]:
    #             print(" c'est la merde")
    #             print(ensIndex)
    for i in ensIndex:
        ensTrans.append(data[i])
    return ensTrans,ensIndex


def calculPoidsD(data,tailleI,ensCT,tailleD):
    poidsD=[]
    comp=[]

    for i in range(len(data)):

        d=set(data[i])
        tailleC=0
        j=0
        while j <len(ensCT):
            if set(ensCT[j][0]) in d:
                tailleC+=-1*math.log2(tailleD[j])
                for k in ensCT[j][0]:
                    d.remove(k)
            if len(ensCT[j][0])<2  and len(d)==len(data[i]):
                j=len(ensCT)       
                tailleC=tailleI[i]
            j+=1
        
        p=round(tailleC)/round(tailleI[i])
        if p >1:
            comp.append(1)#len(data[i]))
        else:
            comp.append(round(p,3))#*len(data[i]))
        
    somC=sum(comp)
    for i in comp:
        poidsD.append((i/somC))

    return poidsD


def userFunction(dataI,motif,dim):
    # score = freqMotif - Prod(freqItems) for item in motif
    score=np.round(motif[0]/dim[0],5)
    tmp=1
    for i in motif[1]:
        tmp=tmp*np.round((len(dataI[i])/dim[0]),5)
    score-=np.round(tmp,5)
    return score

def get_Feedbacks(request,ensMotif):
    feedbacks=[]
    scores=[ensMotif[request[i]][0][4] for i in range(len(request))]
    decreaseOrder= np.argsort(scores)[::-1]
    # print(decreaseOrder)
    # print(scores)
    # i=0
    # if scores[decreaseOrder[0]] > scores[decreaseOrder[-1]]+ 0.00001:
    #     while i<len(scores)-1:
    #         tmp=[request[decreaseOrder[i]]]
    #         while i<len(scores)-1 and scores[decreaseOrder[i]] == scores[decreaseOrder[i+1]]:
    #             tmp.append(request[decreaseOrder[i+1]])
    #             i+=1
    #         for j in tmp:
    #             if i ==len(request)-1:
    #                 if tmp.index(j)!=len(request)-1:
    #                     feedbacks.append([j,request[decreaseOrder[i]]])
    #             else:
    #                 feedbacks.append([j,request[decreaseOrder[i+1]]])
    #         i+=1
    for i in range(len(request)-1):
        for j in range(i+1,len(request)):
            scoreI=ensMotif[request[i]][0][4]
            scoreJ=ensMotif[request[j]][0][4]
            if scoreI> scoreJ:
                feedbacks.append([request[i],request[j]])
            else: 
                if scoreI!=scoreJ:
                    feedbacks.append([request[j],request[i]])
    return feedbacks

def learn_From_Feedbacks(feedbacks,ensMotif,n_features):
    model = pulp.LpProblem("RankLearn", pulp.LpMinimize)
    w = [pulp.LpVariable(f"w_{k}", lowBound=-1, upBound=1) for k in range(n_features)]
    violations = []
    epsilon = 1e-3
    lambda_l1=0.1

    # L1 constraints : u_k ≥ |w_k|
    # u = [pulp.LpVariable(f"u_{k}", lowBound=0, cat="Continuous") for k in range(n_features)]
    # for k in range(n_features):
    #     model += u[k] >= w[k]
    #     model += u[k] >= -w[k]
    for i, j in feedbacks:
        vij = pulp.LpVariable(f"v_{i}_{j}", cat="Binary")
        diff = np.array(ensMotif[i][1]) - np.array(ensMotif[j][1])
        model += pulp.lpSum([w[k] * diff[k] for k in range(n_features)]) >=  epsilon - n_features*2 * vij
        violations.append(vij)

    model += pulp.lpSum(violations)#+lambda_l1*pulp.lpSum(u)
    # model += pulp.lpSum(violations)

    optionsG = [
    #("TimeLimit" , 60.0),
    ("WLSACCESSID", "146955fd-fb4c-409b-bcd0-7eb0472c8d4b"),
    ("WLSSECRET", "ae52ec3b-abb7-4028-9d66-285e40efdc7d"),
    ("LICENSEID", 2610563)
    ]

    solver = GUROBI(msg=False)
    status = model.solve(solver)

    w_pulp = np.array([pulp.value(wk) for wk in w])
    w_pulp= np.array([np.random.uniform(-1, 1) if x is None else x for x in w_pulp], dtype=float)

    return np.round(w_pulp, 3)


def get_Request(ensMotif,poidsT,weights,k,prevRequest,alreadyShown):
    scores=[]
    ensMotifFeatures=[]

    for i in ensMotif:
        motifFeatures=i[5] # vecteur_features transac et items
        motifFeatures.append(i[0])#freq
        motifFeatures.append(i[1])#len
        # motifFeatures.append(i[4])#compression
        scores.append(np.dot(motifFeatures,weights))
        ensMotifFeatures.append(motifFeatures)

    # get best from preivous feedback :
    if len(prevRequest)!=0:
        prevMotif=[]
        for i in prevRequest:
            prevMotif.append(alreadyShown[i][0][4])
        indReq=np.argsort(prevMotif)[-2:][::-1]
        request=[]
        for i in indReq:
            request.append(prevRequest[i])
        i=0
        tmp=np.argsort(scores)[::-1]
        while len(request)<k and i <len(tmp):
            dedans=True
            for ind in alreadyShown:
                if set(ensMotif[tmp[i]][2]) == set(alreadyShown[ind][0][2]):
                    dedans=False
            if dedans:
                newInd=len(alreadyShown)
                alreadyShown[newInd]=[ensMotif[tmp[i]],ensMotifFeatures[tmp[i]]]
                request.append(newInd)
            i+=1
    
    else:  
        request=np.argsort(scores)[-k:][::-1] # retourne les indices des k plus grands learned score
        for i in range(len(request)):
            newInd=len(alreadyShown)
            alreadyShown[newInd]=[ensMotif[request[i]],ensMotifFeatures[request[i]]]
            request[i]=newInd

    return request,alreadyShown


def show_GraphFeedback(dataFile,prefs,indRun):
    maxInd=0
    for i,j in prefs:
        maxInd=max(maxInd,i)
        maxInd=max(maxInd,j)
    
    # Création du graphe dirigé
    G = nx.DiGraph()
    G.add_nodes_from(range(maxInd+1))  
    G.add_edges_from(prefs)

    # Affichage
    plt.figure(figsize=(15, 10))
    pos = nx.spring_layout(G, seed=42)  # Positionnement stable
    nx.draw_networkx_nodes(G, pos, node_size=500, node_color='lightblue')
    nx.draw_networkx_edges(G, pos, arrowstyle='->', arrowsize=20, edge_color='gray')
    nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold')

    plt.title(f"Graphe des préférences ({maxInd} motifs, {len(prefs)} arcs)")
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(f"Res/{dataFile}/GraphPrefs/run{indRun}.eps", format='eps', dpi=300, bbox_inches='tight')
    plt.clf()
    plt.close()

def main(dataFile,nbR,nbT,typeM,seuilF,dirCurr,k):
    
    tempsTotal=time.time()

    # Initialization data and weights
    data=fetch_file(f'data/{dataFile}.dat',int_values=True)
    poidsT=[1/len(data)]*len(data)
    
    minSupp=int(math.ceil(nbT*seuilF))

    ### infoD dict_keys(['n_items', 'avg_transaction_size', 'n_transactions', 'density'])
    #infoD=describe(data)
    #infoD["isRect"]=isRectangular(data)

    # Transposing the dataset for further manipulation
    dataI={}
    for t in range(len(data)):
        for i in data[t]:
            if i in dataI:
                dataI[i].add(t)
            else:
                dataI[i]={t}
    # Gettings initial lenght of each transaction
    sumInit=0
    for i in dataI:
        sumInit+=len(dataI[i])
    tailleI=[]
    for i in data:
        tailleTmp=0
        for j in i:
            tailleTmp+=-1*math.log2(len(dataI[j])/sumInit)
        tailleI.append(tailleTmp)
    
    # POids appris ( un par transaction + un par item + 1 pour la fréquence + +1 pour la taille + 1 pour la compression)
    weightLearned=[1]*(len(data)+len(dataI)+1+1)
    feedbacks=[]
    # Time initialization
    tempsTotalF=[] # initialize at the start to consider previous code into account
    # tempsTotal=time.time()
    tempsSampleF=[]
    tempsSample=0
    tempsKrimpF=[]
    tempsKrimp=0
    tempsWeightF=[]
    tempsWeighT=0
    tempsReqF=[]
    tempsReq=0
    tempsLearnF=[]
    tempsLearn=0
    #print(f" on mine dans {dataFile} les motifs {typeM} avec un seuil de {minSupp}")

    i=0
    ensMotif=[]
    ensRemove=[]
    tmpTaille=[]
    tmpTemps=[]
    request=[]
    alreadyShown=dict()
    while i<nbR:
        # two-step pattern sampling 
        tempsSample=time.time()
        sampleT,indexTransaction=newSampler(data,poidsT,nbT)
        ensMotif=extract(dataFile,data,dataI,sampleT,ensMotif,typeM,minSupp)
        while ensMotif==[]:
            sampleT,indexTransaction=newSampler(data,poidsT,nbT)
            ensMotif=extract(dataFile,data,dataI,sampleT,ensMotif,typeM,minSupp)

        tempsSampleF.append(time.time()-tempsSample)

        print(f"   Sample {i+1} de motifs faits")
        ### Récuparation des transactions sélectionnés
        # f=open(f"Res/{dataFile}/Candidates/runT{i}.txt","a")
        # f.write(f"{indexTransaction},{sampleT}")
        # f.write("\n")
        # f.close()
        
        # Calling Krimp each time # Calling Krimp every x samples
        if i>-1:#7 and i%9==0:

            # Krimp laucnh
            tempsKrimp=time.time()
            f=open("Krimp/bin/datadir.conf","w")
            f.write(f"dataDir = {dirCurr}/data/\n")
            f.write(f"expDir = {dirCurr}/experiments/{dataFile}{i}/\n")
            f.close()
            os.system("Krimp/bin/krimp > result.txt")
            # print(" Krimp Fait")
            # Krimp extraction
            #ensMotif,ensCT,couvKrimp,tailleComp,ensCTFC=extractCT(dataFile,dataI,i)

            # Here ensMotif now has 6 entries freq,len,itemset,cover,user_score, vector_Features
            ensMotif,ensCT,tailleComp=extractCT(dataFile,dataI,i,ensMotif,[len(data),len(dataI)])
            tempsKrimpF.append(time.time() - tempsKrimp)

            if i<nbR-1:
                # Calculating new compression weight
                tempsWeight=time.time()
                poidsT=calculPoidsD(data,tailleI,ensCT,tailleComp)
                tempsWeightF.append(time.time()-tempsWeight)

              
            else:
                tempsWeightF.append(0)
            # User feedback and learning from it
            tempsReq=time.time()

            # return k motifs de ensCT
            request,alreadyShown=get_Request(ensMotif,poidsT,weightLearned,k,request,alreadyShown)
            # print(f"request:{request}")
            
            newFeedbacks=get_Feedbacks(request,alreadyShown)
            
            
            # print(f"newFeedbacks:{newFeedbacks}")
            f=open(f"Res/{dataFile}/feedbacks.txt","a") 
            for j in newFeedbacks:
                # print(j)
                feedbacks.append(j)
                f.write(f"{j}\n")
            f.close() 
            # print(f"feedbacks:{feedbacks}")
            tempsReqF.append(time.time()-tempsReq)
            tempsLearn=time.time()
            weightLearned=learn_From_Feedbacks(feedbacks,alreadyShown,len(weightLearned))
            # print(f"weightLearned:{weightLearned}")
            for j in range(len(poidsT)):
                poidsT[j]=poidsT[j]*(1.1+weightLearned[j])
            sommeP=sum(poidsT)
            for j in range(len(poidsT)):
                poidsT[j]=poidsT[j]/sommeP

            # print(f"PoidsT:{poidsT}")
            tempsLearnF.append(time.time()-tempsLearn)
            show_GraphFeedback(dataFile,feedbacks,i)
        #print(ensMotif)
            f=open("result.txt")
            stop=True
            while stop :
                s=f.readline()
                splitS=s.split()
                if len(splitS)>1:
                    if splitS[1]=='Time:':
                        stop=False
            tmpTemps.append(float(splitS[6]))
            stop=True
            while stop:
                s=f.readline()
                splitS=s.split()
                if len(splitS)>1:
                    if splitS[1]=='Result:':
                        stop=False
            tmpTaille.append(int(splitS[2].split(',')[5].split(')')[0]))
            f.close()
        os.system(f"mv data/candidates/{dataFile}-{typeM}-{minSupp}d.isc Res/{dataFile}/Candidates/run{i}.isc")
        i=i+1
        tempsTotalF.append(time.time()-tempsTotal)
        tempsTotal=time.time()
    
    f=open(f"Res/{dataFile}/tempsSample.txt","a")
    f.write(str(tempsSampleF))
    f.write("\n")
    f.close()
    f=open(f"Res/{dataFile}/tempsTotal.txt","a")
    f.write(str(tempsTotalF))
    f.write("\n")
    f.close()
    f=open(f"Res/{dataFile}/tempsWeight.txt","a")
    f.write(str(tempsWeightF))
    f.write("\n")
    f.close()
    f=open(f"Res/{dataFile}/tempsKrimp.txt","a")
    f.write(str(tempsKrimpF))
    f.write("\n")
    f.close()
    f=open(f"Res/{dataFile}/tempsLearn.txt","a")
    f.write(str(tempsLearnF))
    f.write("\n")
    f.close()
    f=open(f"Res/{dataFile}/tempsReq.txt","a")
    f.write(str(tempsReqF))
    f.write("\n")
    f.close()
    f=open(f"Res/{dataFile}/taille.txt","a")
    f.write(str(tmpTaille))
    f.write("\n")
    f.close()
    f=open(f"Res/{dataFile}/weightLearned.txt","a")
    f.write(f"[{weightLearned[1]}")
    for i in range(1,len(weightLearned)):
        f.write(f",{weightLearned[i]}")
    f.write("]")
    f.write("\n")
    f.close() 
    print(len(alreadyShown[1]))
    f=open(f"Res/{dataFile}/patternShown.txt",'a')
    for i in range(len(alreadyShown)):
       f.write(f"[{i},{alreadyShown[i][1]}]\n")
    f.close()
            
if __name__ == "__main__":
    #np.random.seed(4111998)
    # Récupération des paramètres
    data=str(sys.argv[1])
    tailleS=int(sys.argv[2])
    nbRun=int(sys.argv[3])
    seuilFreq=float(sys.argv[4]) 
    
    k=5
    # Initialisation Krimp
    f=open("compress.conf","r")
    r=open("compress.txt","w")
    s=f.readline()
    while(s!=""):
        if s[0]=="i" and s[1] =="s" and s[2]=="c" and s[3]=="N":
            r.write(f"iscName={data}-closed-"+str(math.ceil(tailleS*seuilFreq))+"d\n")
        elif s[0]=="d" and s[1] =="a" and s[2]=="t" and s[3]=="a" and s[3]=="T" and s[3]=="y":
            r.write(f"dataType= {dataTypeKrimp[data]}")
        else:
            r.write(s)
        s=f.readline()
    r.close()
    f.close()
    os.system("rm compress.conf")
    f=open("compress.txt","r")
    r=open("compress.conf","w")
    s=f.readline()
    while(s!=""):
       r.write(s)
       s=f.readline()
    r.close()
    f.close()
    os.system("rm compress.txt")
    os.system("rm -r experiments/*")

    # Création des dossiers de résultats
    os.system(f"mkdir Res/{data}")
    os.system(f"mkdir Res/{data}/Candidates")
    os.system(f"mkdir Res/{data}/Codetables")
    os.system(f"mkdir Res/{data}/GraphPrefs")
    f=open(f"Res/{data}/ajout.txt","w")
    #f.write(f"{tailleS}-{nbRun}-{seuilFreq}:")
    f.close()
    f=open(f"Res/{data}/trouve.txt","w")
    #f.write(f"{tailleS}-{nbRun}-{seuilFreq}:")
    f.close()


    currDir=os.getcwd()
    main(data,nbRun,tailleS,"closed",seuilFreq,currDir,k)