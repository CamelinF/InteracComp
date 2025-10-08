import os,sys,math,time,random,bisect
from skmine.datasets.fimi import fetch_file
import igraph as ig
import matplotlib.pyplot as plt

from graph import *

# get cover and frequency of pattern m from dataset d 
def getFreq(d,m):
    couv=d[m[0]]
    for i in m:
        couv=couv.intersection(d[i])
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


def dataTransToDataItem(data):
    dataI={}
    for t in range(len(data)):
        for i in data[t]:
            if i in dataI:
                dataI[i].add(t)
            else:
                dataI[i]={t}
    return dataI
# getting codetable CT from fileNamCT for dataset fileNamD
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
        motif.sort()
        tmp=splitS[-1].split(",")[0]
        tmpp=int(tmp.split("(")[1])
        tmp=splitS[-1].split(",")[1]
        freq=int(tmp.split(")")[0])
        if len(motif)>1:
            CTK.append([freq,len(itemset),itemset])
            CTD.append([freq,len(itemset),motif])
        s=f.readline()
        splitS=s.split()
    f.close()
    return CTK,CTD


def getDicoCouv(dataI,listMotifs):
    dictCouv=dict()
    cpt=0
    for i in listMotifs:
        couv,freqM=getFreq(dataI,i[2])
        dictCouv[cpt]=(i[2],couv)
        cpt+=1
        # print(f"frequence K {i[0]} vs true frequence {freqM}")
    nbMotifs=len(listMotifs)
    print(nbMotifs," vs " , cpt )
    return dictCouv

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


if __name__ == "__main__":
    timeDebut=time.time()   

    fileNameData=str(sys.argv[1])
    fileNameCT=str(sys.argv[2])
    seuilJacc=float(sys.argv[3])
    support=int(sys.argv[4])
    data=fetch_file(f'datasets/{fileNameData}.dat',int_values=True)
    dataI=dataTransToDataItem(data)
    CTK,CTD=getCT(fileNameData,fileNameCT)
    # CTk,CTD=getCand(fileNameData,fileNameCT)
    dictCouv=getDicoCouv(dataI,CTD)
   
    linksBetweenMotifs=getLinkFromDico(dictCouv,seuilJacc)
    #linksBetweenMotifs=getLinkFromPattern(CTD,support)
    # linksBetweenMotifs,weights=getAllLinkFromDico(dictCouv)
    # li nksBetweenMotifs,weights=getIntersecLinkFromPattern(CTD)
    g = ig.Graph(n=len(CTD), edges=linksBetweenMotifs, directed=False)
    # Enlève les noeuds isolés
    # non_isolated_vertices = [v.index for v in g.vs if g.degree(v.index) > 0]
    # print(f" on affiche seulements {len(non_isolated_vertices)} sommets")
    # g = g.subgraph(non_isolated_vertices)
    # g.es["weight"]=weights
    


    # clusters = g.community_multilevel()
    # g.vs["cluster"] = clusters.membership  # Assigner l'indice de cluster à chaque sommet
    # print(f"Nombre de communautés : {len(clusters)}")


    pageranks = g.pagerank(directed=True)
    print(f"on  a mis {time.time() - timeDebut} secondes")
   
    # Affichage des scores de PageRank
    top_indices = np.argsort(pageranks)[-15:][::-1]

    # Affichage des 10 meilleurs
    print("Top 15 des sommets par PageRank :")
    for rank, idx in enumerate(top_indices, start=1):
        print(f" Sommet {idx} {dictCouv[idx][0]}  — PageRank: {pageranks[idx]:.5f}")


    fig,ax=makeGraph(g,pageranks)
    
    
    plt.show()