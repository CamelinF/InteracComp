import os,sys,math,time,random
from skmine.datasets.fimi import fetch_file


import afficheCT as affCT




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

# getting codetable CT from fileNamCT for dataset fileNamD
def getCT(fileNamD,fileNamCT):
    dicoK= getDecodage(fileNamD,True)
    f=open(f"codetables/{fileNamCT}.ct")
    s=f.readline()
    s=f.readline()
    s=f.readline()
    splitS=s.split()
    CT=[]
    while s!="":
        itemset=[]
        motif=[]
        for indIt in range(len(splitS)-1):
            item=splitS[indIt]
            itemset.append(int(item))
            motif.append(dicoK[int(item)])
        itemset.sort()
        tmp=splitS[-1].split(",")[0]
        tmpp=int(tmp.split("(")[1])
        tmp=splitS[-1].split(",")[1]
        freq=int(tmp.split(")")[0])
        #couv,freq=getFreq(dataI,motif)
        if len(motif)>1:
            CT.append([freq,len(itemset),itemset])
        s=f.readline()
        splitS=s.split()
    f.close()
    return CT
if __name__ == "__main__":

    fileNameData=str(sys.argv[1])
    fileNameCT=str(sys.argv[2])
    data=fetch_file(f'datasets/{fileNameData}.dat',int_values=True)
    CT=affCT.getCT(fileNameData,fileNameCT)
    for i in CT:
        print(i)