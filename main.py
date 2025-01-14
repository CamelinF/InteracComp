import os,sys,math,time,random,bisect
from skmine.datasets.fimi import fetch_file

# get cover and frequency of pattern m from dataset d 
def getFreq(d,m):
    couv=d[m[0]]
    for i in m:
        couv=couv.intersection(d[i])
        #print(f"{m} à uen couv de {len(couv)}")
    return couv,len(couv)

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
        motif.sort
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

if __name__ == "__main__":

    fileNameData=str(sys.argv[1])
    fileNameCT=str(sys.argv[2])
    data=fetch_file(f'datasets/{fileNameData}.dat',int_values=True)
    dataI=dataTransToDataItem(data)
    CTK,CTD=getCT(fileNameData,fileNameCT)
    for i in CTD:
        couv,freqM=getFreq(dataI,i[2])
        print(f"frequence K {i[0]} vs true frequence {freqM}")
       # print(i)