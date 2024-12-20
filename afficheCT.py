from skmine.datasets.fimi import fetch_file



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

def getCT(fileNamD,fileNamCT):
    dicoK= getDecodage(fileNamD,True)
    f=open(fileNamCT)
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
        tmpp=int(tmp.split("(")[1])
        tmp=splitS[-1].split(",")[1]
        freq=int(tmp.split(")")[0])
        #couv,freq=getFreq(dataI,motif)
        if len(motif)>1:
            ensMotif.append([freq,len(itemset),itemset])


        sumTot+=tmpp
        poids.append(tmpp)
        ensCTFC.append(motif)
        ensCT.append([motif,freq])
        s=f.readline()
        splitS=s.split()
    f.close()
    return CT


if __name__ == "__main__":
    fileNameData=str(sys.argv[1])
    fileNameCT=int(sys.argv[2])

    data=fetch_file(f'datasets/{fileNameData}.dat',int_values=True)