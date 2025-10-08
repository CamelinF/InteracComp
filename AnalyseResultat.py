import os
import numpy as np
import pandas as pd 
import time 
import matplotlib.pyplot as plt
from skmine.datasets.fimi import fetch_file
from sklearn.linear_model import SGDClassifier
from scipy.stats import spearmanr,kendalltau

data=["hepatitis","chess","mushroom","retail","splice1","eisen","pumsb","pumsb_star","connect","weatherAUS","twitter","dota"]
data=["hepatitis"]
# config [c nb Trans,n nb iter,f seuil freq, k taille Requete]
config=[[10,25,20,5]]#,[10,25,20,10]]
temps=["Total","Sample","Krimp","Learn","Req","Weight"]
Paires=["Partiel","Total"]
regularization=["AvecL1"]#"SansL1",
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


## Obtention du temps 

# table_dict=dict()
# first_table_line="datasets;"
# second_table_line=";"
# for c in config:
#     first_table_line+=f"c{c[0]}n{c[1]}f{c[2]}k{c[3]};;"
#     for p in Paires:
#         second_table_line+=f"{p};"

# for t in temps:
#     temps_lines=[]
#     for d in data:
#         #os.system(f"mkdir Res/{d}")
#         plt.figure(figsize=(15, 10))
#         plt.title(f" Evolution temps {t} à chaque itération pour {d}")

#         line_data=f"{d};"
#         for c in config:
#             for p in Paires:
#                 f=fetch_file(f"Res/SizeRequest{c[3]}/Paires_{p}/c{c[0]}n{c[1]}f{c[2]}/{d}/temps{t}.txt",separator=',')
#                 tempsF=[]
#                 for l in f:
#                     for indT in range(len(l)):
#                         if indT==0:
#                             tempsF.append(float(l[indT].split('[')[1]))
#                         else:
#                             if indT==len(l)-1:
#                                 tempsF.append(float(l[indT].split(']')[0]))
#                             else:
#                                 tempsF.append(float(l[indT]))
#                 plt.plot(range(len(tempsF)),tempsF,label=f"{c}-{p}")
#                 line_data+=f"{np.round(sum(tempsF),1)};"
#         temps_lines.append(line_data)
#         plt.legend()
#         #plt.savefig(f"Res/{d}/Gtemps{t}.eps", format='eps', dpi=300, bbox_inches='tight')
#         plt.clf() 
#     table_dict[t]=temps_lines
#     table_dict[t].insert(0,second_table_line)
#     table_dict[t].insert(0,first_table_line)
   
# for i in table_dict:
#     print(type(table_dict[i]))
#     print(type(table_dict[i][0]))
#     f=open(f"Res/temps{i}.txt",'a')
#     for l in table_dict[i]:
#         f.write(f"{l}\n")
#     f.close()

## Check des poids appris 

for d in data:
    print(f"Pour {d}:")
    for regu in regularization:
        print(f" {regu} :")
        for c in config:
            #print(f"    avec la config {c}")
            for p in Paires:
                dicoK=getDecodage(d,True)
                dataI={}
                data=fetch_file(f'data/{d}.dat',int_values=True)
                for t in range(len(data)):
                    for i in data[t]:
                        if i in dataI:
                            dataI[i].add(t)
                        else:
                            dataI[i]={t}
                #print(f"        Get Feedbacks with pairs {p}")
                feedbacks=[]
                f=fetch_file(f"Res/SizeRequest{c[3]}/{regu}/Paires_{p}/c{c[0]}n{c[1]}f{c[2]}/{d}/feedbacks.txt",separator=',')
                for i in f:
                    x=int(i[0].split('[')[1])
                    y=int(i[1].split(']')[0])
                    feedbacks.append([x,y])
                
                f=open(f"Res/SizeRequest{c[3]}/{regu}/Paires_{p}/c{c[0]}n{c[1]}f{c[2]}/{d}/patternShown.txt")

                s=f.readline()
                motifsFeatures=dict()
                while s!='':
                    splitS=s.split('[')
                    ind=int(splitS[1].split(',')[0])
                    splitS=splitS[2].split(']')[0]
                    features=[int(i) for i in splitS.split(',')]
                    features[-1]=features[-1]#/len(dataI)
                    features[-2]=features[-2]#/len(data)
                    motifsFeatures[ind]=features
                    s=f.readline()
                f.close()
                f=open(f"Res/SizeRequest{c[3]}/{regu}/Paires_{p}/c{c[0]}n{c[1]}f{c[2]}/{d}/weightLearned.txt")
                s=f.readline()
                splitS=s.split('[')[1].split(']')[0].split(',')
                weightLearned=[]
                for i in splitS:
                    weightLearned.append(float(i))
                f.close()

                print("             Training with SGD")
                timeTrain=time.time()
                X_pairs=[]
                y_pairs=[]
                for i, j in feedbacks:
                    xi, xj = np.array(motifsFeatures[i]), np.array(motifsFeatures[j])
                    X_pairs.append(xi - xj)   # Si i > j → xi - xj doit avoir score > 0
                    y_pairs.append(1)
                    # On ajoute aussi l'inverse pour équilibrer
                    X_pairs.append(xj - xi)
                    y_pairs.append(0)
                model = SGDClassifier(loss="log_loss", penalty="l1", max_iter=1000)
                model.fit(X_pairs, y_pairs)
                print(f"Temps learn = {time.time()-timeTrain}")
                motifs=[]
                
                print("             Gettings Test set")
                timeAnalyse=time.time()
                allItemsets=[]
                for run in [8,16,24]:#range(c[1]):
                    f=open(f"Res/SizeRequest{c[3]}/{regu}/Paires_{p}/c{c[0]}n{c[1]}f{c[2]}/{d}/Candidates/run{run}.isc")
                    s=f.readline()
                    s=f.readline()
                    s=f.readline()
                    splitS=s.split(":")[1].split()
                    while len(splitS)>2:
                        itemset=[]
                        for ind in range(len(splitS)-1):
                            itemset.append(dicoK[int(splitS[ind])])
                        itemset=sorted(itemset)
                        if not(itemset in allItemsets): 
                            surprise=int(splitS[-1].split('(')[1].split(')')[0])/len(data)
                            estimation=1
                            for i in itemset:
                                estimation*=len(dataI[i])/len(data)
                            surprise-=estimation
                            couv,freq=getFreq(dataI,itemset)
                            itemsetFeatures=[]
                            for i in range(len(data)):
                                if i in couv:
                                    itemsetFeatures.append(1)
                                else:
                                    itemsetFeatures.append(0)    
                            for i in range(len(dataI)):
                                if i in itemset:
                                    itemsetFeatures.append(1)
                                else:
                                    itemsetFeatures.append(0)  
                            itemsetFeatures.append(freq)
                            itemsetFeatures.append(len(itemset))
                            motifs.append([surprise,np.array(itemsetFeatures) @ np.array(weightLearned),itemsetFeatures])
                            allItemsets.append(itemset)
                        s=f.readline()
                        if s!="":
                            splitS=s.split(":")[1].split()
                        else:
                            splitS=[""]
                    f.close()
                
                # print(motifs[0])
                
                print("     Start Testing")
                print(f"        nb Motifs: {len(motifs)}")
                # score ML = X[i] @ clf.coef_.T + clf.intercept_
                # Testing as list
                # scoreTrue=np.array([-motifs[i][0] for i in range(len(motifs))])
                # scoreLCS=np.array([-motifs[i][1] for i in range(len(motifs))])
                # listTrue=np.argsort(scoreTrue)
                # listLCS=np.argsort(scoreLCS)
                # scoreML=[]
                # for m in range(len(motifs)):
                #     x=-(motifs[m][2] @ model.coef_.T+model.intercept_)
                #     scoreML.append(x[0])
                # listML=np.argsort(scoreML)

                # print(f"listTrue :{listTrue[:10]}")
                # print(f"listLCS :{listLCS[:10]}")
                # print(f"listML :{listML[:10]}")
                # rho, _ = spearmanr(listTrue,listLCS)
                # print(f"Spearman LCS : {rho:.4f}")
                # rho, _ = spearmanr(listTrue, listML)
                # print(f"Spearman ML : {rho:.4f}")
                # tau, _ = kendalltau(listTrue, listLCS)
                # print(f"Kendal LCS : {tau:.4f}")
                # tau, _ = kendalltau(listTrue, listML)
                # print(f"Kendall ML : {tau:.4f}")
                
                # # Testing with paires by sublist
                # nbFauxML=0
                # nbVraiML=0
                # nbVraiLCS=0
                # nbFauxLCS=0
                # listTrue=list(listTrue)
                # listLCS=list(listLCS)
                # listML=list(listML)
                # listDeajVu=set()
                # for i in range(len(listTrue)):
                #     motifCheck=listTrue[i]
                    
                #     motifsPlusGrand=set(listTrue[0:i])
                #     motifsPlusPetit=set(listTrue[i+1:len(listTrue)])

                #     indLCS=listLCS.index(motifCheck)
                #     MPGLCS=set(listLCS[0:indLCS])
                #     MPPLCS=set(listLCS[indLCS+1:len(listLCS)])
                #     indML=listML.index(motifCheck)
                #     MPGML=set(listML[0:indML])
                #     MPPML=set(listML[indML+1:len(listML)])
                #     vraiML=(motifsPlusGrand & MPGML).union(motifsPlusPetit & MPPML)-listDeajVu
                #     fauxML=(motifsPlusGrand ^ MPGML).union(motifsPlusPetit ^ MPPML)-listDeajVu
                #     vraiLCS=(motifsPlusGrand & MPGLCS).union(motifsPlusPetit & MPPLCS)-listDeajVu
                #     fauxLCS=(motifsPlusGrand ^ MPGLCS).union(motifsPlusPetit ^ MPPLCS)-listDeajVu
                #     nbVraiML+=len(vraiML)
                #     nbFauxML+=len(fauxML)
                #     nbVraiLCS+=len(vraiLCS)
                #     nbFauxLCS+=len(fauxLCS)
                #     listDeajVu.add(motifCheck)

                # print("nbVraiML: ",nbVraiML)
                # print("nbFauxML: ",nbFauxML)
                # print("nbVraiLCS: ",nbVraiLCS)
                # print("nbFauxLCS: ",nbFauxLCS)
                # print(f" temps passé : {time.time()-timeAnalyse}")
                # Testing with paires
                timeAnalyse=time.time()
                nbFauxML=0
                nbVraiML=0
                nbVraiLCS=0
                nbFauxLCS=0
                
                nbp=np.round(len(motifs)*(len(motifs)+1)/2)
                upd=5
                print(f"nb Motifs : {len(motifs)}")
                print(f"nb Paires à test {nbp}")
                compt=0
                dictPairesPredLCS=dict()
                dictPairesPredML=dict()
                for i in range(len(motifs)):
                    dictPairesPredLCS[i]=[]
                    dictPairesPredML[i]=[]
                for i in range(len(motifs)):
                    for j in range(i+1,len(motifs)):
                        compt+=1
                        if compt % np.round(nbp/upd) ==0:
                            print(f"un {upd}ieme de fait")
                            print(f"    temps passé : {time.time()-timeAnalyse}")
                        y_true=motifs[i][0]>=motifs[j][0]
                        
                        Xi=np.array(motifs[i][2])
                        Xj=np.array(motifs[j][2])
                        # print(f"Xi :{Xi}, {len(Xi)}")
                        # print(f"Xj :{Xj}, {len(Xj)}")
                        
                        pairesTest=Xi-Xj
                        y_predLCS=(np.array(pairesTest) @ np.array(weightLearned))>0
                        if not(y_predLCS):
                             dictPairesPredLCS[i].append(j)
                        else:
                            dictPairesPredLCS[j].append(i)
                        # print(f"Paires test: {pairesTest}, {len(pairesTest)}")
                        y_predML=model.predict([pairesTest])
                        if not(y_predML):
                             dictPairesPredML[i].append(j)
                        else:
                            dictPairesPredML[j].append(i)
                        # fLCS=open(f"Res/SizeRequest{c[3]}/{regu}/Paires_{p}/c{c[0]}n{c[1]}f{c[2]}/{d}/LCSFaux","w")
                        # fML=open(f"Res/SizeRequest{c[3]}/{regu}/Paires_{p}/c{c[0]}n{c[1]}f{c[2]}/{d}/MLFaux","w")
                        if y_true == y_predLCS:
                            nbVraiLCS+=1
                        else:
                            nbFauxLCS+=1
                            # fLCS.write(f" LCS ce trompe sur {i} et {j}")
                        if y_true == y_predML[0]:
                            nbVraiML+=1
                        else:
                            nbFauxML+=1
                            # fML.write(f" ML ce trompe sur {i} et {j}")
                        # fML.close()
                        # fLCS.close()
                scoreTrue=np.array([-motifs[i][0] for i in range(len(motifs))])
                listTrue=np.argsort(scoreTrue)
                listLCS=np.zeros(len(motifs),int)
                listML=np.zeros(len(motifs),int)
                for i in dictPairesPredLCS:
                    listLCS[len(dictPairesPredLCS[i])]=i
                for i in dictPairesPredML:
                    listML[len(dictPairesPredML[i])]=i

                print(f"Pour {d} {regu} et config {c} sur les paires {p} :")
                print(f"    listTrue :{listTrue[:10]}")
                print(f"    listLCS :{listLCS[:10]}")
                print(f"    listML :{listML[:10]}")
                rho, _ = spearmanr(listTrue,listLCS)
                print(f"    Spearman LCS : {rho:.4f}")
                rho, _ = spearmanr(listTrue, listML)
                print(f"    Spearman ML : {rho:.4f}")
                tau, _ = kendalltau(listTrue, listLCS)
                print(f"    Kendal LCS : {tau:.4f}")
                tau, _ = kendalltau(listTrue, listML)
                tau, _ = kendalltau(listTrue, listML)
                print(f"    Kendall ML: {tau:.4f}")

                rho, _ = spearmanr(listTrue[:10],listLCS[:10])
                print(f"    Spearman Top 10LCS : {rho:.4f}")
                rho, _ = spearmanr(listTrue[:10], listML[:10])
                print(f"    Spearman top 10 ML : {rho:.4f}")
                tau, _ = kendalltau(listTrue[:10], listLCS[:10])
                print(f"    Kendal Top 10 LCS : {tau:.4f}")
                tau, _ = kendalltau(listTrue[:10], listML[:10])
                print(f"    Kendall Top 10 ML : {tau:.4f}")
                print("    nbVraiML: ",nbVraiML)
                print("    nbFauxML: ",nbFauxML)
                print("    nbVraiLCS: ",nbVraiLCS)
                print("    nbFauxLCS: ",nbFauxLCS)
                print(f"temps passé : {time.time()-timeAnalyse}")
                