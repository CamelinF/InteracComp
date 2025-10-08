import os

data=["hepatitis","chess","mushroom","retail","splice1","eisen","pumsb","pumsb_star","connect","weatherAUS","twitter"]
data=["hepatitis"]
for i in data:
    print(i)
    # nombre de transactions C=10,  nombre d'itération n=25, Fréquence dans les transactions échantillonné f= 0.2 , nombre de motifs à présenter k=5
    os.system(f"python3 LCS-IH.py {i} 10 25 0.2 5")