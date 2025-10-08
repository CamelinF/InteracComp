import os

data=["hepatitis","chess","mushroom","retail","splice1","eisen","pumsb","pumsb_star","connect","weatherAUS","twitter"]

for i in data:
    print(i)
    os.system(f"python3 LCS-IH.py {i} 10 25 0.2")