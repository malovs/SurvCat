#!/usr/bin/python
# SNP generator for database
n=1500
N=10
q=5
d=0.1
lw=0.1  
up=0.5
proj=3 # lengh of the project
deg=0.3 # defect at follow up
mfa=0.2 # missing at follow up intensity
lmbase=0.05
lmmaxD=0.2
lmmaxCD=0.1

#variables
X = list( list(0 for a in range(0, n)) for b in range(0,N)) #  SNP1
Y = list( list(0 for a in range(0, n)) for b in range(0,N)) # SNP2
T = list(0 for a in range(0,n)) #  Time of failure
dl = list(0 for a in range(0,n)) #  Time of censoring
# functions
def rvar(y,p):
        u=0
        if y<p:
            u=1
        return u

# generator
import random
import math
lmdifD=lmmaxD-lmbase
lmdifCD=lmmaxCD-lmbase
delta = list(rvar(random.random(),q/float(N)) for a in range(0,N))
for j in range(n):
        for i in range(N):
                p=random.uniform(lw,up)
                lm=lmbase
                u=random.random()
                X[i][j]=rvar(u,p)
                u=random.random()
                Y[i][j]=rvar(u,p)
                if delta[i]==1:
                        if (X[i][j]==1 and Y[i][j]==1):
                                lm<-lm+lmdifD/q
                        elif (X[i][j]==1 and Y[i][j]==1):
                                lm<-lm+lmdifCD/q                        
        t1=-math.log(random.random())/lm
        u1=-math.log(random.random())/mfa
        if u1>proj:
                lww=proj*(1-deg)
                u1=random.uniform(lww,proj)
        if t1<=u1:
                T[j]=int(t1*365)
                dl[j]=1
        else:
                T[j]=int(u1*365)
                dl[j]=0
# for database
f1=open("/Users/malovs/Desktop/DATASRV.txt",'w')
f1.write("ID,TIME,Ind \n")

for i in range(n):
        f1.write(str(i+1)+","+str(T[i])+","+str(dl[i])+"\n")
f1.close()

f2=open("/Users/malovs/Desktop/DATASNP.txt",'w')
f2.write("ID,SNPID,SNPA1,SNPA2 \n")
for i in range(n):
        for j in range(N):
                f2.write(str(i+1)+","+str(j+1)+","+str(X[j][i])+","+str(Y[j][i])+"\n")                
f2.close()


