nA=186
nB=38
nAB=13
nO=284
n=nA+nB+nAB+nO

#Step 1:Initialize the parameters
p=0.3  # allele frequencies of the A
q=0.2  # allele frequencies of the B
r=0.5  # allele frequencies of the O
p;q;r

# Step 2: (E) :Calculate approximate values for nAA, nA0, nBB , nB0  
exp1<-function(p,q,r){
    nAA=nA*p/(p+2*r)
    nAO=(nA*2*r)/(p+2*r)
    nBB=nB*q/(q+2*r)
    nBO=(nB*2*r)/(q+2*r)
    c(nAA,nAO,nBB,nBO)
}
#Visualize nAA, nAO,nBB,nBO obtained in the first iteration:
nAA=exp1(p,q,r)[1]
nAO=exp1(p,q,r)[2]
nBB=exp1(p,q,r)[3]
nBO=exp1(p,q,r)[4]

#Step 3: (M):Update p, q and r 
param=function(nAA,nAO,nBB,nBO){
    p=(2*nAA+nAO+nAB)/(2*n)
    q=(2*nBB+nBO+nAB)/(2*n)
    r=(2*nO+nAO+nBO)/(2*n)
    c(p,q,r)
}
#Visualize ˆp, ˆq and ˆr, obtained in the first iteration:
param(nAA,nAO,nBB,nBO)

#Steps 4 and 5:Iterative procedure 
i=0; er=1
while(sum(er>=0.00001)>0){
  e=exp1(p,q,r)                           # Step E
  par=param(e[1],e[2],e[3],e[4])            # Step M
  er=abs(c(p,q,r)-c(par[1],par[2],par[3])) # STOP criteria
  i=i+1
  p=par[1]; q=par[2]; r=par[3]
  cat(i,p,q,r,"\n")           
}
cat(" estimated alle frequency of A=",par[1],"\n",
    "estimated alle frequency of B=",par[2],"\n",
    "estimated alle frequency of C=",par[3])



