gras.func<-function(X0,maxiter=5000,eps=0.1e-5,autofilling=FALSE){
# PURPOSE: estimate a new matrix X with exogenously given row and column 
# totals that is a close as possible to a given original matrix X0 using 
# the Generalized RAS (GRAS) approach. The function is adapted to give 
# warnings for most common sources of errors for I-O practisioners
# Adapted R version from Matlab 
# REFERENCES:
# Temurshoev, U., R.E. Miller and M.C. Bouwmeester (2013), A note on the
# GRAS method, Economic Systems Research, 25, pp. 361-367.
# ------------------------------------------------------------------------
# USAGE: X = gras(X0) with or without maxiter, eps and autofilling
# INPUT:
# -> X0 = benchmark (base) matrix, not necessarily square, with row and column
#  totals in the last column and row
# -> maxiter = maximum number of iterations after which the algorithm returns 
# aprroximate solution to  speed up calculations
# -> eps = convergence tolerance level; if empty, the default threshold 
# is 0.1e-5 (=0.000001)
# -> autofilling = fills initial values for zero base matrix rows or columns 
# if row or column totals are non zero

# OUTPUT:
  # -> X = estimated/adjusted/updated matrix

# -------------------------------------------------------------------------
# Original Matlab function written by:   
#   Umed Temurshoev, 07/10/2010 with later adjustments
#   Current e-mail: umed.temurshoev@ec.europa.eu
# Adapted for use in R by Peter Horvat
# email: peter.horvat@oecd.org
# -------------------------------------------------------------------------
# example:
# X0<-c(7,2,-2,3,9,0,5,8,2,-3,1,1)
# dim(X0)<-c(3,4)
# u<-t(t(c(15,26,-1)))
# v<-c(9,16,17,-2)
# X0<-cbind(rbind(X0,v),rbind(u,0))
# gras.func(X0)
# [,1]      [,2]     [,3]       [,4]
# [1,]  8.976439  3.743161 5.721729 -3.4413294
# [2,]  2.799334 12.256841 9.992313  0.9515109
# [3,] -2.775776  0.000000 1.285958  0.4898178
#

X0<-as.matrix(X0)
u<-t(t(X0[-dim(X0)[1],dim(X0)[2]]))
v<-X0[dim(X0)[1],-dim(X0)[2]]
X0<-X0[-dim(X0)[1],-dim(X0)[2]]
rownamesX0<-dimnames(X0)[[1]]
colnamesX0<-dimnames(X0)[[2]]

if(autofilling){
  if(any(colSums(X0)==0))
    if(any(v[which(colSums(X0)==0)]!=0)){
        X0[,v[which(colSums(X0)==0)]!=0]<-min(abs(X0[X0!=0]))/1000
      }

  if(any(rowSums(X0)==0))
      if(any(u[which(rowSums(X0)==0)]!=0)){
          X0[u[which(rowSums(X0)==0)]!=0,]<-min(abs(X0[X0!=0]))/1000
      }
}
for(i in 1:length(u))
  if(u[i]==0)
    if(sum(X0[i,]!=0))
       cat("Zero constraint, but intermediate matrix contains non-zero value for row:",i,"\n")

for(i in 1:length(v))
   if(v[i]==0)
     if(sum(X0[,i]!=0))
        {cat("Zero constraint, but intermediate matrix contains non-zero value for column:",i,"\n")
      #print(X0[,40:44]) 
     }
              

invd<-function(x){
    invd = 1./x
    invd[x==0] = 1;
    invd = diag(c(as.matrix(invd)))
    return(invd)
}        

library(Matrix)
m<- dim(X0)[1]
n<- dim(X0)[2]
N<-matrix(0,nrow=m,ncol=n)

N[X0<0]<--X0[X0<0]
if(sum(N)!=0){
    r0<-NULL
    c0<-NULL
    x0<-NULL
    for(i in 1:dim(N)[1])
        for(j in 1:dim(N)[2]){
            if(N[i,j]>0){
                r0<-c(r0,i)
                c0<-c(c0,j)
                x0<-c(x0,N[i,j])
            }
        }
        


    N = sparseMatrix(i=r0,j=c0,x=x0,dims=dim(X0))      #could save memory with large-scale matrices
    P = X0+N
    r1<-NULL
    c1<-NULL
    x1<-NULL
    for(i in 1:dim(N)[1])
        for(j in 1:dim(N)[2]){
            if(P[i,j]>0){
                r1<-c(r1,i)
                c1<-c(c1,j)
                x1<-c(x1,P[i,j])
            }
        }
    P = sparseMatrix(i=r1,j=c1,x=x1,dims=dim(X0))
}
P = X0+N
r = matrix(1,nrow=m,ncol=1)      #initial guess for r (suggested by J&O, 2003)
pr = t(P)%*%r

nr = t(N)%*%invd(r)%*%matrix(1,nrow=m,ncol=1)
s1 = invd(2*pr)%*%(v+sqrt(v^2+4*pr*nr))    #first step s
ss = -invd(v)%*%nr
s1[pr==0] = ss[pr==0]
ps = P%*%s1
ns = N%*%invd(s1)%*%matrix(1,nrow=n,ncol=1)
r = invd(2*ps)%*%(u+sqrt(u^2+4*ps*ns))     #first step r
rr = - invd(u)%*%ns
r[ps==0] = rr[ps==0]

pr = t(P)%*%r
nr = t(N)%*%invd(r)%*%matrix(1,nrow=m,ncol=1)
s2 = invd(2*pr)%*%(v+sqrt(v^2+4*pr*nr))    #second step s  
ss = -invd(v)%*%nr
s2[pr==0] = ss[pr==0]
dif = s2-s1
iter = 1                #first iteration
M = max(abs(dif))
while (M > eps){
    s1 = s2
    ps = P%*%s1
    ns = N%*%invd(s1)%*%matrix(1,nrow=n,ncol=1)
    r = invd(2*ps)%*%(u+sqrt(u^2+4*ps*ns))  #previous step r
    rr = -invd(u)%*%ns
    r[ps==0] = rr[ps==0]
    pr = t(P)%*%r
    nr = t(N)%*%invd(r)%*%matrix(1,nrow=m,ncol=1)
    s2 = invd(2*pr)%*%(v+sqrt(v^2+4*pr*nr))  #current step s
    ss = -invd(v)%*%nr
    s2[pr==0] = ss[pr==0]
    dif = s2-s1
    iter = iter+1
    M = max(abs(dif))
    if(iter==maxiter) break
    if(M>1e100){
      stop(paste0("Algorithm is diverging from solution, check input data (mainly constraints) if it is possible to find a solution, problematic column is: ",dimnames(X0)[[2]][which(dif==M)] ))
    } 
}
s = s2                                        #final step s
ps = P%*%s
ns = N%*%invd(s)%*%matrix(1,nrow=n,ncol=1)
r = invd(2*ps)%*%(u+sqrt(u^2+4*ps*ns))        #final step r
rr = -invd(u)%*%ns
r[ps==0] = rr[ps==0]
X = diag(c(as.matrix(r)))%*%P%*%diag(c(as.matrix(s)))-invd(r)%*%N%*%invd(s)       #updated matrix
dimnames(X)<-list(rownamesX0,colnamesX0)
return(as.matrix(X))
}
