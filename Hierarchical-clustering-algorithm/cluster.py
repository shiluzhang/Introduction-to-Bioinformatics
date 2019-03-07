####################################################################################################
## bottom-up hierarchical clustering algorithm 
import argparse
import numpy as np
import math
def Euclid(x,y):
    dist=0
    for i in range(len(x)):
        dist=dist+math.pow((x[i]-y[i]),2)
    return(dist)

def SingleLink(D,C):
    Dp=np.zeros((len(C),len(C)))
    for i in range(len(C)-1):
        for j in range(i+1,len(C)):
            dist=np.inf
            for a in C[i]:
                for b in C[j]:
                    #print(C[i],C[j],a,b)
                    if D[a,b] < dist:
                        dist=D[a,b]
                        
            Dp[i,j]=dist 

    return(Dp)

def CompleteLink(D,C):
    Dp=np.zeros((len(C),len(C)))
    for i in range(len(C)-1):
        for j in range(i+1,len(C)):
            dist=0
            for a in C[i]:
                for b in C[j]:
                    #print(C[i],C[j],a,b)
                    if D[a,b] > dist:
                        dist=D[a,b]
                        
            Dp[i,j]=dist 
    return(Dp)

def AverageLink(D,C):
    Dp=np.zeros((len(C),len(C)))
    #mindist=np.inf
    for i in range(len(C)-1):
        for j in range(i+1,len(C)):
            dist=0
            for a in C[i]:
                for b in C[j]:
                        dist=dist+D[a,b]
                        
            dist=dist/(len(C[i])*len(C[j]))
            Dp[i,j]=dist 
    return(Dp)

def distancefun(link_type,D,C):
    if link_type=='S':
        #print('Single Link')
        Dp=SingleLink(D,C)
    elif link_type=='C':
        Dp=CompleteLink(D,C)
    elif link_type=='A':
        Dp=AverageLink(D,C)
    return(Dp)

def hierarchical(D,k,link_type):
    c={}
    C=[]
    n=D.shape[0]
    Dp=D
    for i in range(0,n):
        c[i]=[]
        c[i].append(i)
        C.append(c[i])
    #print(C)    
    j=n-1
    while len(C)>k:
        j=j+1
        mindist=np.inf
        a=0
        b=0
        for u in range(Dp.shape[1]-1):
            for v in range(u+1,Dp.shape[1]):
                    if  Dp[u,v] < mindist:
                        mindist=Dp[u,v]
                        a=u
                        b=v

        c[j]=C[a]+C[b]
        C.remove(C[a])
        C.remove(C[b-1])
        C.append(c[j])
        Dp=distancefun(link_type,D,C)
        #print(Dp)
        #print(C)  
        
    return(C)
        
def main(args):   
    data_file = args.data_file
    link_type = args.link_type
    k = args.k
    measure=[]  
    genename=[]
    geneid=[]
    with open(data_file, "r") as in_file:
        for line in in_file:
            lin= line.strip()
            #print(lin)
            value = lin.split('\t')
            geneid.append(value[0])
            genename.append(value[1])
            n=len(value)-2
            for i in range(2,len(value)):
                #print(value[i])
                measure.append(float(value[i]))

    m=int(len(measure)/n)
    data=np.array(measure).reshape(m,n);
    #print(data[0,:])    
    D=np.zeros((m,m))
    for i in range(m):
        for j in range(i+1,m):
            D[i,j]=Euclid(data[i,:],data[j,:])
            D[j,i]=D[i,j]
    
    #print(D)        
    #print(link_type)
    C=hierarchical(D,k,link_type)
    dmean=np.mean(data,axis=1)
    #print(dmean)
    id=np.argsort(dmean)
    #print(id)
    cmean=np.zeros(k)
    for i in range(k):
        for j in C[i]:
            cmean[i] = cmean[i]+dmean[j]
        cmean[i]=cmean[i]/len(C[i])
    
    #print(cmean)   
    cid=np.argsort(cmean)
    #print(cid)
    for i in cid:    
        for j in id:
            if j in C[i]:
                print(geneid[j],end='\t')
                print(genename[j],end='\t')
                print("%.3f" %dmean[j],end='\n')
        print("%.3f" %cmean[i])
        print()


    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(

        description=__doc__,

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--data_file',

                        help='data file',

                        type=str,

                        default='') 
    
    parser.add_argument('--link_type',

                        help='link type: S or C or A',

                        type=str,

                        default='')     
    
    parser.add_argument('--k',

                        help='k',

                        type=int,

                        default=6)
       
    
    args = parser.parse_args()
    main(args)
    
####################################################################################################
    