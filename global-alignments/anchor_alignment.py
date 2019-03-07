####################################################################################################
## euler assemble, that takes as input a set of k-mers and outputs a single super string
## that has a k-mer spectrum equal to the set of input k-mers.
import argparse
import numpy as np


def scorefun(x,y):
    if x==y:
        score=r1
    else:
        score=-r2
    return(score)
 
def gappenalty(k,g,s):
    if k==0:
        w=0
    elif k>=1:
        w=g+s*k
    return(s)

def Mscore(i,j,x,y):
    c1=M[i-1,j-1]+scorefun(x,y)
    c2=Ix[i-1,j-1]+scorefun(x,y)
    c3=Iy[i-1,j-1]+scorefun(x,y)
    scores=[c1,c2,c3]
    ### we only need first occurence:
    #id=np.argmax(scores)
    #key=str(i)+'-'+str(j)
    #Mp[key]=id
    winner = np.argwhere(scores == np.amax(scores))
    id=winner.flatten().tolist()    
    key=str(i)+'-'+str(j)
    Mp[key]=[]
    for ix in range(0,len(id)):
        Mp[key].append(id[ix])
        
    return(max(scores))

def Ixscore(i,j):
    c1=M[i-1,j]+g+s
    c2=Ix[i-1,j]+s
    c=max([c1,c2])
    key=str(i)+'-'+str(j)
    Ixp[key]=[]
    if c1>c2:
        Ixp[key].append(0)
    elif c1<c2:
        Ixp[key].append(1)
    else:
        Ixp[key].append(0)
        Ixp[key].append(1)
        
    return(c)

def Iyscore(i,j):
    c1=M[i,j-1]+g+s
    c2=Iy[i,j-1]+s
    c=max([c1,c2])
    key=str(i)+'-'+str(j)
    Iyp[key]=[]
    if c1>c2:
        Iyp[key].append(0)
    elif c1<c2:
        Iyp[key].append(2)
    else:
        Iyp[key].append(0)
        Iyp[key].append(2)        
    return(c)

def getPointer(id):
    if id==0:
        return(Mp)
    elif id==1:
        return(Ixp)
    elif id==2:
        return(Iyp)
    
def globalalign(x,y):
    global M,Ix,Iy
    M=np.zeros((m+1,n+1))
    Ix=np.zeros((m+1,n+1))
    Iy=np.zeros((m+1,n+1))
    global Mp,Ixp,Iyp
    Mp={}
    Ixp={}
    Iyp={}
    # initialization:
    M[0,0]=0
    Ix[0,0]=g
    Iy[0,0]=g
    for i in range(1,m+1):
        Ix[i,0]=g+s*i
        Iy[i,0]=-float('inf')   
        M[i,0]=-float('inf')   
        key=str(i)+'-'+str(0)
        Ixp[key]=[]
        Ixp[key].append(1)
        
    for j in range(1,n+1):
        Iy[0,j]=g+s*j    
        Ix[0,j]=-float('inf') 
        M[0,j]=-float('inf')  
        key=str(0)+'-'+str(j)
        Iyp[key]=[]
        Iyp[key].append(2)        
     
    ## fill in the matrix:
    for i in range(1,m+1):
        for j in range(1,n+1):
            M[i,j]=Mscore(i,j,x[i-1],y[j-1])
            Ix[i,j]=Ixscore(i,j)
            Iy[i,j]=Iyscore(i,j)
     
     
    #print(M)
    #print(Ix)
    #print(Iy)
    #print(Mp)
    #print(Ixp)
    #print(Iyp)    
    
    ## traceback:
    alignseq1=[]
    alignseq2=[]
    kx=m-1
    ky=n-1
    i=m
    j=n
    id=np.argmax([M[i,j],Ix[i,j],Iy[i,j]])
    #print("index: ",id)
    key=str(i)+'-'+str(j)
       
    while i>=0 and j>=0:
        if i==0 and j==0:
            break
        ## only need the first occurence:
        key=str(i)+'-'+str(j)
        #print("key: ",key)
        ## choose which gives the max score:
        Pointer=getPointer(id)
        nextid=Pointer[key][0]  ## only need first one        
        if id==0:
            ## match xi with yi
            alignseq1.append(x[i-1])
            alignseq2.append(y[j-1])
            ## pointer to M[i-1,j-1] Ix[i-1,j-1],Iy[i-1,j-1]]           
            i=i-1
            j=j-1
        elif id==1:
            ## insertion of y
            alignseq1.append(x[i-1])
            alignseq2.append('_')
            i=i-1
        elif id==2:
            ## insertion of x
            alignseq1.append('_')
            alignseq2.append(y[j-1])
            j=j-1            
        id=nextid
      
    outseq=[]    
    outseq.append(''.join(alignseq1[::-1]))
    outseq.append(''.join(alignseq2[::-1]))
    #print(outseq[0])
    #print(outseq[1])
    return(outseq)
        
def main(args):   
    data_file_path = args.inputfile
    ## global variables:
    global K,r1,r2,g,s,m,n    
    K=0
    k=0
    seq=[]
    Xanchor=[]
    Yanchor=[]
    with open(data_file_path, "r") as in_file:
        for line in in_file:
            ## representing K, r1, r2, g and s,
            ## Add one check if there is no \n at the end of line:
            if line[-1]=="\n":
                lin=line[:-1]
            else:
                lin=line
            #print(lin)
            if k==0:
                number_string = lin.split(' ')
                number_string = [int(i) for i in number_string]
                K=number_string[0]
                r1=number_string[1]
                r2=number_string[2]
                g=-number_string[3]
                s=-number_string[4]                
                k=k+1
            elif k==1 or k==2:
                seq.append(lin)
                k=k+1
            elif K>0 and k>=3:
                number = lin.split(' ')
                number = [int(i) for i in number]           
                ## start from position 0:
                Xanchor.append(number[0]-1)
                Yanchor.append(number[1]-1)
                

    
    ##print(seq)
    ##print("K:",K," r1: ",r1," r2: ",r2," g: ",g," s: ",s)
    alignseq1=[]
    alignseq2=[]
    
    if K==0:
        m=len(seq[0])
        n=len(seq[1])        
        outseq=globalalign(seq[0],seq[1])
        alignseq1.append(outseq[0])
        alignseq2.append(outseq[1])         
    elif K>0:
        i=0
        seq1=seq[0][0:Xanchor[i]]
        seq2=seq[1][0:Yanchor[i]]
        m=len(seq1)
        n=len(seq2)                 
        outseq=globalalign(seq1,seq2)
        #print(outseq)
        alignseq1.append(outseq[0])
        alignseq2.append(outseq[1])             
        for i in range(0,K-1):                   
            alignseq1.append(seq[0][Xanchor[i]])
            alignseq2.append(seq[1][Yanchor[i]])            
            seq1=seq[0][(Xanchor[i]+1):Xanchor[i+1]]
            seq2=seq[1][(Yanchor[i]+1):Yanchor[i+1]]      
            m=len(seq1)
            n=len(seq2)                  
            outseq=globalalign(seq1,seq2)
            #print(outseq)
            alignseq1.append(outseq[0])
            alignseq2.append(outseq[1])                
        
        i=K-1
        alignseq1.append(seq[0][Xanchor[i]])
        alignseq2.append(seq[1][Yanchor[i]])            
        seq1=seq[0][(Xanchor[i]+1):len(seq[0])]
        seq2=seq[1][(Yanchor[i]+1):len(seq[1])]      
        m=len(seq1)
        n=len(seq2)                  
        outseq=globalalign(seq1,seq2)
        #print(outseq)
        alignseq1.append(outseq[0])
        alignseq2.append(outseq[1])  
        
    print(''.join(alignseq1))
    print(''.join(alignseq2))

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(

        description=__doc__,

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--inputfile',

                        help='input file',

                        type=str,

                        default='')    
    
    args = parser.parse_args()
    main(args)
    
####################################################################################################
    
