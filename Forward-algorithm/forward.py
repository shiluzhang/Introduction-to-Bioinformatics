####################################################################################################
## Forward algorithm
import argparse
import numpy as np


    
def main(args):   
    n = args.n
    transitions_file = args.transitions_file
    emissions_file = args.emissions_file
    sequence = args.sequence
    
    #transition matrix
    a=np.zeros((n+2,n+2))
    e=np.zeros((n+2,4))
        
    with open(transitions_file, "r") as in_file:
        for line in in_file:
            ## representing K, r1, r2, g and s,
            ## Add one check if there is no \n at the end of line:
            lin= line.strip()
            #print(lin)
            number_string = lin.split(' ')
            #number_string = [int(i) for i in number_string]
            a[int(number_string[0])][int(number_string[1])]=float(number_string[2])
                
    #print("Emission file:")
    with open(emissions_file, "r") as in_file:
        for line in in_file:
            ## representing K, r1, r2, g and s,
            ## Add one check if there is no \n at the end of line:
            lin= line.strip()
            #print(lin)
            number_string = lin.split(' ')
            if number_string[1]=='A':
                e[int(number_string[0])][0]=float(number_string[2])
            elif number_string[1]=='C':
                e[int(number_string[0])][1]=float(number_string[2])    
            elif number_string[1]=='G':
                e[int(number_string[0])][2]=float(number_string[2])
            elif number_string[1]=='T':
                e[int(number_string[0])][3]=float(number_string[2])                        
        

    #print(a)    
    #print(e)    
    #print(sequence)
    sequenceNum=[]
    for s in sequence:
        if s=='A':
            sequenceNum.append(0)
        elif s=='C':
            sequenceNum.append(1)  
        elif s=='G':
            sequenceNum.append(2)  
        elif s=='T':
            sequenceNum.append(3)         
   
    #print(sequenceNum)        
    # initialization
    #f[k][i] be the probability of generating the first i characters and ending in state k    
    # 0:begin  n+1:end
    L=len(sequence)
    f=np.zeros((n+2,L+1))
    f[0][0]=1
    
    # emitting states:
    for i in range(1,L+1):
        for l in range(1,n+1):
            score=0          
            for k in range(0,n+2):
                score=score+f[k][i-1]*a[k][l]            
            
            #print('l:',l,' i:',i,' seqpos: ',sequenceNum[i-1])
            f[l][i]=e[l][sequenceNum[i-1]]*score
        
        #silent state:
        l=n+1
        score=0          
        for k in range(0,n+2):
            score=score+f[k][i]*a[k][l]
            
        f[l][i]=score
            
    #print(f)
    
    for i in range(1,L+1):
        for k in range(1,n+1):
            if f[k][i]==0:
                logprob=-1000000
            else:
                logprob=np.log(f[k][i])
            print("%.2f" %logprob,end='\t')
            
        print()
    
    print("%.2f" %np.log(f[n+1][L]))
    
    #print(''.join(alignseq1))
    #print(''.join(alignseq2))

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(

        description=__doc__,

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--transitions_file',

                        help='transitions file',

                        type=str,

                        default='') 
    
    parser.add_argument('--emissions_file',

                        help='emissions file',

                        type=str,

                        default='')     
    
    parser.add_argument('--n',

                        help='n.',

                        type=int,

                        default=6)
    
    parser.add_argument('--sequence',

                        help='sequence.',

                        type=str,

                        default='')    
    
    args = parser.parse_args()
    main(args)
    
####################################################################################################
    