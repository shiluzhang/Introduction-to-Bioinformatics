####################################################################################################
## euler assemble, that takes as input a set of k-mers and outputs a single super string
## that has a k-mer spectrum equal to the set of input k-mers.
import argparse

class Node(object):
    def __init__(self, node):
        self.node=node
        self.indf=0
        self.outdf=0     
    
    def updateindf(self,indf):
        self.indf=indf
        
    def updateoutdf(self,outdf):
        self.outdf=outdf    
    
    
class Graph(object):
    """Graph data structure, directed"""

    def __init__(self, node, connections):
        self.graph = {}
        self.node={}
        self.balance=True
        self.start=node[0]
        for items in node:
            self.node[items]=Node(items)
            self.graph[items]=[]
            
        self.add_connections(connections)
        ## need to sort the items in dictionary
        for key in node:
            self.graph[key]=sorted(self.graph[key])
        
        for items in self.node:
            #print(self.node[items].node," in df ",self.node[items].indf, " out df ",self.node[items].outdf)
            if self.node[items].outdf>self.node[items].indf:
                self.start=self.node[items].node
                self.balance=False
            elif self.node[items].outdf<self.node[items].indf:
                self.end=self.node[items].node
                self.balance=False
                    
        ## add fake edges:
        if not self.balance:
            self.add(self.end,self.start)
            #print("The graph is not balanced")
            #print("Add fake edge: ",self.end,"->",self.start)

    def add_connections(self, connections):
        """ Add connections (list of tuple pairs) to graph """

        for node1, node2 in connections:
            self.add(node1, node2)
            
            
    def add(self, node1, node2):
        """ Add connection between node1 and node2 """
        self.graph[node1].append(node2)
        Node1=self.node[node1]
        self.node[node1].updateoutdf(Node1.outdf+1)
        Node2=self.node[node2]
        self.node[node2].updateindf(Node2.indf+1)
        
    def remove(self, node1, node2):
        """ Remove connection between node1 and node2 """
        self.graph[node1].remove(node2)
        if len(self.graph[node1])==0:
            del self.graph[node1]
        
        

    
def EulerianCycle(data_file_path):
    kmer=[]  
    connections = []
    node=set()
    
    with open(data_file_path, "r") as in_file:
        for line in in_file:
            k=line[:-1]
            #print(k)
            kmer.append(k)
            n=len(k)
            node1=k[0:(n-1)]
            node2=k[1:]
            node.add(node1)
            node.add(node2)
            connections.append((node1,node2))                        
    
    #print("Kmer:")
    #for k in kmer:
        #print("(%s->%s)" %(k[0:len(k)-1],k[1:]))
                
    node=sorted(node)    
    #print("Nodes are: ", node)        
    g = Graph(node, connections)
    #print(g.graph)
    
    cyclestart=node[0]
    ## need to find start point: start: out-in=1, end: in-out=1
    graphcycle=[]    
    cycle=0
    cyclestartlist=[]
    cyclestartlist.append(cyclestart)
    cyclestartid=0
    while len(g.graph)>0:
        #print("Start: Cycle ",cycle, " ",g.graph)
        #print("Cyclestart ",cyclestart)
        graphpath=[]
        graphpath.append(cyclestart)        
        nextkey=g.graph[cyclestart][0]
        graphpath.append(nextkey)
        g.remove(cyclestart,nextkey)
        #print("Current Graph: Cycle ",cycle, " ",g.graph)
        while nextkey!=cyclestart:
            node1=nextkey
            if nextkey not in g.graph.keys():
                break
            else:
                node2=g.graph[nextkey][0]
                graphpath.append(node2)
                g.remove(node1,node2)
                #print(g.graph)
                nextkey=node2
            
        if cycle==0:
            graphcycle=graphpath
            #print("cycle " , str(cycle) ," " ,graphpath)
        else:
            #print("cycle " , str(cycle) ," " ,graphpath)
            #insertpos=graphcycle.index(cyclestart)
            insertpos=cyclestartid
            for i in range(1,len(graphpath)):
                graphcycle.insert(i+insertpos,graphpath[i])
                                
        #print("Graph Path: ",graphcycle, " cycle ",str(cycle))     
        ## find the next start with unused outgoing edges:
        while len(g.graph)>0 and cyclestart not in g.graph.keys():
            cyclestartid=cyclestartid+1
            #graphid=graphcycle.index(cyclestart)+1
            if cyclestartid==len(graphcycle):
                break
            else:
                cyclestart=graphcycle[cyclestartid] ##empty
                #print("cyclestart ",cyclestart)
        
        cycle=cycle+1    
            
            
    #print(graphcycle)
    euler=[]
    ## remove fake edges:
    if g.balance==False:
        for i in range(0,len(graphcycle)-1):
            if graphcycle[i]==g.end and graphcycle[i+1]==g.start:
                ## this is the fake edge that comes first:
                id=i+1
                #print(id)
                break
            
        ## reorder:
        for i in range(id,len(graphcycle)-1):  ## no need to put the last one
            euler.append(graphcycle[i])
            
        for i in range(0,id):
            euler.append(graphcycle[i])        
        
    else:
        id=0
        euler=graphcycle
    
    
    #print(euler)
    eulerpath=[]   
    eulerpath.append(euler[0])
    for i in range(1,len(euler)):
        eulerpath.append(euler[i][len(euler[i])-1])

    
    ## the sequence is:    
    print(''.join(eulerpath))

    
def main(args):   
    data_file_path = args.kmerfile
    #data_file_path = "3_10.kmers.txt"
    EulerianCycle(data_file_path)

    

        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(

        description=__doc__,

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--kmerfile',

                        help='Kmer input file',

                        type=str,

                        default='')    
    
    args = parser.parse_args()
    main(args)
    
####################################################################################################
    