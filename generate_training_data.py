'''
Not used. Obsolete.

Created on May 26, 2014

@author: HELE
'''
import sys
import glob,os
if __name__=='__main__':
    if len(sys.argv)<3:
        print "Usage: "+sys.argv[0]+" [training data amount] [folder of raw files]"
        exit()
    amount=sys.argv[1]
    each_amount=int(amount)/7
    folder=sys.argv[2]
    
    types={"INV":0,"INS":0,"DEL":0,"FUS":0,"TAN":0,"TRAN":0,"WOFS":0}
    
    header=58 #first 58 lines are header
    
    arffhead=[]
    collection=[]
    final_result=[]
    
    os.chdir(folder)
    for f in glob.glob("*.arff"):
        ctemp=open(f).readlines()
        arffhead=ctemp[0:58]
        collection+=ctemp[58:]
    
    for l in collection:
        ll=l.strip().split(",")
        if types[ll[-1]]<each_amount:
            types[ll[-1]]+=1
            final_result.append(l)
    
    outf=open(amount+".arff","w")
    for l in arffhead:
        outf.write(l)
    for l in final_result:
        outf.write(l);
    
    