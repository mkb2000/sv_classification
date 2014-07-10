'''
merge arff files into a new file
'''
import sys
if __name__=="__main__":
    if len(sys.argv)<2:
        print "Usage: "+sys.argv[0]+" arff1 arff2 [arff3...]"
        exit()
    
    arffs=sys.argv[1:]
    outDIR='/'.join(sys.argv[1].split('/')[0:-1])
#     print outDIR
#     exit()
    
    preArffhead=''
    data=''
    for f in arffs:
        content=open(f).readlines()
        arffhead=''
        for l in content:
            if l[0]=='@':
                arffhead+=l
            else:
                data+=l
        
        #check if head are consistent
        if preArffhead=='':
            preArffhead=arffhead
        else:
            if preArffhead!=arffhead:
                print "arffheads are not consistent at"+f
                exit()
    
    #out put merged file
    outfile=open(outDIR+'/'+'merged.arff','w')
    outfile.write(preArffhead+data)
    outfile.close()
                
                
                
