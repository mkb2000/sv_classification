'''
Created on Jan 10, 2014

@author: HELE

change the position of original breakpoint file:svtype length breakpoint1 breakpoint2 [breakpoint3 breakpoint4]
merge two line events:INV FUS
'''
import sys
from svDetectFuncs import delimeter

def eventFromBreakpointFile(bkf,outputfile):
    breakpointfile=bkf
    breakpoints=open(breakpointfile).readlines()
    doubleline=0; #double line event
    prebp=[];
    preType="";
    outString="";
    newlines=[];
    i=0;
    for bp in breakpoints:
        bp=bp.strip().split(delimeter)
        if bp[2]=="INV" or "FUS" in bp[2]:
            doubleline+=1
            if doubleline==2:
                doubleline=0
                if bp[2]=="INV":
                    newl=[bp[2],bp[3],bp[0],bp[1]]
                else:
                    pks=[int(prebp[0]),int(bp[0]),int(bp[1]),int(prebp[1])]
                    pks.sort()
                    newl=[bp[2],bp[3]]+pks
                newlines[-1]=(delimeter).join(str(x) for x in newl)
                i+=1
            else:
                newl=[bp[2],bp[3],bp[0],bp[1]]
                newlines.append((delimeter).join(newl));
            
        else:
            newl=[bp[2],bp[3],bp[0],bp[1]]
            newlines.append((delimeter).join(newl));
            
        prebp=bp
    
    print str(len(breakpoints)) + " lines in original breakpoint file"
    print str(i)+" lines merged"
    print str(len(newlines)) +" lines after merge"
    outstr=""
    for l in newlines:
        print l
        outstr+=l+'\n'
    open(outputfile,'w').write(outstr)
    
if __name__=='__main__':
    if len(sys.argv)<3:
        print "Usage: "+sys.argv[0]+"breakpoints.txt output"
    
    eventFromBreakpointFile(sys.argv[1], sys.argv[2])
#         exit()
# # breakpointfile="data2/simulated_chr12_1.fa_breakpoints.txt"
# breakpointfile=sys.argv[1]
# breakpoints=open(breakpointfile).readlines()
# doubleline=0; #double line event
# prebp=[];
# preType="";
# outString="";
# newlines=[];
# i=0;
# for bp in breakpoints:
#     bp=bp.strip().split(delimeter)
#     if bp[2]=="INV" or "FUS" in bp[2]:
#         doubleline+=1
#         if doubleline==2:
#             doubleline=0
#             if bp[2]=="INV":
#                 newl=[bp[2],bp[3],bp[0],bp[1]]
#             else:
#                 pks=[int(prebp[0]),int(bp[0]),int(bp[1]),int(prebp[1])]
#                 pks.sort()
#                 newl=[bp[2],bp[3]]+pks
#             newlines[-1]=(delimeter).join(str(x) for x in newl)
#             i+=1
#         else:
#             newl=[bp[2],bp[3],bp[0],bp[1]]
#             newlines.append((delimeter).join(newl));
#         
#     else:
#         newl=[bp[2],bp[3],bp[0],bp[1]]
#         newlines.append((delimeter).join(newl));
#         
#     prebp=bp
# 
# print str(len(breakpoints)) + " lines in original breakpoint file"
# print str(i)+" lines merged"
# print str(len(newlines)) +" lines after merge"
# outstr=""
# for l in newlines:
#     print l
#     outstr+=l+'\n'
# open(sys.argv[2],'w').write(outstr)
        