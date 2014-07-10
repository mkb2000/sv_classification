'''
Rule-based SV classification.

Created on Dec 13, 2013

@author: HELE

Detect structural variation from output of Socrate.

Content of each line at each index:
C1_realign[0]
C1_realign_dir[1]
C1_realign_consensus[2]
C1_anchor[3]
C1_anchor_dir[4]
C1_anchor_consensus[5]
C1_long_support[6]
C1_long_support_bases[7]
C1_short_support[8]
C1_short_support_bases[9]
C1_short_support_max_len[10]
C1_avg_realign_mapq[11]
C2_realign[12]
C2_realign_dir[13]
C2_realign_consensus[14]
C2_anchor[15]
C2_anchor_dir[16]
C2_anchor_consensus[17]
C2_long_support[18]
C2_long_support_bases[19]
C2_short_support[20]
C2_short_support_bases[21]
C2_short_support_max_len[22]
C2_avg_realign_mapq[23]

'''
import sys
from svDetectFuncs import *
from SVResultCompare import SVResultComapreClass
from eventsFromBreakpointFile import eventFromBreakpointFile

if len(sys.argv)<3:
        print "Usage: "+sys.argv[0]+" Socrates_paired_result.txt breakpoints.txt"
        exit()


resultList=[]

# inputfile="data1/results_Socrates_paired_simulated_ecoli_5_simseq_s_long_sc_l25_q5_m5_i95.txt"
# standardfile="data1/simulated_ecoli_5.fa_breakpoints.txt"
inputfile=sys.argv[1] #"simulated_data/set2/results_Socrates_paired_simulated_chr12_1_simseq_s_long_sc_l25_q5_m5_i95.txt"

standardfile=sys.argv[2] #"simulated_data/set2/simulated_chr12_1.fa_breakpoints.txt"

inputcontent=open(inputfile).readlines()

head=inputcontent[0].rstrip().split(delimeter)
content=inputcontent[1:]


prevl=None
prevr=None
i=0

# predict SV for each line of Socrates output.
for l in content:
    result=detect(prevl, prevr, l)
    resultList.append(result)
#     x=resultList[i][1][0]
    # interspersed insertion has two lines, combine the result of previous line
    if result[0]==SVtypes.interspersedInsertion:  #or (result==SVtypes.inversion and prevr==SVtypes.inversion):
        resultList[i-1][0]=-1 #two line => one event. delete the previous result
        breakpoints=[resultList[i-1][1][0],resultList[i][1][0],resultList[i-1][1][1],resultList[i][1][1]]
        breakpoints.sort();
        resultList[i][1]=breakpoints
    # inversion has two lines, combine the result of previous line.
    if result[0]==SVtypes.inversion and prevl!=None and (abs(realignLoci(l)-realignLoci(prevl))<blurbp1 or abs(realignLoci(l)-anchorLoci(prevl)<blurbp1)):
        if resultList[i-1][0]==SVtypes.inversion:
            resultList[i-1][0]=-1
            resultList[i][1][0]=min(resultList[i-1][1][0],resultList[i][1][0])
            resultList[i][1][1]=max(resultList[i-1][1][1],resultList[i][1][1])
    
    i+=1
    prevl=l #previous line
    prevr=result #previous result
#     print "At line:"+str(i)
#     break

## after assay of each line, search for translocation=interspersed insertion + deletion.
i=0
for r in resultList:
    if r[0]==SVtypes.interspersedInsertion:
        detectTransloc(i,resultList)
    i+=1
    
#list to string. merge double line events
resultsString=""
for r in resultList:
    tempStr=getResultString(r,content)
    if tempStr:
        resultsString+=tempStr+"\n"

print "result str:\n",resultsString,"\n"


##compare result with standard.
newbkf=standardfile+"new.txt"
eventFromBreakpointFile(standardfile, newbkf)
compare=SVResultComapreClass(resultsString,standardfile,inputfile,newbkf)
compare.compare()
print compare.comparedResult

#check index.html within this folder for comparison result. 
standardStr=open(standardfile).read()
writeComapreResultToHTML(resultsString, standardStr, compare.comparedResult,compare.recall(),compare.precision())
print compare.recall(),compare.precision()