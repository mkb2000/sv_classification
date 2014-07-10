'''
Collection of functions for SV classification.

Created on Apr 16, 2014

@author: HELE
'''
from svDetectFuncs import *
import sys
blurbp1=15 #for difference of breakpoint loci between c1 and c2. 10pb
blurbp2=100 #for difference of breakpoint loci between two continuous lines. 150bp
blurbp3=15 #for novel insertion detection
delimeter="\t"
delimeter_loci=":"
deout=',' #delimeter_output
feaNum=9#3 #num of features extracted from C1, that 18 for a paired scorates-line.
maxLen=feaNum*4 #two line event length


def loci_compar(l1,l2):
    return '1' if l1<l2 else '0'
def dir_transform(d):
    return '1' if d=='+' else '0'
def build_new_line(*args):
    if len(args)%feaNum!=0:
        print "Err: num of new line"
        return;
    newline=''
    for arg in args:
        newline+=arg
        newline+=deout
    return newline[:-1]
def locus_in_pre_line(prelines,locations):
    for locus in locations:
        for i in range(0,len(prelines)):
            if is_element_in_array_with_tolerance(locus, prelines[i]):
                return i
    return -1
def is_element_in_array_with_tolerance(element,array,tolerance=blurbp1):
    for e in array:
        if abs(int(element)-int(e))<tolerance:
            return True
    return False 
def make_every_line_same_length(lines,length,placeholder='?'):
    for i in range(0,len(lines)):
        ll=len(lines[i].strip().split(','))
        filler=placeholder+','
        diff=length-ll
        if diff>0:
            filler=filler*diff
            lines[i]=lines[i]+deout+filler[:-1]
        elif diff<0:
            print "error in finding the max len:",lines[i]
def data_for_weka(rawdata):
    #TODO: feature string goes here
#     featureList=["Manual_result","c1_loci_compare","c1_realign_dir","c1_anchor_dir", "C1_long_support","C1_long_support_bases","C1_short_support","C1_short_support_bases","C1_short_support_max_len","C1_avg_realign_mapq", "c2_loci_com","c2_realign_dir", "c2_anchor_dir",  "C2_long_support","C2_long_support_bases","C2_short_support","C2_short_support_bases","C2_short_support_max_len","C2_avg_realign_mapq"];
    featureList=["num_lines","c1_loci_compare","c1_realign_dir","c1_anchor_dir", "C1_long_support","C1_long_support_bases","C1_short_support","C1_short_support_bases","C1_short_support_max_len","C1_avg_realign_mapq", "c2_loci_com","c2_realign_dir", "c2_anchor_dir",  "C2_long_support","C2_long_support_bases","C2_short_support","C2_short_support_bases","C2_short_support_max_len","C2_avg_realign_mapq"];
#     featureList=["Manual_result","loci_compare","realign_dir","anchor_dir", "long_support","long_support_bases","short_support","short_support_bases","short_support_max_len","avg_realign_mapq", "loci_com","c2_realign_dir", "c2_anchor_dir",  "C2_long_support","C2_long_support_bases","C2_short_support","C2_short_support_bases","C2_short_support_max_len","C2_avg_realign_mapq"];
    
    datastring=''
    #header
    datastring+="@relation 'sv'\n"
    num_of_attri=len(rawdata[0].strip().split(deout))
    for i in range(0,num_of_attri-1):
        if i==0:
            datastring+="@attribute "+featureList[i]+" real\n"
        else:
            ln=(i-1)/(feaNum*2);
            li=(i-1)%(feaNum*2);
            
            datastring+="@attribute "+"l"+str(ln)+"_"+featureList[li+1]+" real\n"
           
#         datastring+="@attribute F"+str(i)+" real\n"
    datastring+="@attribute svtype {INV,INS,DEL,FUS,TAN,TRAN,WOFS}\n"
    datastring+="@data\n"
    for l in rawdata:
        datastring+=l+'\n'
    return datastring

#change the position of original breakpoint file:svtype length breakpoint1 breakpoint2 [breakpoint3 breakpoint4]
#merge two line events:INV FUS
def eventsFromBreakpointFile(bk_file):
    breakpointfile=bk_file#"data1/simulated_ecoli_5.fa_breakpoints.txt"
    breakpoints=open(breakpointfile).readlines()
    doubleline=0; #double line event
    prebp=[];
    newlines=[];
    i=0;
    for bp in breakpoints:
        bp=bp.strip().split(delimeter)
        if bp[2]=="INV" or "FUS" in bp[2]:
            if "FUS" in bp[2]:
                bp[2]="FUS"
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
    
    #Translocation
    toRemoveList=[]
    for ii in xrange(0,len(newlines)):
        nl=newlines[ii]
        nls=nl.strip().split(delimeter) 
        bkpoints=nls[2:]
        if nls[0]=="FUS":
            for l in newlines:
                ls=l.strip().split(delimeter)
                if ls[0]=="DEL":
                    if is_element_in_array_with_tolerance(ls[2], bkpoints, 3) or is_element_in_array_with_tolerance(ls[3], bkpoints, 3):
                        nls[0]="TRAN"
                        newlines[ii]=(delimeter).join(nls)
                        toRemoveList.append(l)
                        print "removed del for TRAN:" + str(l)
    
    for tr in toRemoveList:
        newlines.remove(tr)
    
    print "bkfile:"+ str(len(breakpoints)) + " lines in original breakpoint file"
    print "bkfile:"+ str(i)+" lines merged"
    print "bkfile:"+ str(len(newlines)) +" lines after merge"
    return newlines

def get_unpaired_file_name(paired_file):
    paired_file=paired_file.split('/')
    pwd=''
    for i in xrange(0,len(paired_file)-1):
        pwd+=paired_file[i]+'/'
#     pwd=paired_file[0]+'/' if len(paired_file)>1 else ''
    paired_file=paired_file[-1].split('_')
    if paired_file[2]=='paired':
        paired_file[2]='unpaired'
    return pwd+'_'.join(paired_file)
def get_bk_file_name(paired_file):
    paired_file=paired_file.split('/')
    pwd=''
    for i in xrange(0,len(paired_file)-1):
        pwd+=paired_file[i]+'/'
#     pwd=paired_file[0]+'/' if len(paired_file)>1 else ''
    paired_file=paired_file[-1].split('_')
    bk_file=paired_file[3:-8]
    return pwd+'_'.join(bk_file)+".fa_breakpoints.txt"
def get_basename(paired_file):
    paired_file=paired_file.split('/')
    paired_file=paired_file[-1].split('_')
    bk_file=paired_file[3:-8]
    return '_'.join(bk_file)

'''
To make each line same length. The max length is 3 scorates-lines, which is translocation.
If any line is longer than 3 socrates-lines, every 3 socrates-lines become a new line for machine learning.
'''
def sparse_long_line_after_paired_data(lines,loci,line_num_count):
    toAdd=[]
    toRemove=[]
    toAdd_loci=[]
    toAdd_lc=[]
    for i in range(0,len(lines)):
        l=lines[i][:]
        l=l.split(',')
        length=len(l)
        if length>feaNum*6:
        #larger than 3-line event
            ls=[None]*(length/(feaNum*2)) #after cat paired file, a single line may contain multiple scorates line. this break down line and contain all scorates line. each scoratesline has c1 and c2
            for j in range(0,length/(feaNum*2)):
                ls[j]=l[feaNum*2*j:feaNum*2*(j+1)]
                ls[j]=','.join(ls[j])
            
            #every 3 scorates-line become a new line
            for j in range(0,len(ls)):
                for k in range(j+1,len(ls)):
                    for m in range(k+1,len(ls)):
                        toAdd.append(ls[j]+ls[k]+ls[m])
                        toAdd_loci.append(loci[i][:])#copy loci, prepare to add to loci[]
                        toAdd_lc.append(line_num_count[i]);
            toRemove.append(i)
    toRemove.sort(cmp=None, key=None, reverse=True)
    for tr in toRemove:
        del lines[tr]
        del loci[tr]
        del line_num_count[tr]
    return [lines+toAdd,loci+toAdd_loci,line_num_count+toAdd_lc]

def sparse_long_line_after_unpaired_data(lines,loci,line_num_count):
    toAdd=[]
    toRemove=[]
    toAdd_lc=[]
    toAdd_loci=[]
    for i in range(0,len(lines)):
        l=lines[i][:]
        l=l.split(',')
        length=len(l)
        if length>feaNum*6:
        #larger than 3-line event
            ls=[None]*(length/feaNum) #after cat paired file, a single line may contain multiple scorates line. this break down line and contain all scorates line. each scoratesline has c1 and c2
            for j in range(0,length/feaNum):
                ls[j]=l[feaNum*j:feaNum*(j+1)]
                ls[j]=','.join(ls[j])
            ls[4]=ls[4]+','+ls[5]
            del ls[5]
            ls[2]=ls[2]+','+ls[3]
            del ls[3]
            ls[0]=ls[0]+','+ls[1]
            del ls[1]
            
            #every 3 scorates-line become a new line
            for j in range(0,len(ls)):
                for k in range(j+1,len(ls)):
                    for m in range(k+1,len(ls)):
                        toAdd.append(ls[j]+deout+ls[k]+deout+ls[m])
                        toAdd_loci.append(loci[i][:])#copy loci, prepare to add to loci[]
                        toAdd_lc.append(line_num_count[i]);
            toRemove.append(i)
    toRemove.sort(cmp=None, key=None, reverse=True)
    for tr in toRemove:
        del lines[tr]
        del loci[tr]
        del line_num_count[tr]
    return [lines+toAdd,loci+toAdd_loci,line_num_count+toAdd_lc]

def check_max_length_of_line(lines):
    maxl=0
    for l in lines:
        ll=len(l.split(','))
        if ll>maxl:
            maxl=ll
    return maxl

def manul_detect(unpaired_file,bk_file):
#     if len(sys.argv)<3:
#         print "Usage: "+sys.argv[0]+" Socrates_paired_result.txt breakpoints.txt"
#         exit()


    resultList=[]
    
    # inputfile="data1/results_Socrates_paired_simulated_ecoli_5_simseq_s_long_sc_l25_q5_m5_i95.txt"
    # standardfile="data1/simulated_ecoli_5.fa_breakpoints.txt"
    inputfile=unpaired_file# sys.argv[1] #"simulated_data/set2/results_Socrates_paired_simulated_chr12_1_simseq_s_long_sc_l25_q5_m5_i95.txt"
    
    standardfile=bk_file #sys.argv[2] #"simulated_data/set2/simulated_chr12_1.fa_breakpoints.txt"
    
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
    return resultList
    #list to string. merge double line events
#     resultsString=""
#     for r in resultList:
#         tempStr=getResultString(r,content)
#         if tempStr:
#             resultsString+=tempStr+"\n"
     