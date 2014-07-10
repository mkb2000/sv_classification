'''
obsolete

Created on Feb 5, 2014

@author: HELE
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
blurbp1=15 #for difference of breakpoint loci between c1 and c2. 10pb
blurbp2=100 #for difference of breakpoint loci between two continuous lines. 150bp
blurbp3=15 #for novel insertion detection
delimeter="\t"
delimeter_loci=":"
deout=',' #delimeter_output
feaNum=9#3 #num of features extracted from C1, that 18 for a paired scorates-line.
maxLen=feaNum*4 #two line event length

def loci_compar(l1,l2):
    return '1' if c1_realign_loci<c1_anchor_loci else '0'
def dir_transform(d):
    return '1' if l[1].strip()=='+' else '0'
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
    datastring=''
    #header
    datastring+="@relation 'sv'\n"
    num_of_attri=len(rawdata[0].strip().split(deout))
    for i in range(0,num_of_attri-1):
        datastring+="@attribute F"+str(i)+" real\n"
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
def sparse_long_line_after_paired_data(lines,loci):
    toAdd=[]
    toRemove=[]
    toAdd_loci=[]
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
            toRemove.append(i)
    toRemove.sort(cmp=None, key=None, reverse=True)
    for tr in toRemove:
        del lines[tr]
        del loci[tr]
    return [lines+toAdd,loci+toAdd_loci]

def sparse_long_line_after_unpaired_data(lines,loci):
    toAdd=[]
    toRemove=[]
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
            toRemove.append(i)
    toRemove.sort(cmp=None, key=None, reverse=True)
    for tr in toRemove:
        del lines[tr]
        del loci[tr]
    return [lines+toAdd,loci+toAdd_loci]

def check_max_length_of_line(lines):
    maxl=0
    for l in lines:
        ll=len(l.split(','))
        if ll>maxl:
            maxl=ll
    return maxl
    
            
        
    
    
if __name__=='__main__':
    if len(sys.argv)<3:
        print "Usage: "+sys.argv[0]+"Socrates_paired_result.txt outputDir [break file] [unpaired file]"
        exit()
    
    sc_paired_file=sys.argv[1]
    outDIR=sys.argv[2]
    sc_unpaired_file=sys.argv[4] if len(sys.argv)>4 else get_unpaired_file_name(sc_paired_file)
    bk_file=sys.argv[3] if len(sys.argv)>3 else get_bk_file_name(sc_paired_file)
        
    file_basename=get_basename(sc_paired_file)
    
#     filename="./data1/results_Socrates_paired_simulated_ecoli_5_simseq_s_long_sc_l25_q5_m5_i95.txt"
    fhandler=open(sc_paired_file)
    filelines=fhandler.readlines()
    content=filelines[1:]
    loci=[]
    newlines=[]
    maxLen=0
    ## Paired content
    for l in content:
        l=l.strip().split(delimeter)
        c1_realign_loci=l[0].strip().split(delimeter_loci)[1]
        c1_anchor_loci=l[3].strip().split(delimeter_loci)[1]
        c2_realign_loci=l[12].strip().split(delimeter_loci)[1]
        c2_anchor_loci=l[15].strip().split(delimeter_loci)[1]
        locations=[c1_realign_loci,c1_anchor_loci,c2_realign_loci,c2_anchor_loci]
    #     loci.append()
        c1_loci_com=loci_compar(c1_realign_loci, c1_anchor_loci)
        c2_loci_com=loci_compar(c2_realign_loci, c2_anchor_loci)
        
        c1_realign_dir=dir_transform(l[1].strip())
        c1_anchor_dir=dir_transform(l[4].strip())
        c2_realign_dir=dir_transform(l[13].strip())
        c2_anchor_dir=dir_transform(l[16].strip())
        
        #new line header:
        #c1_loci_compare c1_realign_dir c1_anchor_dir C1_long_support    C1_long_support_bases    C1_short_support    C1_short_support_bases    C1_short_support_max_len    C1_avg_realign_mapq c2_loci_com c2_realign_dir c2_anchor_dir  C2_long_support    C2_long_support_bases    C2_short_support    C2_short_support_bases    C2_short_support_max_len    C2_avg_realign_mapq
        newline=build_new_line(c1_loci_com,c1_realign_dir,c1_anchor_dir,l[6],l[7],l[8],l[9],l[10],l[11] ,c2_loci_com,c2_realign_dir,c2_anchor_dir,l[18],l[19],l[20],l[21],l[22],l[23])
#         newline=build_new_line(c1_loci_com,c1_realign_dir,c1_anchor_dir,c2_loci_com,c2_realign_dir,c2_anchor_dir)
        
        #Concatenate two lines. If any locus in this line exist in any previous line, concatenate these two.
        i=locus_in_pre_line(loci, locations)
        if i!=-1:
            newlines[i]=newlines[i]+deout+newline
            if len(newlines[i].strip().split(deout))>maxLen:
                maxLen=len(newlines[i].strip().split(deout))
            loci[i]=loci[i]+locations
        else:
            newlines.append(newline)
            loci.append(locations)
    
    [newlines,loci]=sparse_long_line_after_paired_data(newlines,loci)
#     print check_max_length_of_line(newlines)
    
    
    
#    Unpaired content
    fpcontent=open(sc_unpaired_file).readlines()
    for l in fpcontent:
        l=l.strip().split(delimeter)
        c1_realign_loci=l[0].strip().split(delimeter_loci)[1]
        c1_anchor_loci=l[3].strip().split(delimeter_loci)[1]
        locations=[c1_realign_loci,c1_anchor_loci]
        c1_loci_com=loci_compar(c1_realign_loci, c1_anchor_loci)
            
        c1_realign_dir=dir_transform(l[1].strip())
        c1_anchor_dir=dir_transform(l[4].strip())
            
        #new line header:
        #c1_loci_compare c1_realign_dir c1_anchor_dir C1_long_support    C1_long_support_bases    C1_short_support    C1_short_support_bases    C1_short_support_max_len    C1_avg_realign_mapq 
        newline=build_new_line(c1_loci_com,c1_realign_dir,c1_anchor_dir,l[6],l[7],l[8],l[9],l[10],l[11])
#         newline=build_new_line(c1_loci_com,c1_realign_dir,c1_anchor_dir)
            
        #Concatenate two lines. If any locus in this line exist in any previous line, concatenate these two.
        i=locus_in_pre_line(loci, locations)
        if i!=-1:
            newlines[i]=newlines[i]+deout+newline
            if len(newlines[i].strip().split(deout))>maxLen:
                maxLen=len(newlines[i].strip().split(deout))
            loci[i]=loci[i]+locations
        else:
            newlines.append(newline)
            loci.append(locations)
    
    [newlines,loci]=sparse_long_line_after_unpaired_data(newlines,loci)
    print "length of each line:",check_max_length_of_line(newlines)
    maxLen=feaNum*6
    make_every_line_same_length(newlines, maxLen)
    
    #add related event to each line
#     bkfile="./data1/mergedbk_simulated_ecoli_5.fa_breakpoints.txt"
    bkcontent=eventsFromBreakpointFile(bk_file)
    for bk in bkcontent:
        found=False
        bk=bk.strip()
        if bk:
            bkk=bk.split(delimeter)
            bks=bkk[2:]
            sv_class=bkk[0]
            if sv_class=="FUS+":
                sv_class="FUS"
                print "chang fus+ to fus"
            #give the sv type to the line which contains all the locations of an event.
            for i in range(0,len(loci)):
                numBKinLine=0;
                for bkpoint in bks:
                    if is_element_in_array_with_tolerance(bkpoint, loci[i]):
                        numBKinLine+=1
                if numBKinLine==len(bks):
                    found=True
    #             if is_element_in_array_with_tolerance(p1, loci[i]) or is_element_in_array_with_tolerance(p2, loci[i]):
                    if len(newlines[i].split(deout))==maxLen:
                        newlines[i]=newlines[i]+','+sv_class
                    
                
            #impossible to find a line contains all the locations of an event. So give this sv type to any line contains more than one of the locations of an event.
            #this line has not been assigned a sv type.
            if not found:
                for i in range(0,len(loci)):
                    numBKinLine=0;
                    for bkpoint in bks:
                        if is_element_in_array_with_tolerance(bkpoint, loci[i]):
                            numBKinLine+=1
                    if numBKinLine>0:
                        found=True
                        if len(newlines[i].split(deout))==maxLen:
                            newlines[i]=newlines[i]+','+sv_class
            
            if not found:
                print "this event cannot be found in the output file of Socrates",bk
            #TODO: multi sv class may be located to one line
                
    make_every_line_same_length(newlines, maxLen+1, 'WOFS') #wrong of socrates output
    
    for l in newlines:
        print l
    
    outf=outDIR+'/'+file_basename+".arff"
    outfh=open(outf,"w")
               
    outfh.write(data_for_weka(newlines))     
    
#     lof=open(outDIR+'/'+file_basename+"_loci.txt","w")
#     for ll in loci:
#         lof.write(str(ll)+"\n");



