'''
extract features for ML from a Socrates output. Include paired and unpaired output of Socrates.

Three-line for TRANs, two-line for FUS and INV, one line fore DEL and INS, TAN


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
from ML_funcs import *
    
    
if __name__=='__main__':
    if len(sys.argv)<3:
        print "Usage: "+sys.argv[0]+" Socrates_paired_result.txt outputDir [break file] [unpaired file]"
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
    line_num_count=[]
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
            line_num_count[i]+=1
        else:
            newlines.append(newline)
            loci.append(locations)
            line_num_count.append(1)
    
    [newlines,loci,line_num_count]=sparse_long_line_after_paired_data(newlines,loci,line_num_count)
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
            line_num_count[i]+=1
        else:
            newlines.append(newline)
            loci.append(locations)
            line_num_count.append(1)
    
    [newlines,loci,line_num_count]=sparse_long_line_after_unpaired_data(newlines,loci,line_num_count)
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
    
    #add manual classifier result to the front
#     manual_result=manul_detect(sc_paired_file, bk_file);
#     for i in range(0,len(loci)):
#         locations=loci[i]
#         found=0
#         for l in manual_result:
#             sv_type=l[0]
#             locs=l[1]
#             if sv_type<0:
#                 continue
#             match_num=0
#             for loc in locs:
#                 if is_element_in_array_with_tolerance(loc, locations):
#                     match_num+=1
#             if match_num==len(locs):
#                 found=1;
#                 break
#         if found==1:
#             newlines[i]=str(sv_type)+","+newlines[i]
#         else:
#             newlines[i]="?"+","+newlines[i]

            #add line_num_count to the front as a feature
    for i in range(0,len(loci)):
        newlines[i]=str(line_num_count[i])+","+newlines[i]
    
    for l in newlines:
        print l
    
    outf=outDIR+'/'+file_basename+".arff"
    outfh=open(outf,"w")
               
    outfh.write(data_for_weka(newlines))     
    
#     lof=open(outDIR+'/'+file_basename+"_loci.txt","w")
#     for ll in loci:
#         lof.write(str(ll));
    
    




