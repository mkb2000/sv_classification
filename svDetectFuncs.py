'''
collection of functions for rule-based approach.

Created on Jan 8, 2014

@author: HELE
'''

blurbp1=15 #for difference of breakpoint loci between c1 and c2. 10pb
blurbp2=100 #for difference of breakpoint loci between two continuous lines. 150bp
blurbp3=15 #for novel insertion detection
delimeter="\t"
htmlStartDiv='<div>'
htmlEndDiv='</div>\n'

class SVtypes:
    error=0
    tandem=1
    deletion=2
    inversion=4
    interspersedInsertion=tandem+deletion
    translocation=8
    novelInsertion=16

def getResultString(result,content):
    string=""
    bpstr=""
    r=result[0]
    if r==-1: # deleted by merging two lines
        return
    if r==-2:
        return
        string="deleted for translocation"
    if r==SVtypes.tandem:
        string="tandem"
    if r==SVtypes.deletion:
        string="deletion"
    if r==SVtypes.inversion:
        string="inversion"
    if r==SVtypes.interspersedInsertion:
        string="interspersed insertion"
    if r==SVtypes.translocation:
        string="translocation"
    if r==SVtypes.novelInsertion:
        string="novel insertion"
    if r==SVtypes.error:
        string="error"
    for breakpoint in result[1]:
        bpstr+=delimeter+breakpoint
    return string+bpstr
    
def printResultWithLine(result,content):
    string=""
    bpstr=""
    r=result[0]
    if r==-1:
        return
    if r==SVtypes.tandem:
        string="tandem"
    if r==SVtypes.deletion:
        string="deletion"
    if r==SVtypes.inversion:
        string="inversion"
    if r==SVtypes.interspersedInsertion:
        string="interspersed insertion"
    if r==SVtypes.translocation:
        string="translocation"
    if r==SVtypes.novelInsertion:
        string="novel insertion"
    if r==SVtypes.error:
        string="error"
    for breakpoint in result[1]:
        bpstr+=delimeter+breakpoint
    print string+bpstr
#     print string+":"+content[l].split("\t")[0].split(":")[1]

def detect (prevLine,prevResult,line):
    if prevResult:
        prevResult=prevResult[0]
    result=-1
    line=line.strip().split(delimeter)
    c1_realign=line[0].split(':')[1]
    c1_realign_dir=line[1]
    c1_realign_consensus=line[2]
    c1_anchor=line[3].split(':')[1]
    c1_anchor_dir=line[4]
    c1_anchor_consensus=line[5]
    c1_long_support=line[6]
    c1_long_support_bases=line[7]
    c1_short_support=line[8]
    c1_short_support_bases=line[9]
    c1_short_support_max_len=line[10]
    c1_avg_realign_mapq=line[11]
    c2_realign=line[12].split(':')[1]
    c2_realign_dir=line[13]
    c2_realign_consensus=line[14]
    c2_anchor=line[15].split(':')[1]
    c2_anchor_dir=line[16]
    c2_anchor_consensus=line[17]
    c2_long_support=line[18]
    c2_long_support_bases=line[19]
    c2_short_support=line[20]
    c2_short_support_bases=line[21]
    c2_short_support_max_len=line[22]
    c2_avg_realign_mapq=line[23]
    otherInfo=line[24]
    
    if c1_realign_dir!=c1_anchor_dir and c2_realign_dir!=c2_anchor_dir:
#         not inversion
        
        #if +anchor_loci>-realign_loci and -anchor_loci<+realign_loci, is tandem
        #if +anchor_loci<-realign_loci and -anchor_loci>+realign_loci, is deletion
        if (int(c1_anchor)>int(c1_realign) and c1_anchor_dir=="+") and (int(c2_anchor)<int(c2_realign) and c2_anchor_dir=="-"):
            result=SVtypes.tandem
        elif (int(c2_anchor)>int(c2_realign) and c2_anchor_dir=="+") and (int(c1_anchor)<int(c1_realign) and c1_anchor_dir=="-"):
            result=SVtypes.tandem
        elif (int(c2_anchor)>int(c2_realign) and c2_anchor_dir=="-") and (int(c1_anchor)<int(c1_realign) and c1_anchor_dir=="+"):
            result=SVtypes.deletion
        elif (int(c1_anchor)>int(c1_realign) and c1_anchor_dir=="-") and (int(c2_anchor)<int(c2_realign) and c2_anchor_dir=="+"):
            result=SVtypes.deletion
        elif abs(int(c1_anchor)-int(c1_realign))<blurbp3:
            result=SVtypes.novelInsertion
        else:
            result=SVtypes.error
        
        #if two lines, interspersed insertion = del + tandem
        if prevLine:
            prevLine=prevLine.strip().split(delimeter)
            distance1=int(prevLine[0].split(':')[1])-int(c1_realign)
            distance2=int(prevLine[3].split(':')[1])-int(c1_anchor)
            if abs(distance1)<blurbp2 or abs(distance2)<blurbp2:#if this two events are nearby.
                if result+prevResult==SVtypes.interspersedInsertion:
                    result=SVtypes.interspersedInsertion
        
    elif c1_realign_dir==c1_anchor_dir and c2_realign_dir==c2_anchor_dir:
#         inversion
        result=SVtypes.inversion
    else:
        result=SVtypes.error

    startPoint=min(c1_realign,c1_anchor,c2_realign,c2_anchor)
    endPoint=max(c1_realign,c1_anchor,c2_realign,c2_anchor)
    if "nserted sequence" in otherInfo and abs(int(startPoint)-int(endPoint))<5:
        result=SVtypes.novelInsertion
    return [result,[startPoint,endPoint]];

def detectTransloc(i,resultList):
    tolerance=5 #tolerance for blunt.
    #find the mobile part. r:[interInserType,[p1,p2,p3,p4]]
    r=resultList[i]
    p1=int(r[1][0])
    p2=int(r[1][1])
    p3=int(r[1][2])
    p4=int(r[1][3])
    mobilePart=[p1,p2] if abs(int(p1)-int(p2))>tolerance else [p3,p4]
    # try to find if the mobile part is deleted. If so, it is translocation.
    ind=0
    for rr in resultList:
        if rr[0]==SVtypes.deletion:
            p1=int(rr[1][0])
            p2=int(rr[1][1])
            if abs(mobilePart[0]-p1)<tolerance and abs(mobilePart[1]-p2)<tolerance:
                #is translocation
                resultList[ind][0]=-2
                resultList[i][0]=SVtypes.translocation
        ind+=1
            
    
def realignLoci(line):
    return int(line.strip().split(delimeter)[0].split(":")[1])
def anchorLoci(line):
    return int(line.strip().split(delimeter)[3].split(":")[1])
def wrapDIV(str):
    return '<p>'+str+'</p>\n'

def wrapColor(str,kind):
    htmlRightWrapper='<div class="alert alert-success">'
    htmlWrongWrapper='<div class="alert alert-danger">'
    htmlOtherWrapper='<div class="alert alert-warning">'
    htmlOther2Wrapper='<div class="alert alert-info">'
    return str
    if kind=="right":
        return htmlRightWrapper+str+htmlEndDiv
    if kind=="wrong":
        return htmlWrongWrapper+str+htmlEndDiv
    if kind=="other":
        return htmlOtherWrapper+str+htmlEndDiv
    if kind=="other2":
        return htmlOther2Wrapper+str+htmlEndDiv
    
def writeComapreResultToHTML(myResult,stand,compareResult,recal,precision):
    myResult=myResult.strip().split("\n")
    stand=stand.strip().split("\n")
    checkroll=[]
    for x in compareResult:
        for y in x:
            checkroll.append(y)
    
    htmlHead=open("head.html").read();
    htmlFoot=open("foot.html").read();
    htmlBody=""
    if not recal==0 and not precision==0:
        htmlBody+="Recal:"+str(recal)+"; Precision:"+str(precision)
    htmlBody+="<tr><td>Detect result</td><td>Source in breakpoint file</td><td>line No. in breakpoint file</td></tr>";

    
    i=0;
    unusedbk=0;
    lastStandPoint=0;
    for cr in compareResult:
        htmlCol1="";
        htmlCol2="";
        htmlclass="";
        if cr==[]: #my result has extra line
            htmlCol1=wrapColor(wrapDIV(myResult[i]), "other2")
            htmlCol2=wrapColor(wrapDIV("None source in breakpoint file"), "other2")
            htmlclass="other2"
            htmlBody+="<tr class='"+htmlclass+"'><td>"+htmlCol1+"</td>"+"<td>"+htmlCol2+"</td><td>"+"-</td></tr>"
        elif cr[0]==-2:
            htmlCol2=wrapColor(wrapDIV(stand[lastStandPoint+1]), "other");
            htmlCol1=wrapColor(wrapDIV("Deletion is part of translocation"), "other");
            htmlclass="other"
            htmlBody+="<tr class='"+htmlclass+"'><td>"+htmlCol1+"</td><td>"+htmlCol2+"</td><td>"+str(lastStandPoint+1)+"</td></tr>"
            checkroll.append(lastStandPoint+1)
        else:
            jump=cr[0]-lastStandPoint
#             print jump
## This is not very good, may miss breakpoints in standard.
            if jump>1: #standard has extra line, which may be  missed by Socrates.
                for x in range(1,jump):
                    if lastStandPoint+x not in checkroll and (lastStandPoint+x) < len(stand):
                        checkroll.append(lastStandPoint+x)
                        htmlCol2=wrapColor(wrapDIV(stand[lastStandPoint+x]), "other");
                        htmlCol1=wrapColor(wrapDIV("Cannot be detected"), "other");
                        htmlclass="other"
                        htmlBody+="<tr class='"+htmlclass+"'><td>"+htmlCol1+"</td><td>"+htmlCol2+"</td><td>"+str(lastStandPoint+x)+"</td></tr>"
                        unusedbk+=1
#                         print "unused dk:",lastStandPoint+x
                lastStandPoint+=jump
                
            if len(cr)>0:
                col2Str=""
                col1Str=""
                col1Str=wrapDIV(myResult[i])
                for x in cr:
                    col2Str+=wrapDIV(stand[x])
                htmlCol1=wrapColor(col1Str, "right");
                htmlCol2=wrapColor(col2Str, "right");
                lastStandPoint=min(cr)
                htmlclass="right"
                htmlBody+="<tr class='"+htmlclass+"'><td>"+htmlCol1+"</td><td>"+htmlCol2+"</td><td>"+str(cr)+"</td></tr>"
        i+=1
    html=htmlHead+htmlBody+htmlFoot
    open("index.html",'w').write(html)
#     print "unused bk amount:", unusedbk
    
# def writeComapreResultToHTML(myResult,stand,compareResult):
#     myResult=myResult.strip().split("\n")
#     stand=stand.strip().split("\n")
#     checkroll=[]
#     for x in compareResult:
#         for y in x:
#             checkroll.append(y)
#     
#     htmlHead=open("head.html").read();
#     htmlFoot=open("foot.html").read();
#     htmlBody="";
#     htmlCol1="";
#     htmlCol2="";
#     
#     i=0;
#     lastStandPoint=0;
#     for cr in compareResult:
#         
#         if cr==[]: #my result has extra line
# #             htmlBody+=wrapColor(wrapDIV(myResult[i]), "other2")
#             htmlCol1+=wrapColor(wrapDIV(myResult[i]), "other2")
#             htmlCol2+=wrapColor(wrapDIV("None source in breakpoint file"), "other2")
#         else:
#             jump=cr[0]-lastStandPoint
#             if jump>1: #stand has extra line, which may be a del in translo, or missed by Socrates.
#                 for x in range(1,jump):
#                     if lastStandPoint+x not in checkroll:
# #                         htmlBody+=wrapColor(wrapDIV(stand[lastStandPoint+x]), "other");
#                         htmlCol2+=wrapColor(wrapDIV(stand[lastStandPoint+x]), "other");
#                         htmlCol1+=wrapColor(wrapDIV("Cannot be detected"), "other");
#                 lastStandPoint+=jump
#             if len(cr)>0:
#                 col2Str=""
#                 col1Str=""
#                 col1Str=wrapDIV(myResult[i])
#                 for x in cr:
#                     col2Str+=wrapDIV(stand[x])
# #                 htmlBody+=wrapColor(lineStr, "right");
#                 htmlCol1+=wrapColor(col1Str, "right");
#                 htmlCol2+=wrapColor(col2Str, "right");
#                 lastStandPoint=cr[-1]
#         i+=1
#     
#     htmlrow='<div class="col-md-5">'
#     htmlrowx='<div class="col-md-2">'
# #     html=htmlHead+htmlBody+htmlFoot
#     html=htmlHead+htmlrow+htmlCol1+htmlEndDiv+htmlrow+htmlCol2+htmlEndDiv+htmlrowx+htmlEndDiv+htmlFoot
#     open("index.html",'w').write(html)
    