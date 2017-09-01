#!/usr/bin/python
# Define dependent clases
def write_line(linechrom,linestart,lineend,linetype,linelen,linescore,lineprog) :
    global inscnt,delcnt,dupcnt,invcnt,etccnt
    tempstart = abs(linestart)
    tempend = abs(lineend)
    if tempstart<tempend  : 
        templen = tempend-tempstart
    else :
        templen = tempstart-tempend
    if templen != linelen :
        linelen = templen
    outextract='%s\t%i\t%i\t%s\t%i\t%s\t%s\n' % \
                    (linechrom,linestart,lineend,linetype,linelen,linescore,lineprog)
    if tempstart > tempend :
        linestart = tempend
        lineend = tempstart
    else :
        linestart = tempstart
        lineend = tempend
    outbed='%s\t%i\t%i\n' % (linechrom,linestart,lineend)
    if linetype == 'INS' :
        inscnt += 1
        ins_extract.write(outextract)
        ins_bed.write(outbed)
    elif linetype == 'DEL' :
        delcnt += 1
        del_extract.write(outextract)
        del_bed.write(outbed)
    elif linetype == 'DUP' :
        dupcnt += 1
        dup_extract.write(outextract)
        dup_bed.write(outbed)
    elif linetype == 'INV' :
        invcnt += 1
        inv_extract.write(outextract)
        inv_bed.write(outbed)
    else :   
        etccnt += 1
        etc_extract.write(outextract)
    return;

# Get args passed
from sys import argv
if len(argv) > 1 :
    script, inname = argv
    prefix = str(inname).split('.',1)[0] 
else :
    print ("No input file specified")
    exit(2)

# Print Start
print ("Analyzing File: %s" % inname)

input_file = open(inname, 'r')
ins_extract = open(prefix+'_ins_extract.txt','w')
del_extract = open(prefix+'_del_extract.txt','w')
dup_extract = open(prefix+'_dup_extract.txt','w')
inv_extract = open(prefix+'_inv_extract.txt','w')
etc_extract = open(prefix+'_etc_extract.txt','w')
ins_bed = open(prefix+'_ins_extract.bed','w')
del_bed = open(prefix+'_del_extract.bed','w')
dup_bed = open(prefix+'_dup_extract.bed','w')
inv_bed = open(prefix+'_inv_extract.bed','w')

# Initialize variables
inputcnt=outputcnt = 0
inscnt=delcnt=dupcnt=invcnt=etccnt = 0
filetype = None
fileprog = None

# Writing header lines
hdrline="Chrom1\tPOS1\tPOS2\tType\tLen\tScore\tProg\n"
ins_extract.write(hdrline)
del_extract.write(hdrline)
dup_extract.write(hdrline)
inv_extract.write(hdrline)
etc_extract.write(hdrline)

# Process the input file
import csv
csv.field_size_limit(500 * 1024 * 1024)

for line in csv.reader(input_file,delimiter='\t',quoting=csv.QUOTE_NONE) :
    inputcnt += 1
    if (inputcnt%1000==0) :
        print "...read records", inputcnt
    if len(line) > 0 :
        if inputcnt == 1 and line[0][:1] != '#' :
            ans = raw_input("Header record not present, enter program (cnvnator/?): ")
            if ans == 'cnvnator' :
                filetype = 'cnvnator'
                fileprog = 'cnvnator'
                print 'Cnvnator file type selected'
            else :
                print "Invalid program entered: ",ans
                exit(4)
        elif line[0][:1] == '#' :
            if line[0][:18] == '##fileformat=VCFv4' :
                filetype = 'vcf'
                print 'VCF header detected'
            if line[0][:16] == '##source=pindel' :
                fileprog = 'pindel'
                print 'Pindel header detected'
            if line[0][:8] == '#CHROM_A' :
                filetype = 'lumpy'
                fileprog = 'lumpy'
                print 'Lumpy heaader detected'
            if line[0][:21]=='#Command: breakdancer' :
                filetype = 'breakdancer'
                fileprog = 'breakdancer'
                print 'Breakdancer header detected'
            if filetype == 'breakdancer' and line[0] == '#Chr1':
                if len(line) > 7 :
                    print "Detected valid header: ",filetype
                else :
                    print "Invalid header record: ",line
                    exit(3)
        else :
            if filetype == 'breakdancer' :
                outputcnt +=1
                write_line(line[0],int(line[1]),int(line[4]),line[6],int(line[7]),int(line[8]),fileprog)
            elif filetype == 'cnvnator' :
                outputcnt +=1
                a = line[1].split(':')
                b = a[1].split('-')
                write_line(a[0],int(b[0]),int(b[1]),str(line[0][:3]).upper(),float(line[2]),'NA',fileprog)
            elif filetype == 'lumpy' :
                outputcnt +=1
                svlen = int(line[2])-int(line[1])
                write_line(line[0],int(line[1]),int(line[2]),str(line[10]),svlen,'NA',fileprog)
            elif filetype == 'vcf' :
                outputcnt +=1
                svtype = 'NUL'
                svend = 0
                svlen = 0
                a = line[7].split(';')
                d = {}
                for i in a :
                    if "=" in i :
                        b = i.split('=')
                        d[str(b[0])]=b[1]
                if "END" in d :
                    svend = int(d.get("END"))
                elif "SVEND" in d :
                    svend = int(d.get("SVEND"))
                if svend == 0 :
                    print "END not found"
                if "SVTYPE" in d :
                    svtype = d.get("SVTYPE")[:3]
                else :
                    print "SVTYPE not found"
                if "SVLEN" in d :
                    svlen = int(d.get("SVLEN"))
                else :
                    print "SVLEN not found"
                if "SVMETHOD" in d :
                    if 'DELLY' in d.get("SVMETHOD") :
                        fileprog = 'delly'    
                write_line(line[0],int(line[1]),svend,svtype,svlen,'NA',fileprog)
            else :
                print "Un-recognized file format"
                print line
                ans = raw_input("Press enter to continue (y/n): ")
                if ans.strip() == 'n' :
                    exit(99)


# Done, close files
print ("Processed: %s, input: %i output: %i ins: %i del: %i dup: %i inv: %i etc: %i " % \
       (inname,inputcnt,outputcnt,inscnt,delcnt,dupcnt,invcnt,etccnt))
ins_bed.close()
del_bed.close()
dup_bed.close()
inv_bed.close()
ins_extract.close()
del_extract.close()
dup_extract.close()
inv_extract.close()
etc_extract.close()

# Cleanup unused files and sort the BED files
import subprocess
print ("Cleaning up")
if inscnt > 0 :
    print "... sorting ins"
    opt = r"sort -t $'\t' -k1,1 -k2,2g -o "+prefix+"_ins_sort.bed "+prefix+"_ins_extract.bed"
    subprocess.call(opt,shell=True)
else :
    opt = r"rm "+prefix+"_ins_extract.bed "+prefix+"_ins_extract.txt "
    subprocess.call(opt,shell=True)
if delcnt > 0 :
    print "... sorting del"
    opt = r"sort -t $'\t' -k1,1 -k2,2g -o "+prefix+"_del_sort.bed "+prefix+"_del_extract.bed"
    subprocess.call(opt,shell=True)
else :
    opt = r"rm "+prefix+"_del_extract.bed "+prefix+"_del_extract.txt "
    subprocess.call(opt,shell=True)
if dupcnt > 0 :
    print "... sorting dup"
    opt = r"sort -t $'\t' -k1,1 -k2,2g -o "+prefix+"_dup_sort.bed "+prefix+"_dup_extract.bed"
    subprocess.call(opt,shell=True)
else :
    opt = r"rm "+prefix+"_dup_extract.bed "+prefix+"_dup_extract.txt "
    subprocess.call(opt,shell=True)
if invcnt > 0 :
    print "... sorting inv"
    opt = r"sort -t $'\t' -k1,1 -k2,2g -o "+prefix+"_inv_sort.bed "+prefix+"_inv_extract.bed"
    subprocess.call(opt,shell=True)
else :
    opt = r"rm "+prefix+"_inv_extract.bed "+prefix+"_inv_extract.txt "
    subprocess.call(opt,shell=True)
if etccnt == 0 :
    opt = r"rm "+prefix+"_etc_extract.txt "
    subprocess.call(opt,shell=True)

