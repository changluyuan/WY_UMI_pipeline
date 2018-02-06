import pysam,os
import sys
import gzip
import datetime,time
from Levenshtein import *
from collections import OrderedDict
if len(sys.argv) <> 4:
    print('Usage: '+sys.argv[0] + '  bamfile  outfile.bam  minmem')
    exit(-1)
bamfile = pysam.AlignmentFile(sys.argv[1],'rb')
outname = sys.argv[2][:len(sys.argv[2])-4]
outfile_tmp = pysam.AlignmentFile(outname+'.tmp.bam','wb',template=bamfile)
outfile_mk = pysam.AlignmentFile(sys.argv[2],'wb',template=bamfile)
def consensus_maker(grouped_reads_list, read_length):
    # The consensus maker uses a simple "majority rules" algorithm to qmake a consensus at each base position.  If no
    # nucleotide majority reaches above the minimum theshold (--cut_off), the position is considered undefined and an 'N'
    # is placed at that position in the read.'''
    nuc_identity_list = [0, 0, 0, 0, 0, 0]  # In the order of T, C, G, A, N, Total
    nuc_key_dict = {0: 'T', 1: 'C', 2: 'G', 3: 'A', 4: 'N'}
    consensus_read = ''
    cut_off = 0.7
    for i in xrange(read_length):  # Count the types of nucleotides at a position in a read. i is the nucleotide index
        # within a read in grouped_reads_list
        for j in xrange(len(grouped_reads_list)):  # Do this for every read that comprises a SMI group. j is the read
            # index within grouped_reads_list
            try:
                if grouped_reads_list[j][i] == 'T':
                    nuc_identity_list[0] += 1
                elif grouped_reads_list[j][i] == 'C':
                    nuc_identity_list[1] += 1
                elif grouped_reads_list[j][i] == 'G':
                    nuc_identity_list[2] += 1
                elif grouped_reads_list[j][i] == 'A':
                    nuc_identity_list[3] += 1
                elif grouped_reads_list[j][i] == 'N':
                    nuc_identity_list[4] += 1
                else:
                    nuc_identity_list[4] += 1
                nuc_identity_list[5] += 1
            except:
                break
        try:
            for j in [0, 1, 2, 3, 4]:
                if float(nuc_identity_list[j])/float(nuc_identity_list[5]) >= cut_off:
                    consensus_read += nuc_key_dict[j]
                    break
                elif j == 4:
                    consensus_read += 'N'
        except:
            consensus_read += 'N'
        nuc_identity_list = [0, 0, 0, 0, 0, 0]  # Reset for the next nucleotide position
    return consensus_read, len(grouped_reads_list)

def get_consensus_reads(read_dict):
### minmem is the minimal number of repilicate reads which will be used to error suppression
    minmem = int(sys.argv[3])
    polyN = 'N'*75
    for dict_tag in read_dict.keys():
        cigar_string_set = {}
        for cigar_string in read_dict[dict_tag][7].keys():
            cigar_string_set[cigar_string] = read_dict[dict_tag][7][cigar_string][0]
        max_cigar = max(cigar_string_set,key=cigar_string_set.get)
        print(cigar_string_set[max_cigar])
        if cigar_string_set[max_cigar] >= minmem:
            consensus, fam_size = consensus_maker(read_dict[dict_tag][7][max_cigar][2:],
                                                    len(read_dict[dict_tag][7][max_cigar][2]))
#            if (consensus.count("N")/float(len(consensus))) <= 0.9:  ##the value can be changed
            read = pysam.AlignedRead()
            read.qname = read_dict[dict_tag][0]
            read.flag = read_dict[dict_tag][1]
            read.query_sequence = consensus
            read.reference_id = read_dict[dict_tag][2]
            read.reference_start = read_dict[dict_tag][3]
            read.mapq = 60
            read.cigarstring = read_dict[dict_tag][7][max_cigar][1]
            read.next_reference_id = read_dict[dict_tag][4]
            read.next_reference_start = read_dict[dict_tag][5]
            read.template_length = read_dict[dict_tag][6]
            read.query_qualities = 'J' * len(consensus)
            outfile_tmp.write(read)
            outfile_mk.write(read)
        else:
            read = pysam.AlignedRead()
            read.qname = read_dict[dict_tag][0]
            read.flag = read_dict[dict_tag][1] + 1024
            read.query_sequence = polyN
            read.reference_id = read_dict[dict_tag][2]
            read.reference_start = read_dict[dict_tag][3]
            read.mapq = 0
            read.cigarstring = '*'
            read.next_reference_id = read_dict[dict_tag][4]
            read.next_reference_start = read_dict[dict_tag][5]
            read.template_length = read_dict[dict_tag][6]
            read.query_qualities = 'J' * 75
            outfile_tmp.write(read)
            outfile_mk.write(read)
pattern=['ACTGCATA','GAGCCTTA','TATCCTCT','TCTCTCCG','AAGGCTAT','AAGGAGTA','CTAAGCCT','TTATGCGA','CGTCTAAT','CTCTCTAT','TCGACTAG','TTCTAGCT','CTATTAAG','GTAAGGAG','CCTAGAGT']
umi_start_reverse={}
umi_start={}
pos_to_name={}
name_to_position = {}
write_name = {}
read_number_count=0
wrong_umi = 0
umi_1 = 0
umi_2 = 0
for read in bamfile.fetch():
    #quality = 'E' * len(read.query_sequence)
    qname = read.query_name
    mindist = 100
    if read_number_count % 1000000 == 0:
        sys.stderr.write("Reads processed:" + str(read_number_count) + "\t")
        sys.stderr.write(time.strftime("%Y-%m-%d %H:%M:%S") + "\n")
    if read.is_read2 and not read.is_unmapped and read.mapping_quality >=55:
        umi = qname[len(qname)-8:]
        for key in pattern:
            dist=distance(umi,key)
            if dist < mindist:
                mindist = dist
                minpattern=key
        if umi not in pattern:
            wrong_umi += 1
        if mindist == 2:
            umi_2 += 1
        if mindist == 1:
            umi_1 +=1
            fixumi = minpattern
        elif mindist == 0:
            fixumi = minpattern
        elif mindist > 1:
            read_number_count+=1
            continue
        chrom = read.reference_name
        if read.is_proper_pair:
            flag = '1'
        else:
            flag = '0'
        if read.is_reverse:
            start = read.reference_end
            if (fixumi,chrom,start,'reverse') not in pos_to_name:
                pos_to_name[(fixumi,chrom,start,'reverse')] = [(qname,flag)] 
            else:
                pos_to_name[(fixumi,chrom,start,'reverse')].append((qname,flag))
        else:
            start = read.reference_start
            if (fixumi,chrom,start,'forward') not in pos_to_name:
                pos_to_name[(fixumi,chrom,start,'forward')] = [(qname,flag)]
            else:
                pos_to_name[(fixumi,chrom,start,'forward')].append((qname,flag))
    read_number_count+=1
#print(wrong_umi,umi_1,umi_2)
choose = 0
for (pos,names) in pos_to_name.iteritems():
    for name in names:
        name_to_position[name[0]] = pos
        if name[1] == '1' and not choose:
            write_name[name[0]] = 1
            choose = 1
    if not choose:
        write_name[names[0][0]] = 1
    choose = 0
bamfile.close()
bamfile_1 = pysam.AlignmentFile(sys.argv[1],'rb') 
read_dict_R1 = {}
read_dict_R2 = {}
read_number_count = 0
#until_eof=True can read 77/147 unmappped reads
for read in bamfile_1.fetch(until_eof=True):
    qname = read.query_name
    if read_number_count % 1000000 == 0:
        sys.stderr.write("Reads processed:" + str(read_number_count) + "\t")
        sys.stderr.write(time.strftime("%Y-%m-%d %H:%M:%S") + "\n")
    if read.is_read1:
        if qname not in name_to_position:
            if not read.is_unmapped:
                read.flag += 1024
                outfile_mk.write(read)
            else:
                outfile_mk.write(read)
            continue
        unique_pos = name_to_position[qname] 
        if unique_pos not in read_dict_R1:
            read_dict_R1[unique_pos] = [qname,read.flag,read.reference_id,read.reference_start,read.next_reference_id,read.next_reference_start,read.template_length,
                                               {str(read.cigarstring):[0,read.cigarstring]}]
        if not read.is_unmapped:
            if str(read.cigarstring) not in read_dict_R1[unique_pos][7]:
                read_dict_R1[unique_pos][7][str(read.cigarstring)] = [0, read.cigarstring]
            read_dict_R1[unique_pos][7][str(read.cigarstring)].append(read.query_sequence)
            read_dict_R1[unique_pos][7][str(read.cigarstring)][0] += 1
        if qname in write_name:
            read_dict_R1[unique_pos][0] = qname
            read_dict_R1[unique_pos][1] = read.flag
            read_dict_R1[unique_pos][2] = read.reference_id
            read_dict_R1[unique_pos][3] = read.reference_start
            read_dict_R1[unique_pos][4] = read.next_reference_id
            read_dict_R1[unique_pos][5] = read.next_reference_start
            read_dict_R1[unique_pos][6] = read.template_length
    if read.is_read2:
        #chrom = read.reference_name
        if qname not in name_to_position:
            if not read.is_unmapped:
                read.flag += 1024
                outfile_mk.write(read)
            else:
                outfile_mk.write(read)
            continue
        unique_pos = name_to_position[qname]
        if unique_pos not in read_dict_R2:
            read_dict_R2[unique_pos] = [qname,read.flag,read.reference_id,read.reference_start,read.next_reference_id,read.next_reference_start,read.template_length,
                                                   {str(read.cigarstring):[0,read.cigarstring]}]
        if str(read.cigarstring) not in read_dict_R2[unique_pos][7]:
            read_dict_R2[unique_pos][7][str(read.cigarstring)] = [0, read.cigarstring]
        read_dict_R2[unique_pos][7][str(read.cigarstring)].append(read.query_sequence)
        read_dict_R2[unique_pos][7][str(read.cigarstring)][0] += 1
        if qname in write_name:
            read_dict_R2[unique_pos][0] = qname
            read_dict_R2[unique_pos][1] = read.flag
            read_dict_R2[unique_pos][2] = read.reference_id
            read_dict_R2[unique_pos][3] = read.reference_start
            read_dict_R2[unique_pos][4] = read.next_reference_id
            read_dict_R2[unique_pos][5] = read.next_reference_start
            read_dict_R2[unique_pos][6] = read.template_length
    if qname not in write_name:
        if not read.is_unmapped:
            read.flag += 1024
            outfile_mk.write(read)
        else:
            outfile_mk.write(read)
    read_number_count+=1
get_consensus_reads(read_dict_R1)
get_consensus_reads(read_dict_R2)  
bamfile_1.close()
outfile_tmp.close()
os.system('samtools sort -n {} > {}'.format(outname+'.tmp.bam',outname+'.sorted.bam'))
os.system('bamToFastq -i {} -fq {} -fq2 {}'.format(outname+'.sorted.bam',outname+'.R1.fastq',outname+'.R2.fastq'))
os.system('gzip {} {}'.format(outname+'.R1.fastq',outname+'.R2.fastq'))
os.system('rm {} {}'.format(outname+'.tmp.bam',outname+'.sorted.bam'))

