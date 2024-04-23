import subprocess
import os
import datetime
import json

# Generate reverse strand of a DNA sequence
def reverse(string):
    length=len(string)
    string=string.upper()
    outstring = ''
    for pos in range(length):
        if string[length-pos-1] == 'A':
            outstring = outstring + 'T'
        elif string[length-pos-1] == 'T':
            outstring = outstring + 'A'
        elif string[length-pos-1] == 'G':
            outstring = outstring + 'C'
        elif string[length-pos-1] == 'C':
            outstring = outstring + 'G'
    return outstring

# Load genomic sequence into a dictionary
# Inputs a fasta of genomic sequence
# Returns a dictionary, key is contig name, value is sequence of the contig.
def load_sequence(file):
    sequence=open(file,'r').read()
    contigs=sequence.split('>')[1:]
    outputdic={}
    for contig in contigs:
        contig_split=contig.split('\n')
        name=contig_split[0]
        sequence=''.join(contig_split[1:])
        outputdic[name]=sequence
    return outputdic

# Reomve GIs that are part of other GIs
def remove_small_GI(dict):
    out_dic={}
    for key in dict.keys():
        keep=1
        for key2 in dict.keys():
            s1s=dict[key][0][1]
            s1e=dict[key][0][2]
            s2s=dict[key][1][1]
            s2e=dict[key][1][2]
            s1s2=dict[key2][0][1]
            s1e2=dict[key2][0][2]
            s2s2=dict[key2][1][1]
            s2e2=dict[key2][1][2]
            if key!=key2:
                if s1s>=s1s2 and s1e<=s1e2:
                    if s2s>=s2s2 and s2e<=s2e2:
                        keep=0
        if keep==1:
            out_dic[key]=dict[key]
    return out_dic

#Identifies integrases. Uses pre-calculated integrase library. Inputs a fasta file of prodigal predicted proteins.
#Returns a dictionary of predicted integrases.
#Each record is a list, containing contig number, start base, end base.
def identify_integrase(fasta,gff,threads,evalue,id,output,working):
    integrasedb=working+'db/recombinaseDB/'
    if 'integrase.phr' not in os.listdir(integrasedb):
        print('==========FIRST TIME USE, PREPARING INTEGRASE DATABASE==========')
        print(datetime.datetime.now())
        subprocess.run('makeblastdb -dbtype prot -in '+integrasedb+'GO_0004803_0008907.fasta -out '+integrasedb+'integrase',shell=True)
        print('==========INTEGRASE DATABASE PREPARED==========')
        print(datetime.datetime.now())
    if threads==1:
        subprocess.run('blastp -db '+integrasedb+'integrase -query '+fasta+' -out '+output+'/temp/blastresult.txt -evalue '+str(evalue)+" -outfmt \'6 qseqid pident\'",shell=True)
    else:
        subprocess.run('blastp -db '+integrasedb+'integrase -query '+fasta+' -out '+output+'/temp/blastresult.txt -evalue '+str(evalue)+" -outfmt \'6 qseqid pident\' -mt_mode 1 -num_threads "+str(threads),shell=True)
    blast=open(output+'temp/blastresult.txt','r').readlines()
    out = []
    for item in blast:
        list = item.split()
        if float(list[1]) >= float(id):
            if list[0] not in out:
                out.append(list[0])
    out_dic={}
    input=open(gff,'r').readlines()
    for line in input:
        line_list=line.split('\t')
        if len(line_list)>1:
            contig=line_list[0]
            start=int(line_list[3])
            end=int(line_list[4])
            info=line_list[-1]
            gene_number=info[3:].split(';')[0]
            if gene_number in out:
                out_dic[gene_number]=[]
                out_dic[gene_number].append(contig)
                out_dic[gene_number].append(start)
                out_dic[gene_number].append(end)
    return out_dic

#Create sequences to be searched for repeat sequences.
def create_seqs(integrase,size,size_spacer,genome_sequence,start,end):
    contig_length=len(genome_sequence[integrase[0]])
    if start==1:
        seq1a=''
        start_a=start
    elif start-size_spacer<=0:
        start_a=1
        seq1a=get_sequence(genome_sequence,integrase[0],1,start-1)
    else:
        start_a=start-size_spacer
        seq1a=get_sequence(genome_sequence,integrase[0],start_a,start-1)
    if end==contig_length:
        seq2a=''
        end_a=end
    elif end+size>=contig_length:
        end_a=contig_length
        seq2a=get_sequence(genome_sequence,integrase[0],end+1,contig_length)
    else:
        end_a=end+size
        seq2a = get_sequence(genome_sequence, integrase[0], end + 1,end_a)
    seq2a_rev=reverse(seq2a)
    if start==1:
        seq1b=''
        start_b=start
    elif start-size<=0:
        start_b=1
        seq1b=get_sequence(genome_sequence,integrase[0],1,start-1)
    else:
        start_b=start-size
        seq1b=get_sequence(genome_sequence,integrase[0],start_b,start-1)
    if end==contig_length:
        seq2b=''
        end_b=end
    elif end+size_spacer>=contig_length:
        end_b=contig_length
        seq2b=get_sequence(genome_sequence,integrase[0],end+1,contig_length)
    else:
        end_b=end+size_spacer
        seq2b = get_sequence(genome_sequence, integrase[0], end + 1,end_b)
    seq2b_rev=reverse(seq2b)
    return seq1a,seq2a,seq2a_rev,start_a,end_a,seq1b,seq2b,seq2b_rev,start_b,end_b

#Construt blast library for each genomic island, for pre-caculation when first time use.
def construct_blast_library(protein_dic, gi,output,working):
    with open(working+'db/blastdb/'+gi+'.fasta','w') as f:
        for protein in protein_dic.keys():
            f.write('>'+protein+'\n')
            f.write(protein_dic[protein]+'\n')
    subprocess.run('makeblastdb -dbtype prot -in '+working+'db/blastdb/' + gi + '.fasta -out '+working+'db/blastdb/' + gi,
                   shell=True)

# Calculate comparison scores with blast
# hit_score: The percentage of input GI genes found in reference GI
# completeness_score: The percentage of reference GI genes found in input GI
def giscore_blast(inputlist):
    id=float(inputlist[0])
    db_GI=inputlist[1]
    inputdic=inputlist[2]
    refdic=inputlist[3]
    evalue=float(inputlist[4])
    GI=inputlist[5]
    threads=int(inputlist[6])
    output=inputlist[7]
    term=inputlist[8]
    acc=inputlist[9]
    working=inputlist[10]
    gene_list_input=inputdic.keys()
    gene_list_ref=refdic.keys()
    total_input_length=len(gene_list_input)
    total_hit_score=0
    total_ref_length=len(gene_list_ref)
    total_completeness_score=0
    if not os.path.exists(output+'temp/input/'+db_GI):
        os.makedirs(output+'temp/input/'+db_GI)
    if not os.path.exists(output+'temp/ref/'+db_GI):
        os.makedirs(output+'temp/ref/'+db_GI)
    with open(output+'temp/proteins/'+GI,'r') as fgiscore:
        if len(fgiscore.read())>0:
            subprocess.run('blastp -db '+working+'db/blastdb/' + db_GI + ' -query '+output+'temp/proteins/'+ GI+
                           ' -out '+output+'temp/input/' +db_GI + '/temp.txt -evalue ' + str(evalue) +
                           ' -outfmt "6 qaccver pident"', shell=True)
    if 'temp.txt' in os.listdir(output+'temp/input/'+db_GI):
        with open(output+'temp/input/'+db_GI+'/temp.txt','r') as tb:
            tb_result=tb.readlines()
            if len(tb_result)!=0:
                tb_list=[]
                tb_done_list=[]
                for line in tb_result:
                    line_name=line.split()[0]
                    if line_name not in tb_list:
                        tb_list.append(line_name)
                for line2 in tb_result:
                    line_list=line2.split()
                    if float(line_list[1])>=id and (line_list[0] not in tb_done_list):
                        total_hit_score+=1
                        tb_done_list.append(line_list[0])
        subprocess.run('rm '+output+'temp/input/'+db_GI+'/temp.txt',shell=True)
    with open(working+'db/blastdb/'+db_GI + '.fasta','r') as fgiscore2:
        if len(fgiscore2.read())>0:
            subprocess.run('blastp -db '+output+'temp/' + GI + ' -query '+working+'db/blastdb/'+db_GI +
                           '.fasta -out '+output+'temp/ref/' + db_GI + '/temp.txt -evalue ' + str(evalue)
                           + ' -outfmt "6 qaccver pident"', shell=True)
    if 'temp.txt' in os.listdir(output+'temp/ref/'+db_GI):
        with open(output+'temp/ref/' + db_GI + '/temp.txt', 'r') as tb2:
            tb_result2 = tb2.readlines()
            if len(tb_result2) != 0:
                tb_list2 = []
                tb_done_list2 = []
                for line4 in tb_result2:
                    line_name2 = line4.split()[0]
                    if line_name2 not in tb_list2:
                        tb_list2.append(line_name2)
                for line3 in tb_result2:
                    line_list2 = line3.split()
                    if float(line_list2[1]) >= id and (line_list2[0] not in tb_done_list2):
                        total_completeness_score += 1
                        tb_done_list2.append(line_list2[0])
        subprocess.run('rm '+output+'temp/ref/' + db_GI + '/temp.txt',shell=True)
    if total_input_length==0:
        hit_score=0
    else:
        hit_score=100*total_hit_score/total_input_length
    if total_ref_length==0:
        completeness_score=0
    else:
        completeness_score=100*total_completeness_score/total_ref_length
    total_score=hit_score+completeness_score
    subprocess.run('rm -rf '+output+'temp/ref/'+db_GI,shell=True)
    subprocess.run('rm -rf '+output+'temp/input/'+db_GI,shell=True)
    outputresult={}
    outputresult['name'] = term
    outputresult['accession'] = acc
    outputresult['hit_score'] = hit_score
    outputresult['complete_score'] = completeness_score
    outputresult['total_score'] = total_score
    with open(output+'temp/compare/'+GI+'/'+db_GI+'.json','w') as jsonout:
        json.dump(outputresult,jsonout)

# Generates a list of subsequences
# Takes in the full sequence, and the length of subsequences
# Outputs a list, each member indicates a subsequences
# Each subsequence is a list, in the form of [sequence, Start of mother_sequence, End of mother_sequence]
def repeat_generator(seq,length,inverted):
    number=len(seq)-length+1
    outseq=[]
    if number<1:
        return outseq
    else:
        if inverted==0:
            for n in range(number):
                n2=n+length
                outseq.append([seq[n:n2],n+1,n2])
        else:
            for n in range(number):
                n2=n+length
                outseq.append([seq[n:n2],reverse_position(len(seq),n2),reverse_position(len(seq),n+1)])
        return outseq

# Caculates the number of mismatches
# inputs two sequences of the same length
# Outputs the number of mismatches
def mismatch(seq1,seq2):
    length=len(seq1)
    mm=0
    for n in range(length):
        if seq1[n]!=seq2[n]:
            mm+=1
    return mm

# Calculates matched subsequences of a desired length from two sequences
# Allows mismatches
# outputs a list of sequence matches, each sequence object is defined in repeat_generator
def rm(f1,f2,length,threshold,inverted):
    seq_f1=repeat_generator(f1,length,0)
    seq_f2=repeat_generator(f2,length,inverted)
    match_out=[]
    matches=[]
    if seq_f1==[] or seq_f2==[]:
        return match_out
    else:
        for s1 in seq_f1:
            for s2 in seq_f2:
                if mismatch(s1[0],s2[0])<=threshold:
                    if s1[-1]!='seq1':
                        s1.append('seq1')
                        matches.append(s1)
                    if s2[-1]!='seq2':
                        s2.append('seq2')
                        matches.append(s2)
        for match in matches:
            if match not in match_out:
                match_out.append(match)
        return match_out

# Performs clustering of matched repeats
# Inputs matches, and allowed mismatches
# Outputs a dictionary of clusters, format dict[cluster1]=[seq1,seq2,seq3], sequence object is defined in repeat_generator
def clustering(matches,threshold):
    cluster_number=1
    clusters={}
    if matches==[]:
        return clusters
    else:
        clusters['cluster_'+str(cluster_number)]=[]
        clusters['cluster_' + str(cluster_number)].append(matches[0])
        for match in matches[1:]:
            clustered_1=0
            for key in clusters.keys():
                if mismatch(match[0],clusters[key][0][0])<=threshold:
                    clusters[key].append(match)
                    clustered_1=1
                    break
            if clustered_1==0:
                cluster_number+=1
                clusters['cluster_'+str(cluster_number)]=[]
                clusters['cluster_' + str(cluster_number)].append(match)
    return clusters

# Finds closest repeat pairs for each cluster
# Inputs clustering results
# Outputs closest repeat pairs for each cluster, sequence object is defined in repeat_generator
def best_repeat_pairs(f1,f2,length,threshold,inverted):
    matches=rm(f1,f2,length,threshold,inverted)
    clusters=clustering(matches,threshold)
    repeat_pairs={}
    if clusters=={}:
        return repeat_pairs
    else:
        for key2 in clusters.keys():
            group_s1=[]
            group_s2=[]
            for item in clusters[key2]:
                if item[-1]=='seq1':
                    group_s1.append(item)
                else:
                    group_s2.append(item)
            if len(group_s1)>=1 and len(group_s2)>=1:
                s1_closest=group_s1[0]
                s2_closest=group_s2[0]
                if len(group_s1)>1:
                    for item_s1 in group_s1[1:]:
                        if item_s1[1]>s1_closest[1]:
                            s1_closest=item_s1
                if len(group_s2)>1:
                    for item_s2 in group_s2[1:]:
                        if item_s2[1]<s2_closest[1]:
                            s2_closest=item_s2
                repeat_pairs[key2]=[s1_closest,s2_closest]
        return repeat_pairs

# Identify GIs, used for iteration, outputs json files.
def identify_GI(arg):
    integrases=arg[0]
    integrase_key=arg[1]
    size=int(arg[2])
    genome_sequence=arg[3]
    minlen=int(arg[4])
    mismatch=int(arg[5])
    size_spacer=int(arg[6])
    output=arg[7]
    minlen_gi=int(arg[8])
    repeat_forward_a={}
    repeat_reverse_a={}
    repeat_forward_b={}
    repeat_reverse_b={}
    integrase=integrases[integrase_key]
    start_integrase=int(integrase[1])
    end_integrase=int(integrase[2])
    seqs=create_seqs(integrase,size,size_spacer,genome_sequence,start_integrase,end_integrase)
    seq1a=seqs[0]
    seq2a=seqs[1]
    seq2a_rev=seqs[2]
    start_fragment_a=int(seqs[3])
    end_fragment_a=int(seqs[4])
    seq1b=seqs[5]
    seq2b=seqs[6]
    seq2b_rev=seqs[7]
    start_fragment_b=int(seqs[8])
    end_fragment_b=int(seqs[9])
    # Identify repeat pairs
    rep_pairs_forward_a=best_repeat_pairs(seq1a,seq2a,minlen,mismatch,0)
    rep_pairs_forward_b=best_repeat_pairs(seq1b,seq2b,minlen,mismatch,0)
    rep_pairs_reverse_a=best_repeat_pairs(seq1a,seq2a_rev,minlen,mismatch,1)
    rep_pairs_reverse_b=best_repeat_pairs(seq1b,seq2b_rev,minlen,mismatch,1)
    # Get information on each genomic island
    if len(rep_pairs_forward_a.keys())!=0:
        for key_forward_a in rep_pairs_forward_a.keys():
            forw_a=[]
            forw_a.append(rep_pairs_forward_a[key_forward_a][0][:-1])
            forw_a.append(rep_pairs_forward_a[key_forward_a][1][:-1])
            forw_out_a=[]
            forw_out_a.append('direct_repeat')
            forw_out_a.append(integrase_key)
            forw_out_a.append(str(integrase[0]))
            forw_out_a.append(str(start_integrase))
            forw_out_a.append(str(end_integrase))
            forw_out_a.append(str(start_fragment_a+int(forw_a[0][1])-1))
            forw_out_a.append(str(end_integrase+int(forw_a[1][2])))
            forw_a.append(forw_out_a)
            if end_integrase+int(forw_a[1][2])-start_fragment_a-int(forw_a[0][1])+1>=minlen_gi:
                repeat_forward_a[key_forward_a] = forw_a
    if len(rep_pairs_forward_b.keys())!=0:
        for key_forward_b in rep_pairs_forward_b.keys():
            forw_b=[]
            forw_b.append(rep_pairs_forward_b[key_forward_b][0][:-1])
            forw_b.append(rep_pairs_forward_b[key_forward_b][1][:-1])
            forw_out_b=[]
            forw_out_b.append('direct_repeat')
            forw_out_b.append(integrase_key)
            forw_out_b.append(str(integrase[0]))
            forw_out_b.append(str(start_integrase))
            forw_out_b.append(str(end_integrase))
            forw_out_b.append(str(start_fragment_b+int(forw_b[0][1])-1))
            forw_out_b.append(str(end_integrase+int(forw_b[1][2])))
            forw_b.append(forw_out_b)
            if end_integrase+int(forw_b[1][2])-start_fragment_b-int(forw_b[0][1])+1>=minlen_gi:
                repeat_forward_b[key_forward_b] = forw_b
    if len(rep_pairs_reverse_a.keys())!=0:
        for key_reverse_a in rep_pairs_reverse_a.keys():
            rev_a=[]
            rev_a.append(rep_pairs_reverse_a[key_reverse_a][0][:-1])
            rev_a.append(rep_pairs_reverse_a[key_reverse_a][1][:-1])
            rev_out_a=[]
            rev_out_a.append('inverted_repeat')
            rev_out_a.append(integrase_key)
            rev_out_a.append(str(integrase[0]))
            rev_out_a.append(str(start_integrase))
            rev_out_a.append(str(end_integrase))
            rev_out_a.append(str(start_fragment_a+int(rev_a[0][1])-1))
            rev_out_a.append(str(end_integrase+int(rev_a[1][2])))
            rev_a.append(rev_out_a)
            if end_integrase + int(rev_a[1][2]) - start_fragment_a - int(rev_a[0][1]) + 1 >= minlen_gi:
                repeat_reverse_a[key_reverse_a]=rev_a
    if len(rep_pairs_reverse_b.keys())!=0:
        for key_reverse_b in rep_pairs_reverse_b.keys():
            rev_b=[]
            rev_b.append(rep_pairs_reverse_b[key_reverse_b][0][:-1])
            rev_b.append(rep_pairs_reverse_b[key_reverse_b][1][:-1])
            rev_out_b=[]
            rev_out_b.append('inverted_repeat')
            rev_out_b.append(integrase_key)
            rev_out_b.append(str(integrase[0]))
            rev_out_b.append(str(start_integrase))
            rev_out_b.append(str(end_integrase))
            rev_out_b.append(str(start_fragment_b+int(rev_b[0][1])-1))
            rev_out_b.append(str(end_integrase+int(rev_b[1][2])))
            rev_b.append(rev_out_b)
            if end_integrase+int(rev_b[1][2]) - start_fragment_b-int(rev_b[0][1])+1 >= minlen_gi:
                repeat_reverse_b[key_reverse_b]=rev_b
    with open(output+'temp/json/'+integrase_key+'_forward_a.json','w') as ff_a:
        json.dump(repeat_forward_a,ff_a)
    with open(output+'temp/json/'+integrase_key+'_reverse_a.json','w') as fr_a:
        json.dump(repeat_reverse_a,fr_a)
    with open(output+'temp/json/'+integrase_key+'_forward_b.json','w') as ff_b:
        json.dump(repeat_forward_b,ff_b)
    with open(output+'temp/json/'+integrase_key+'_reverse_b.json','w') as fr_b:
        json.dump(repeat_reverse_b,fr_b)

# Connect two strings that are 'step' bases different
def connect(str1,str2,step):
    return str1+str2[-1*step:]

def smaller(num1,num2):
    if num1<num2:
        return num1
    else:
        return num2

def bigger(num1,num2):
    if num1>num2:
        return num1
    else:
        return num2

def bubblesort(repeat):
    n=len(repeat)
    for i in range(n):
        for j in range(n-i-1):
            if repeat[j][0][2]>repeat[j+1][0][2]:
                repeat[j],repeat[j+1]=repeat[j+1],repeat[j]
    return repeat

def repeat_connector(repeat_list2,step):
    repeat_connection_indicator = 0
    while repeat_connection_indicator == 0:
        repeat_connection_indicator = 1
        repeat_list3 = bubblesort(repeat_list2)
        Flag=0
        for repeat in repeat_list3:
            for repeat1 in repeat_list3:
                if int(repeat[0][2]) == int(repeat1[0][2]) - step and repeat[2][0] == repeat1[2][0]:
                    repeat_connection_indicator = 0
                    seq_seq1 = connect(repeat[0][0], repeat1[0][0], step)
                    seq_seq2 = connect(repeat[1][0], repeat1[1][0], step)
                    start_seq1 = smaller(repeat[0][1], repeat1[0][1])
                    end_seq1 = bigger(repeat[0][2], repeat1[0][2])
                    start_seq2 = smaller(repeat[1][1], repeat1[1][1])
                    end_seq2 = bigger(repeat[1][2], repeat1[1][2])
                    gi_start = smaller(int(repeat[2][5]), int(repeat1[2][5]))
                    gi_end = bigger(int(repeat[2][6]), int(repeat1[2][6]))
                    seq1 = []
                    seq1.append(seq_seq1)
                    seq1.append(start_seq1)
                    seq1.append(end_seq1)
                    seq2 = []
                    seq2.append(seq_seq2)
                    seq2.append(start_seq2)
                    seq2.append(end_seq2)
                    info = []
                    info.append(repeat[2][0])
                    info.append(repeat[2][1])
                    info.append(repeat[2][2])
                    info.append(repeat[2][3])
                    info.append(repeat[2][4])
                    info.append(str(gi_start))
                    info.append(str(gi_end))
                    new_repeat = []
                    new_repeat.append(seq1)
                    new_repeat.append(seq2)
                    new_repeat.append(info)
                    repeat_list2.append(new_repeat)
                    if repeat in repeat_list2:
                        repeat_list2.remove(repeat)
                    if repeat1 in repeat_list2:
                        repeat_list2.remove(repeat1)
                    Flag=1
                    break
            if Flag==1:
                break
    return repeat_list2

# Connect adjacent repeats so that maximum repeats can be achieved
def maximize_repeats(repeats,mismatch):
    integrase_list=[]
    repeat_out={}
    for key in repeats.keys():
        if repeats[key][2][1] not in integrase_list:
            integrase_list.append(repeats[key][2][1])
    repeats_separated={}
    for integrase in integrase_list:
        repeats_separated[integrase]=[]
        for key2 in repeats.keys():
            if repeats[key2][2][1]==integrase:
                repeats_separated[integrase].append(repeats[key2])
    for key4 in repeats_separated.keys():
        repeat_list2=repeats_separated[key4]
        repeat_list2=repeat_connector(repeat_list2,1)
        if mismatch>=2:
            step1=2
            while step1<=mismatch:
                repeat_list2=repeat_connector(repeat_list2,step1)
                step1+=1
        repeat_key4_num=1
        for item in repeat_list2:
            repeat_out['repeat_'+key4+'_'+str(repeat_key4_num)]=item
            repeat_key4_num+=1
    repeat_out=remove_small_GI(repeat_out)
    return repeat_out

# Obtain sequences from the contig number, first base number, last base number
def get_sequence(genome,contig,start,end):
    if start>end:
        seq=''
    elif start==end:
        seq=genome[contig][-1]
    else:
        if start>len(genome[contig]):
            seq=''
        elif start==len(genome[contig]):
            seq=genome[contig][-1]
        else:
            if end>=len(genome[contig]):
                seq=genome[contig][start-1:]
            else:
                seq=genome[contig][start-1:end]
    return seq

# Calculate forward position of a reverse strand
def reverse_position(length,position):
    return length+1-position