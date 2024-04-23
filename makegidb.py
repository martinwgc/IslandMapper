# -*- coding: UTF-8 -*-
# Make database for IslandMapper
# Requires seqkit, prodigal, blast

from Bio import Entrez
from Bio import SeqIO
import datetime
import json
import subprocess
import os

Entrez.api_key = '2a3dfe4f955ca4d470c99bfad2e800175f09'
Entrez.email = 'wangmingyu@sdu.edu.cn'
Entrez.tool = 'AMR_Mining'
file = './db/GI_database.txt'


def remove_repeat(inputlist):
    temp_list = []
    for item in inputlist:
        if item not in temp_list:
            temp_list.append(item)
    return temp_list


# Download protein or nucleotide sequences, for Wang lab.
# Usage: get_sequence(acc, type, start, end), acc: accession number, type: sequence type, ='nuc' or 'prot',
# start: starting position, end:ending position. Start and end have to be both provided or not provided.
def get_sequence(acc, dbtype, start=0, end=0):
    if dbtype != 'nuc' and dbtype != 'prot':
        sequence = 'error'
        print('Please indicate whether you need protein or nucleotide sequences')
    else:
        if dbtype == 'nuc':
            database = 'nuccore'
        elif dbtype == 'prot':
            database = 'prot'
        try:
            if int(start) != 0:
                handle = Entrez.efetch(db=database, id=acc, rettype='fasta', retmode='text', seq_start=int(start),
                                       seq_stop=int(end))
                sequence = str(SeqIO.read(handle, format='fasta').seq)
            else:
                handle = Entrez.efetch(db=database, id=acc, rettype='fasta', retmode='text')
                sequence = str(SeqIO.read(handle, format='fasta').seq)
        except IOError:
            sequence = 'error'
    return sequence


def protnblast(giacc, gi_seq):
    temp = open('./temp/temp2.fasta', 'w')
    temp.write('>temp\n')
    temp.write(gi_seq)
    temp.close()
    subprocess.run('prodigal -c -p meta -i ./temp/temp2.fasta -a ./temp/temp2.faa', shell=True)
    subprocess.run('seqkit split -i ./temp/temp2.faa', shell=True)
    gene_list2 = os.listdir('./temp/temp2.faa.split/')
    counter3 = 0
    proteins = {}
    for gene2 in gene_list2:
        counter3 += 1
        gene_name2 = 'Protein' + str(counter3)
        gene_line = open('./temp/temp2.faa.split/' + gene2, 'r').readlines()[1:]
        gene_line2 = []
        for line3 in gene_line:
            gene_line2.append(line3.strip('\n'))
        proteins[gene_name2] = ''.join(gene_line2)
    subprocess.run('rm ./temp/temp2.fasta', shell=True)
    subprocess.run('rm ./temp/temp2.faa', shell=True)
    subprocess.run('rm -rf ./temp/temp2.faa.split/', shell=True)
    return proteins


gi = open(file, 'r').readlines()
gidb = {}
for line in gi[1:]:
    properties = line[:-1].split('\t')
    repeat_list = []
    if properties[5] != 'NA':
        repeat_list = repeat_list+properties[5].split('|')
    if properties[6] != 'NA':
        repeat_acc_list = properties[6].split('|')
        for repeat_acc in repeat_acc_list:
            if ':' in repeat_acc:
                repeat_ncbi_acc = repeat_acc.split(':')[0]
                repeat_ncbi_start = repeat_acc.split(':')[1].split('-')[0]
                repeat_ncbi_end = repeat_acc.split(':')[1].split('-')[1]
            else:
                repeat_ncbi_acc = repeat_acc
                repeat_ncbi_start = 0
                repeat_ncbi_end = 0
            repeat_seq = get_sequence(repeat_ncbi_acc, 'nuc', repeat_ncbi_start, repeat_ncbi_end)
            repeat_list.append(repeat_seq)
    repeat_list = remove_repeat(repeat_list)

    giacc = properties[0]
    print(giacc)
    print(datetime.datetime.now())
    gidb[giacc] = {}
    gidb[giacc]['acc'] = properties[1]
    gidb[giacc]['term'] = properties[2]
    gidb[giacc]['genebank'] = properties[3]
    gidb[giacc]['pmid'] = properties[4]
    gidb[giacc]['repeat'] = repeat_list
    gidb[giacc]['transferability'] = properties[7]
    gidb[giacc]['evidence'] = properties[7]
    gidb[giacc]['antibiotic-resistance'] = properties[9]
    gidb[giacc]['metal-resistance'] = properties[10]
    gidb[giacc]['pathogenicity'] = properties[11]
    if ':' in properties[3]:
        gi_ncbi_acc = properties[3].split(':')[0]
        gi_ncbi_start = properties[3].split(':')[1].split('-')[0]
        gi_ncbi_end = properties[3].split(':')[1].split('-')[1]
    else:
        gi_ncbi_acc = properties[3]
        gi_ncbi_start = 0
        gi_ncbi_end = 0
    gi_seq = get_sequence(gi_ncbi_acc, 'nuc', gi_ncbi_start, gi_ncbi_end)
    gidb[giacc]['sequence'] = gi_seq
    proteins = protnblast(giacc, gi_seq)
    gidb[giacc]['proteins'] = proteins
    with open('./db/blastdb/'+giacc+'.fasta','w') as fo:
        for protein_key in gidb[giacc]['proteins'].keys():
            fo.write('>'+protein_key+'\n')
            fo.write(gidb[giacc]['proteins'][protein_key]+'\n')
with open('./db/gidb.json', 'w') as jsondump:
    json.dump(gidb, jsondump)