# !/usr/bin/env python
# coding=utf-8
# author=Mingyu Wang and Mengge Zhang
# Dependencies seqkit, prodigal,biopython,

import datetime
import json
import subprocess
import os
from Bio import SeqIO
import argparse

#processing of a GI record
def RecordProcess(giacc,dbseq_folder):
    gi_file=dbseq_folder+giacc+'.gb'
    gi=list(SeqIO.parse(gi_file,'genbank'))[0]
    gi_seq=str(gi.seq)
    genes={}
    protein_counter=0
    for feature in gi.features:
        if feature.type=='CDS':
            protein_counter+=1
            name='Protein'+str(protein_counter)
            protein=''
            if 'translation' in feature.qualifiers.keys():
                protein=feature.qualifiers['translation'][0]
            elif 'pseudogene' not in feature.qualifiers.keys() and 'pseudo' not in feature.qualifiers.keys():
                print(giacc+' may have a missing translation \n')
                error_file.write(giacc+' may have a missing translation \n')
            if protein!='':
                genes[name]=protein
    if len(genes.keys())==0:
        print(giacc)
        print(gi_seq)
        temp=open('./temp/temp2.fasta','w')
        temp.write('>temp\n')
        temp.write(gi_seq)
        temp.close()
        subprocess.run('prodigal -c -p meta -i ./temp/temp2.fasta -a ./temp/temp2.faa',shell=True)
        subprocess.run('seqkit split -i ./temp/temp2.faa', shell=True)
        gene_list2 = os.listdir('./temp/temp2.faa.split/')
        counter3=0
        for gene2 in gene_list2:
            counter3+=1
            gene_name2 = 'Protein'+str(counter3)
            gene_line=open('./temp/temp2.faa.split/' + gene2, 'r').readlines()[1:]
            gene_line2=[]
            for line3 in gene_line:
                gene_line2.append(line3.strip('\n'))
            genes[gene_name2] = ''.join(gene_line2)
        subprocess.run('rm ./temp/temp2.fasta',shell=True)
        subprocess.run('rm ./temp/temp2.faa',shell=True)
        subprocess.run('rm -rf ./temp/temp2.faa.split/',shell=True)
    with open('./db/blastdb/'+giacc+'.fasta','w') as tp:
        for gene3 in genes.keys():
            tp.write('>'+gene3+'\n')
            tp.write(genes[gene3]+'\n')
    subprocess.run('makeblastdb -dbtype prot -in ./db/blastdb/'+giacc+'.fasta -out ./db/blastdb/'+giacc,shell=True)
    return(gi_seq,genes)

#Parameters
parser=argparse.ArgumentParser()
parser.add_argument('-d',type=str,required=True,metavar='Path_to_DB_file',help='Path to the file containing database information')

args=parser.parse_args()
dbtxt=args.d
outputdb='./db/gidb.json'
dbseq_folder='./db/GI_database_seqs/'

#Open error log files.
subprocess.run('rm -f makedb_error_log.txt',shell=True)
error_file=open('makedb_error_log.txt','w')

#main program
gidb={}
GIinput=open(dbtxt,'r').readlines()
counter=0
for GIrecord in GIinput[1:]:
    counter+=1
    GIrecord_list=GIrecord.strip('\n').split('\t')
    gitype=GIrecord_list[2]
    if not gitype:
        gitype='Untyped'
    giacc=GIrecord_list[0]
    name=GIrecord_list[1]
    accession=GIrecord_list[3]
    processed_GI=RecordProcess(giacc,dbseq_folder)
    gidb[giacc]={}
    gidb[giacc]['gitype']=gitype
    gidb[giacc]['name']=name
    gidb[giacc]['accession']=accession
    gidb[giacc]['seq']=processed_GI[0]
    gidb[giacc]['genes']=processed_GI[1]

#dumping into json object
with open(outputdb,'w') as jf:
    json.dump(gidb,jf)