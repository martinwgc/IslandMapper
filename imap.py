# !/usr/bin/env python
# coding=utf-8
# IslandMapper v0.1
# Author=Mingyu Wang and Mengge Zhang
# Dependencies: seqkit, EMBOSS suite (seqret,needle), prodigal,blast
# Need to be able to create folders where output_folder is

import subprocess
import json
import os
import argparse
import datetime

# Parameters
parser=argparse.ArgumentParser()
parser.add_argument('-i',type=str,required=True,metavar='Input_file_path',help='Path for input file,must be fasta format of nucleotide sequences')
parser.add_argument('-o',type=str,required=True,metavar='Output_file_path',help='Path for output file')
parser.add_argument('-a',type=str,required=True,choices=['blast','needle'],metavar='Algorithm_of_choice',\
                    help='Choose between EMBOSS needle and blast',default='blast')
parser.add_argument('-id',type=int,required=True,metavar='Identity_Threshold',help='Identity threshold for sequence comparison', default=60)
parser.add_argument('-cov',type=str,required=True,metavar='Coverage_Threshold',help='sequence coverage threshold when using blast algorithm',default=90)
parser.add_argument('-t',type=int,metavar='Threads_for_BLAST',help='Number of threads for blast function',default=1)
parser.add_argument('-db',type=str,metavar='Json_path',help='Path for GI database file, in json format',default='./db/gidb.json')
parser.add_argument('-d',type=str,metavar='Output_path',help='Output directory',default='./GIannot/')
parser.add_argument('-e',type=float,metavar='evalue',help='evalue threshold when using blast algorithm',default=1e-5)
args=parser.parse_args()

threads=args.t
identity_threshold=args.id
inputs=args.i
db=args.db
output=args.o
output_folder=args.d
method=args.a
coverage_threshold=args.cov
evalue=args.e

# Parses EMBOSS Needleman-Wunsch results
def ResultParser(result_file):
    list=open(result_file,'r').readlines()
    length=len(list)
    identity=0
    for num in range(length):
        if list[num][:4]=='# Id':
            identity=list[num][-8:-3]
            if identity[0]=='(':
                identity=identity[1:]
            elif identity[:2]==' (':
                identity=identity[2:]
        else:
            continue
    return identity

# Calculates identities-for multiprocessing of needle
def calc_identity(db_GI3,gene_input,gene_ref):
    subprocess.run('needle -asequence ./temp/input/'+db_GI3+'/' + gene_input + ' -bsequence ./temp/ref/' + db_GI3+'/'+gene_ref + \
                            ' -outfile ./temp/tempneedle/temp_' + db_GI3+'_'+gene_ref + '_' + gene_input + '.txt -brief Y \
                                      -gapopen 10.0 -gapextend 0.5 -auto', shell=True)
    identity = float(ResultParser('./temp/tempneedle/temp_' + db_GI3+'_'+gene_ref + '_' + gene_input + '.txt'))
    subprocess.run('rm ./temp/tempneedle/temp_' + db_GI3+'_'+gene_ref + '_' + gene_input + '.txt', shell=True)
    return identity

# Calculate comparison scores with needle
# hit_score: The percentage of input GI genes found in reference GI
# completeness_score: The percentage of reference GI genes found in input GI
def giscore(inputlist):
    id=float(inputlist[0])
    db_GI2=inputlist[1]
    inputdic=inputlist[2]
    refdic=inputlist[3]
    gene_list_input=inputdic.keys()
    gene_list_ref=refdic.keys()
    subprocess.run('mkdir ./temp/input/'+db_GI2+'/',shell=True)
    subprocess.run('mkdir ./temp/ref/'+db_GI2+'/',shell=True)
    for gene_input in gene_list_input:
        with open('./temp/input/'+db_GI2+'/'+gene_input,'w') as f_input:
            f_input.write(inputdic[gene_input])
    for gene_ref in gene_list_ref:
        with open('./temp/ref/'+db_GI2+'/'+gene_ref,'w') as f_ref:
            f_ref.write(refdic[gene_ref])
    comparison_matrix={}
    for gene_input2 in gene_list_input:
        comparison_matrix[gene_input2]={}
        for gene_ref2 in gene_list_ref:
            gene_identity_result=calc_identity(db_GI2,gene_input2,gene_ref2)
            comparison_matrix[gene_input2][gene_ref2]=gene_identity_result
    total_input_length=len(gene_list_input)
    total_hit_score=0
    for gene_input3 in gene_list_input:
        counter=0
        for gene_ref3 in gene_list_ref:
            if float(comparison_matrix[gene_input3][gene_ref3])>float(id):
                counter=1
        total_hit_score=total_hit_score+counter
    hit_score=100*total_hit_score/total_input_length
    total_ref_length=len(gene_list_ref)
    total_completeness_score=0
    for gene_ref4 in gene_list_ref:
        counter2=0
        for gene_input4 in gene_list_input:
            if float(comparison_matrix[gene_input4][gene_ref4])>float(id):
                counter2=1
        total_completeness_score=total_completeness_score+counter2
    completeness_score=100*total_completeness_score/total_ref_length
    total_score=hit_score+completeness_score
    subprocess.run('rm -rf ./temp/ref/'+db_GI2,shell=True)
    subprocess.run('rm -rf ./temp/input/'+db_GI2,shell=True)
    return hit_score,completeness_score,total_score

# Calculate comparison scores with blast
# hit_score: The percentage of input GI genes found in reference GI
# completeness_score: The percentage of reference GI genes found in input GI
def giscore_blast(inputlist):
    id=float(inputlist[0])
    db_GI2=inputlist[1]
    inputdic=inputlist[2]
    refdic=inputlist[3]
    cover=float(inputlist[4])
    evalue=float(inputlist[5])
    GI2=inputlist[6]
    threads=int(inputlist[7])
    gene_list_input=inputdic.keys()
    gene_list_ref=refdic.keys()
    total_input_length=len(gene_list_input)
    total_hit_score=0
    total_ref_length=len(gene_list_ref)
    total_completeness_score=0
    subprocess.run('mkdir ./temp/input/'+db_GI2,shell=True)
    subprocess.run('mkdir ./temp/ref/'+db_GI2,shell=True)
    subprocess.run('blastp -db ./db/blastdb/' + db_GI2 + ' -query '+output_folder+'faa/'+GI2+'.fa -out ./temp/input/' +\
                   db_GI2 + '/temp.txt -evalue ' + str(evalue) + ' -outfmt "6 qaccver pident qcovs"  -num_threads ' + \
                   str(threads), shell=True)
    with open('./temp/input/'+db_GI2+'/temp.txt','r') as tb:
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
                if float(line_list[1])>=id and float(line_list[2])>=cover and (line_list[0] not in tb_done_list):
                    total_hit_score+=1
                    tb_done_list.append(line_list[0])
    subprocess.run('rm ./temp/input/'+db_GI2+'/temp.txt',shell=True)
    subprocess.run('blastp -db ./temp/' + GI2 + ' -query ./db/blastdb/'+db_GI2 + '.fasta -out ./temp/ref/' + db_GI2 \
                   + '/temp.txt -evalue ' + str(evalue) + ' -outfmt "6 qaccver pident qcovs" -num_threads ' \
                   + str(threads), shell=True)
    with open('./temp/ref/' + db_GI2 + '/temp.txt', 'r') as tb2:
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
                if float(line_list2[1]) >= id and float(line_list2[2]) >= cover and (line_list2[0] not in tb_done_list2):
                    total_completeness_score += 1
                    tb_done_list2.append(line_list2[0])
    subprocess.run('rm ./temp/ref/' + db_GI2 + '/temp.txt',shell=True)
    hit_score=100*total_hit_score/total_input_length
    completeness_score=100*total_completeness_score/total_ref_length
    total_score=hit_score+completeness_score
    subprocess.run('rm -rf ./temp/ref/'+db_GI2,shell=True)
    subprocess.run('rm -rf ./temp/input/'+db_GI2,shell=True)
    return hit_score,completeness_score,total_score


# hit_score: The percentage of input GI genes found in reference GI
# completeness_score: The percentage of reference GI genes found in input GI
#Preparation
with open(db,'r') as fi:
    database=json.load(fi)
datafile=open(inputs,'r').readlines()
subprocess.run('mkdir '+output_folder,shell=True)
subprocess.run('mkdir '+output_folder+'faa/',shell=True)
subprocess.run('mkdir '+output_folder+'gb/',shell=True)
subprocess.run('mkdir '+output_folder+'gff/',shell=True)
subprocess.run('mkdir '+output_folder+'results/',shell=True)
subprocess.run('rm -rf ./temp/*',shell=True)
subprocess.run('mkdir ./temp/input/',shell=True)
subprocess.run('mkdir ./temp/ref/',shell=True)
subprocess.run('mkdir ./temp/tempneedle/',shell=True)
subprocess.run('rm error.txt',shell=True)
error_file=open('error.txt','a')
GIs=[]
for data in datafile:
    if data[0]=='>':
        GIs.append(data[1:-1])
# Main program
try:
    for GI in GIs:
        genes={}
        seq_record=subprocess.run('seqkit grep -p '+GI+' '+inputs, shell=True,universal_newlines=True,stdout=subprocess.PIPE)
        seq=seq_record.stdout
        print('---------------------Working on '+GI+'----------------------')
        print('Current time')
        print(datetime.datetime.now())
        with open('./temp/temp.fasta','w') as f:
            f.write(seq)
        # Uses prodigal to predict genes in a given GI.
        if len(seq)>20000:
            subprocess.run('prodigal -c -i ./temp/temp.fasta -a '+output_folder+'faa/'+GI+'.fa -f gff -o ./temp/'+GI+'.gff',shell=True)
        else:
            subprocess.run('prodigal -c -p meta -i ./temp/temp.fasta -a '+output_folder+'faa/'+GI+'.fa -f gff -o ./temp/'+GI+'.gff',shell=True)
        subprocess.run('seqkit split -i '+output_folder+'faa/'+GI+'.fa',shell=True)
        gene_list=os.listdir(output_folder+'faa/'+GI+'.fa.split/')
        for gene in gene_list:
            gene_name=gene.split('.')[-2][3:]
            genes[gene_name]=open(output_folder+'faa/'+GI+'.fa.split/'+gene,'r').read()
        subprocess.run('rm -rf '+output_folder+'faa/'+GI+'.fa.split/',shell=True)
        # Adds /product_id= in gff, to be compatible with Easyfig.
        gff=open('./temp/'+GI+'.gff','r').readlines()
        with open(output_folder+'gff/'+GI+'.gff','w') as f_gff:
            line_gff_o=[]
            for line_gff in gff:
                if line_gff[0]=='#':
                    line_gff_o.append(line_gff)
                else:
                    line_gff_list=line_gff.split('\t')
                    product_id=line_gff_list[-1].split(';')[0].split('_')[-1]
                    line_gff_list[-1]='product='+line_gff_list[0]+'_'+str(product_id)+';'+line_gff_list[-1]
                    line_gff_o.append('\t'.join(line_gff_list))
            f_gff.writelines(line_gff_o)
        subprocess.run('rm ./temp/'+GI+'.gff',shell=True)
        subprocess.run('seqret -sequence ./temp/temp.fasta -feature -fformat gff -fopenfile '+output_folder+'gff/'+GI+'.gff \
            -osformat genbank -osname_outseq '+GI+' -auto -osdirectory2 '+output_folder+'gb/',shell=True)
        subprocess.run('rm ./temp/temp.fasta',shell=True)
        compare_result={}
        compare_result[GI]={}
        database_GIs=database.keys()
        if len(genes.keys())!=0:
            if method=='blast':
                with open(output_folder+'faa/'+GI+'.fa','w') as tg:
                    for gene_key in genes.keys():
                        tg.write(genes[gene_key]+'\n')
                subprocess.run('makeblastdb -dbtype prot -in '+output_folder+'faa/'+GI+'.fa -out ./temp/' + GI,
                           shell=True)
            for db_GI in database_GIs:
                db_genes=database[db_GI]['genes']
                args=[identity_threshold,db_GI,genes,db_genes,coverage_threshold,evalue,GI,threads]
                if method=='needle':
                    compare_scores=giscore(args)
                elif method=='blast':
                    compare_scores=giscore_blast(args)
                compare_result[GI][db_GI]={}
                compare_result[GI][db_GI]['gitype']=database[db_GI]['gitype']
                compare_result[GI][db_GI]['name']=database[db_GI]['name']
                compare_result[GI][db_GI]['accession']=database[db_GI]['accession']
                compare_result[GI][db_GI]['hit_score']=compare_scores[0]
                compare_result[GI][db_GI]['complete_score']=compare_scores[1]
                compare_result[GI][db_GI]['total_score']=compare_scores[2]
            if method == 'blast':
                subprocess.run('rm ./temp/'+GI+'*',shell=True)
            with open(output_folder+'results/'+GI+'_result.txt','w') as fo:
                fo.write('GI\tgitype\tname\taccession\thit_score\tcomplete_score\ttotal_score\n')
                for key in compare_result[GI].keys():
                    fo.write(key+'\t'+compare_result[GI][key]['gitype']+'\t'+compare_result[GI][key]['name']\
                     +'\t'+compare_result[GI][key]['accession']+'\t'+str(compare_result[GI][key]['hit_score'])+\
                     '\t'+str(compare_result[GI][key]['complete_score'])+'\t'+str(compare_result[GI][key]['total_score'])+'\n')
        else:
            continue
except Exception:
    error_file.write(Exception.args)
    error_file.write('\n')
finally:
    if output_folder[-1]=='/':
        result_folder=output_folder+'results/'
    else:
        result_folder=output_folder+'/results/'

    list1 = os.listdir(result_folder)
    output = open(output, 'w')
    output.write('GI_input\tGI\tgitype\tname\taccession\thit_score\tcomplete_score\ttotal_score\n')
    for item in list1:
        print(item)
        name = item.split('_')[0]
        with open(result_folder + item, 'r') as f_o:
            result = f_o.readlines()[1:]
            for line_result in result:
                line_list = line_result.split()
                if float(line_list[-1]) != 0:
                    output.write(name + '\t' + line)
    print('FINALLY...ALL DONE. HAVE A NICE DAY. THANKS A BUNCH')


