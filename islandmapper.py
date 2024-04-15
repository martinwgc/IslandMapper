# -*- coding: UTF-8 -*-
# Search for potential Genomic Islands.
# IslandMapper v0.2.2
# Author=Mingyu Wang, Mengge Zhang, Shuyao Wang
# Requirement: blast >=2.15.0+, prodigal >=2.5.3, tqdm >=4.64.1

import argparse
import os
from tools import *
import datetime
from multiprocessing import Pool
from tqdm import tqdm
import json


# Prints error message
def err_callback(err):
    print(f'error:{str(err)}')


#Arguments and preparation
parser=argparse.ArgumentParser()
parser.add_argument('-i', type=str, required=True, metavar='Path to the input genomic sequence file')
parser.add_argument('-w',type=str,required=False,default='./',metavar='Path to folder of IslandMapper')
parser.add_argument('--minlen',type=int,required=False,default=15,metavar='Minimum length of repeats')
parser.add_argument('--minlen_gi',type=int,required=False,default=10000,metavar='Minimum length of genomic islands')
parser.add_argument('-t',type=int,required=False,default=20,metavar='Number of threads to use')
parser.add_argument('--size',type=int,required=False,default=100000,metavar='Size to search on either side of integrase gene')
parser.add_argument('--sizespacer',type=int,required=False,default=500,metavar='Maximum size between repeats and integrases')
parser.add_argument('--mismatch',type=int,required=False,default=2,metavar='Mismatches allowed in defining repeats')
parser.add_argument('--id',type=int,required=False,metavar='Identity_Threshold',help='Identity threshold for sequence comparison', default=60)
parser.add_argument('--db',type=str,required=False,metavar='Json_path',help='Path for GI database file, in json format',default='./db/gidb.json')
parser.add_argument('-e',type=float,required=False,metavar='evalue',help='evalue threshold when using blast algorithm',default=1e-10)
args=parser.parse_args()
genome=args.i
working=args.w
if working[-1]!='/':
    working=working+'/'
output=working+'output/'
minlen=args.minlen
minlen_gi=args.minlen_gi
threads=args.t
size=args.size
size_spacer=args.sizespacer
mismatch=args.mismatch
id=args.id
db=args.db
evalue=args.e
if not os.path.exists(output):
    os.makedirs(output)
subprocess.run('rm -rf '+output,shell=True)
if not os.path.exists(output+'temp'):
    os.makedirs(output+'temp')
if not os.path.exists(output+'temp/compare'):
    os.makedirs(output+'temp/compare')
if not os.path.exists(output+'temp/proteins'):
    os.makedirs(output+'temp/proteins')
if not os.path.exists(output+'temp/GI'):
    os.makedirs(output+'temp/GI')
if not os.path.exists(output+'temp/json'):
    os.makedirs(output+'temp/json')
if not os.path.exists(output+'genome'):
    os.makedirs(output+'genome')
#output GI table
outputtable=output+'GI.tsv'
outputtable_removed=output+'removed_GI.tsv'

#Main program
if __name__=='__main__':
    # Load GI database
    with open(db, 'r') as fi:
        gidb2 = json.load(fi)
    gidb={}
    for key_gidb in gidb2.keys():
        if gidb2[key_gidb]['proteins']!={}:
            gidb[key_gidb]=gidb2[key_gidb]

    # Creates blast library for genomic islands
    if 'GI000001.phr' not in os.listdir(working+'db/blastdb/'):
        print('==========FIRST TIME USE, PREPARING GI DATABASES==========')
        print(datetime.datetime.now())
        for gi_key in gidb.keys():
            gi_proteins=gidb[gi_key]['proteins']
            construct_blast_library(gi_proteins,gi_key,output,working)
        print('==========GI DATABASES PREPARED==========')
        print(datetime.datetime.now())

    # Predict genes in the genome, outputs protein fasta file
    print('==========STARTING GENE PREDICTION==========')
    print(datetime.datetime.now())
    subprocess.run('prokka --quiet --force --outdir '+output+'genome'+' --cpus '+ str(threads)+' '+genome, shell=True)

    # Identify integrases in the genome, dictionary format: integrases[gene]=[contig,start,end]
    print('==========GENES PREDICTED==========')
    print('==========IDENTIFYING INTEGRASES/TRANSPOSASES==========')
    print(datetime.datetime.now())
    fastafolder=output+'genome/'
    fastafolderfiles=os.listdir(fastafolder)
    for fastafile in fastafolderfiles:
        if fastafile[-4:]=='.faa':
            fasta=output+'genome/'+fastafile
        if fastafile[-4:]=='.gff':
            gff=output+'genome/'+fastafile
    integrases=identify_integrase(fasta,gff,threads,evalue,id,output,working)
    print('==========INTEGRASES/TRANSPOSASES PREDICTED==========')
    print('==========CALCULATING GIs==========')
    print(datetime.datetime.now())

    # Creates a dictionary of contig versus sequences, dictionary format: genome_sequence[contig]=sequence
    genome_sequence=load_sequence(genome)
    contigs=genome_sequence.keys()

    # Identify genomes
    pool=Pool(threads)
    iterated=[]
    for integrase_key in integrases.keys():
        integrase_args=[]
        integrase_args.append(integrases)
        integrase_args.append(integrase_key)
        integrase_args.append(size)
        integrase_args.append(genome_sequence)
        integrase_args.append(minlen)
        integrase_args.append(mismatch)
        integrase_args.append(size_spacer)
        integrase_args.append(output)
        integrase_args.append(minlen_gi)
        iterated.append(integrase_args)
    pbar=tqdm(total=len(iterated))
    pbar.set_description('Identifying GI')
    update=lambda *args: pbar.update()
    for analyzed_integrase in iterated:
        pool.apply_async(func=identify_GI,args=(analyzed_integrase,),callback=update,error_callback=err_callback)
    pool.close()
    pool.join()

    #Removes GIs that are part of other GIs
    jsonlist=os.listdir(output+'temp/json')
    repeats={}
    repeat_count=1
    for jsonfile in jsonlist:
        with open(output+'temp/json/'+jsonfile,'r') as fj:
            json_dict=json.load(fj)
            for key in json_dict:
                repeats['repeat'+str(repeat_count)]=json_dict[key]
                repeat_count+=1
    repeats=maximize_repeats(repeats,mismatch)
    repeats=remove_small_GI(repeats)
    separate_results=separate_repeats(repeats)
    repeats_removed=separate_results[1]
    repeats=separate_results[0]
    print('==========GIs calculated==========')
    print(datetime.datetime.now())

    # Outputs GI tables
    print('==========STARTING TO OUTPUT PUTATIVE GENOMIC ISLANDS==========')
    with open(outputtable,'w') as fo:
        if len(repeats.keys())==0:
            fo.write('nothing found')
        else:
            fo.write('GI_num\trepeat_type\tintegrase_name\tcontig\tintegrase_start\tintegrase_end\tGI_start\tGI_end\n')
            gi=1
            for repfkey in repeats.keys():
                fo.write('GI_'+str(gi)+'\t'+'\t'.join(repeats[repfkey][2])+'\n')
                with open(output+'temp/GI/GI_'+str(gi)+'.fasta','w') as fastf:
                    start_gi=int(repeats[repfkey][2][5])
                    end_gi=int(repeats[repfkey][2][6])
                    fastf.write('>'+'GI_'+str(gi)+'\n')
                    fastf.write(get_sequence(genome_sequence,repeats[repfkey][2][2],start_gi,end_gi))
                gi+=1
    with open(outputtable_removed,'w') as fo_removed:
        if len(repeats_removed.keys())==0:
            fo_removed.write('nothing found')
        else:
            fo_removed.write('GI_num\trepeat_type\tintegrase_name\tcontig\tintegrase_start\tintegrase_end\tGI_start\tGI_end\n')
            gi_removed=1
            for repfkey_removed in repeats_removed.keys():
                fo_removed.write('GI_removed_'+str(gi_removed)+'\t'+'\t'.join(repeats_removed[repfkey_removed][2])+'\n')
    print('==========PUTATIVE GENOMIC ISLANDS CALCULATED==========')
    print(datetime.datetime.now())

    # Starting to compare with known GIs
    print('==========STARTING COMPARISON WITH KNOWN GENOMIC ISLANDS==========')
    GI_list=[]
    for item in os.listdir(output+'temp/GI/'):
        GI_list.append(item)
    if not os.path.exists(output+'temp/input/'):
        os.makedirs(output+'temp/input/')
    if not os.path.exists(output+'temp/ref/'):
        os.makedirs(output+'temp/ref/')
    print('==========WORKING ON GI GENE PREDICTION AND BLAST DATABASE==========')
    print(datetime.datetime.now())
    for GI in GI_list:
        print(datetime.datetime.now())
        with open(output + 'temp/GI/' + GI, 'r') as gi_count:
            gi_list = gi_count.readlines()[1]
            subprocess.run(
                    'prodigal -i ' + output + 'temp/GI/' + GI + ' -a ' + output + 'temp/proteins/' + GI + ' -o ' +
                    output + 'temp/proteins/' + GI + '.gff -q -p meta', shell=True)
        if os.path.exists(output + 'temp/proteins/' + GI):
            with open(output + 'temp/proteins/' + GI, 'r') as fgi:
                fgi_string = fgi.read()
                if fgi_string.count('>') != 0:
                    subprocess.run(
                        'makeblastdb -dbtype prot -in ' + output + 'temp/proteins/' + GI + ' -out ' + output + 'temp/' + GI,
                        shell=True)
    print('==========GI BLAST DATABASES MADE==========')
    print(datetime.datetime.now())
    print('==========STARTING TO COMPARE WITH KNOWN GIS==========')
    for GI2 in tqdm(GI_list, desc='Comparing with GI'):
        if os.path.exists(output + 'temp/proteins/' + GI2):
            with open(output + 'temp/proteins/' + GI2, 'r') as fgi2:
                genes = {}
                fgi_string2 = fgi2.read()
                fgi_list = fgi_string2.split('>')[1:]
                for fgi_line in fgi_list:
                    fgi_name = fgi_line.split('\n')[0]
                    fgi_seq = ''.join(fgi_line.split('\n')[1:])
                    genes[fgi_name] = fgi_seq
        compare_result = {}
        compare_result[GI2] = {}
        database_GIs = gidb.keys()
        if not os.path.exists(output + 'temp/compare/'+GI2):
            os.makedirs(output + 'temp/compare/'+GI2)
        iterated2=[]
        for db_GI in database_GIs:
            db_genes = gidb[db_GI]['proteins']
            args=[]
            args.append(id)
            args.append(db_GI)
            args.append(genes)
            args.append(db_genes)
            args.append(evalue)
            args.append(GI2)
            args.append(threads)
            args.append(output)
            args.append(gidb[db_GI]['term'])
            args.append(gidb[db_GI]['acc'])
            args.append(working)
            iterated2.append(args)
        pool2 = Pool(threads)
        for compared_gi in iterated2:
            pool2.apply_async(func=giscore_blast, args=(compared_gi,), error_callback=err_callback)
        pool2.close()
        pool2.join()
        db_GI_list=os.listdir(output+'temp/compare/'+GI2)
        for db_GI2 in db_GI_list:
            with open(output+'temp/compare/'+GI2+'/'+db_GI2,'r') as jsoninf:
                jsonin=json.load(jsoninf)
                compare_result[GI2][db_GI2] = {}
                compare_result[GI2][db_GI2]['name'] = jsonin['name']
                compare_result[GI2][db_GI2]['accession'] = jsonin['accession']
                compare_result[GI2][db_GI2]['hit_score'] = jsonin['hit_score']
                compare_result[GI2][db_GI2]['complete_score'] = jsonin['complete_score']
                compare_result[GI2][db_GI2]['total_score'] = jsonin['total_score']
        fo2=open(output + GI2[:-6] + '_result.txt', 'a')
        fo2.write('GI\tname\taccession\thit_score\tcomplete_score\ttotal_score\n')
        for key_o in compare_result[GI2].keys():
            fo2.write(key_o[:-5] + '\t'  + compare_result[GI2][key_o]['name'] \
                                 + '\t' + compare_result[GI2][key_o]['accession'] + '\t' + str(
                            compare_result[GI2][key_o]['hit_score']) + \
                                 '\t' + str(compare_result[GI2][key_o]['complete_score']) + '\t' + str(
                            compare_result[GI2][key_o]['total_score']) + '\n')
    list1 = os.listdir(output)
    output_final = open(output+'compiled_results.tsv', 'w')
    output_final.write('GI_input\tGI\tname\taccession\thit_score\tcomplete_score\ttotal_score\n')
    output_list=[]
    for item1 in list1:
        if '_result.txt' in item1:
            output_list.append(item1)
    for item2 in output_list:
        name = '_'.join(item2.split('_')[:2])
        with open(output + item2, 'r') as f_o:
            result = f_o.readlines()[1:]
            for line_result in result:
                line_list = line_result.split()
                if float(line_list[-1]) != 0:
                    output_final.write(name + '\t' + line_result)
    print(datetime.datetime.now())
    print('==========ALL DONE, THANK YOU FOR USING ISLANDMAPPER==========')