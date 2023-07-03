# IslandMapper
Compares a file of putative genomic island nucleotide fasta sequences against a genomic island database to identify similar islands

System requirement: Linux, Python 3 (3.10 tested)
Dependencies-need to be installed with your python: BioPython
Dependencies-need to run from your command line:EMBOSS suite (needle, seqret), seqkit, NCBI-blast+, prodigal

Installation:
No installation needed, simply use python makegidb.py or python imap.py to run from your commandline.

Database:
We do have a preliminary GI database version 2023-07-03 insalled in the ./db folder. If you prefer to use your own GI database, simply using the makegidb.py script to create your database. Please make sure all Genbank sequendce files (.gb/.genbank) for each of the GI are located in ./db/GI_Database_seqs. Of cause you also need to locate your database description file. Just copy the format of the ./db/GI_database.txt preinstalled. Everything is self-explanary (we hope).

Usage:
python imap.py [-h] -i Input_file_path -o Output_file_path -a Algorithm_of_choice -id Identity_Threshold -cov Coverage_Threshold [-t Threads_for_BLAST] [-db Json_path] [-d Output_path]
               [-e evalue]

options:
  -h, --help            show this help message and exit
  
  -i Input_file_path    Path for input file,must be fasta format of nucleotide sequences
  
  -o Output_file_path   Path for output file
  
  -a Algorithm_of_choice
                        Choose between EMBOSS needle and blast
                        
  -id Identity_Threshold
                        Identity threshold for sequence comparison
                        
  -cov Coverage_Threshold
                        sequence coverage threshold when using blast algorithm
                        
  -t Threads_for_BLAST  Number of threads for blast function
  
  -db Json_path         Path for GI database file, in json format
  
  -d Output_path        Output directory
  
  -e evalue             evalue threshold when using blast algorithm
  

Warning:
You can choose between EMBOSS needle or blast algorithms. In theory, EMBOSS needle has better accuracies when comparing proteins. However, it is ~100 fold slower than blast. It's your call. With my test dataset, it takes ~2 minutes to calculate one GI with blast, and ~2 hours to calculate one GI with needle.


