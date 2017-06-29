# sage
## script for converting Sage FASTQ files into coding and non-coding RNA count tables
### To download the newest version enter the following command in a terminal window
    git clone https://github.com/jelber2/sage.git
### To use, invoke:
    ./convert_SAGE_FASTQ_to_RNA_count_table.sh 
### If you get an error such as:
    bash: ./convert_SAGE_FASTQ_to_RNA_count_table.sh: permission denied
### Then enter
    chmod u+x ./convert_SAGE_FASTQ_to_RNA_count_table.sh
### Now run the program again
    ./convert_SAGE_FASTQ_to_RNA_count_table.sh
### You should get the following:
    
    You need to provide a config file

    For example:
    ./convert_SAGE_FASTQ_to_RNA_count_table.sh config.txt

### If you provide a config file, you should get the following output
    
    config.txt is the config file you want to use
    
    22 June 2017                    # Date
    Jean P. Elbers                  # Analyst
    .fastq.gz                       # extension of FASTQ files to demultiplex
    sage.fastq.gz                   # name of concatenated FASTQ file for Step 1: Concatenates raw FASTQ
    gettys-barcodes.txt             # name and path to the barcode file for Step 2: Demultiplexes
    demux.log                       # name for demultiplexing log for Step 2: Demultiplexes
    .fq                             # extension that you gave the demultiplexed fastq files from in Step 2
    /home/jms/References_SAGE/SAGE28_ref.udb                        # location of coding RNA reference
    /home/jms/References_SAGE/SAGE28_ncRNA_ref.udb          # location of non-coding RNA reference
    coding-rna-count-table.txt      # name of coding RNA output count table
    non-coding-rna-count-table.txt  # name of non-coding RNA output count table
    /usr/local/bin/usearch          # location of usearch executable

    AND ANY ERRORS WOULD APPEAR HERE

### If no errors, then you will get the prompt:
    Are you satisfied with the config file settings?
    y or n
### If you enter "n", then you will get
    Please edit config file
### If you enter "y", then you will get
    Do you want to demultiplex the reads?
    y or n
### If "y", then you will get
    Starting Step 1: Concatenate reads
### If "n", then you will get
    Starting Step 3: FASTQ -> FASTA -> reverse complement
