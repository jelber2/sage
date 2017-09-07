# sage
## script for converting Sage FASTQ files into coding and non-coding RNA count tables
### To download enter the following command in a terminal window
    git clone https://github.com/jelber2/sage.git
    # this may take more than 10 minutes
    # this will create a folder called sage in the current directory
    # you should see something like this in the terminal:
    Cloning into 'sage'...
    remote: Counting objects: 52, done.
    remote: Compressing objects: 100% (3/3), done.
    remote: Total 52 (delta 0), reused 1 (delta 0), pack-reused 49
    Unpacking objects: 100% (52/52), done.
    Checking connectivity... done.
## NOTE YOU MUST UNZIP THE SAGE REFS
    cd /home/jms/sage/
    unzip SAGE28_ncRNA_ensembl_ref.zip
    unzip SAGE28_ncRNA_ncbi_ref.zip
    unzip SAGE28_ncRNA_nonencode_2016_ref.zip
    unzip SAGE28_ref_GCF_000001635.25_GRCm38.p5_NCBI_2017_07_20.zip
## You must also delete the example FASTQ.GZ files if you copied the FASTQ files for analysis on the current experiment
    rm ?.fastq.gz
### To update to the newest version,
### go to the directory housing the script
    # example:
    cd /home/jms/sage/
### Then type the following command
    git pull origin master
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
    SAGE28_ref_GCF_000001635.25_GRCm38.p5_NCBI_2017_07_20.udb                        # location of coding RNA reference
    SAGE28_ncRNA_nonencode_2016_ref.udb		# location of nonencode non-coding RNA reference
    SAGE28_ncRNA_ncbi_ref.udb		# location of ncbi non-coding RNA reference
    SAGE28_ncRNA_ensembl_ref.udb		# location of ensembl non-coding RNA reference
    coding-rna-count-table.txt      # name of coding RNA output count table
    nonencode-ncrna-count-table.txt	# name of nonencode non-coding RNA output count table
    ncbi-ncrna-count-table.txt	# name of ncbi non-coding RNA output count table
    ensembl-ncrna-count-table.txt	# name of ensembl non-coding RNA output count table
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
### Output will be the following count tables:
    coding-rna-count-table.txt
    nonencode-ncrna-count-table.txt
    ncbi-ncrna-count-table.txt
    ensembl-ncrna-count-table.txt
### Output will also be:
    demux.log
    undemultiplexed.fastq.gz
    sample1.fq.gz
    sample2.fq.gz
    sampleN.fq.gz
