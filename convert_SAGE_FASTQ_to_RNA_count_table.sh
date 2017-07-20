#! /bin/bash

###############################################################################
# "convert_SAGE_FASTQ_to_RNA_count_table.sh"
# created by: Jean P. Elbers
# last modified: 29 June 2017, 11:46h
###############################################################################
# Description
#
# Takes FASTQ files from the Illumina NextSeq (single end 50bp reads) created 
# by converting the .bcl files to FASTQ using Illumina BaseSpace
# and
#
# Step 1: Concatenates raw FASTQ reads from lanes 1,2,3,4 of the NextSeq
# Step 2: Demultiplexes the reads allowing for no mismatches
# Step 3: Converts FASTQ to FASTA then reverse complements
# Step 4: Get tags starting with CATG followed by 21 to 23 bases
#         (tag length=25-27bp)
# Step 5: Convert file of tags into FASTA format
# Step 6: Dereplicates FASTA file and performs usearch_global
# Step 7: Process usearch output, keeping only tags with matches
#         of 25, 26, or 27 bases to the coding RNA virtual reference of 28 bp,
#         gets rid of tags that match to more than one gene, uses the
#         "size" information from the dereplicated FASTA to add that many counts
#         for the particular tag, creates count table for each sample,
#         Note that some genes are repeated, this will be addressed in Step 8
# Step 8: Uses R and package "plyr" to process each raw count table for a
#         sample, first creates a list of inputs called "files", 
#         then reads in all raw count tables as a list, then for converts 
#         raw counts into counts by adding up counts in genes that repeat in 
#         the table, then merges all count tables- putting zeros where a gene
#         is missing in one sample, and finally saves the resulting count table
# Step 9: Remove intermediate files
# 
###############################################################################
# Prerequisites
# Remove "#" comments to copy and paste code into a Linux terminal
# A.Setup ~/bin/ directory
#    # after logging into gcflinux-vm1.pbrc.edu
#    pwd
#    # should be /home/jms/
#    mkdir bin
# B.Get and install seqtk
#    cd ~/bin/
#    mkdir seqtk
#    cd seqtk
#    wget https://raw.githubusercontent.com/lh3/seqtk/32e7903e8fd36cf8975a05295156cc69ca57c82b/Makefile
#    wget https://raw.githubusercontent.com/lh3/seqtk/32e7903e8fd36cf8975a05295156cc69ca57c82b/khash.h
#    wget https://raw.githubusercontent.com/lh3/seqtk/32e7903e8fd36cf8975a05295156cc69ca57c82b/kseq.h
#    wget https://raw.githubusercontent.com/lh3/seqtk/32e7903e8fd36cf8975a05295156cc69ca57c82b/seqtk.c
#    make
#    # commit/version is 32e7903e8fd36cf8975a05295156cc69ca57c82b
# C.Get and install sabre
#    cd ~/bin/
#    mkdir sabre
#    cd sabre/
#    wget https://raw.githubusercontent.com/najoshi/sabre/master/Makefile
#    mkdir src
#    cd src/
#    wget https://raw.githubusercontent.com/najoshi/sabre/039a55e500ba07b7e6432ea6ec2ddcfb3471d949/src/barcode.c
#    wget https://raw.githubusercontent.com/najoshi/sabre/039a55e500ba07b7e6432ea6ec2ddcfb3471d949/src/demulti_paired.c
#    wget https://raw.githubusercontent.com/najoshi/sabre/039a55e500ba07b7e6432ea6ec2ddcfb3471d949/src/demulti_single.c
#    wget https://raw.githubusercontent.com/najoshi/sabre/039a55e500ba07b7e6432ea6ec2ddcfb3471d949/src/kseq.h
#    wget https://raw.githubusercontent.com/najoshi/sabre/039a55e500ba07b7e6432ea6ec2ddcfb3471d949/src/sabre.c
#    wget https://raw.githubusercontent.com/najoshi/sabre/039a55e500ba07b7e6432ea6ec2ddcfb3471d949/src/sabre.h
#    cd ..
#    make
#    # commit/version is 039a55e500ba07b7e6432ea6ec2ddcfb3471d949
#
# D.Get and install Parallel
#    cd ~/bin/
#    wget http://ftp.gnu.org/gnu/parallel/parallel-20130922.tar.bz2
#    tar xjf parallel-20130922.tar.bz2
#    cd parallel-20130922/
#    ./configure
#    make
#    cd ..
#    rm parallel-20130922.tar.bz2
#    # version is 20130922
# E.Install R
#    sudo yum install R
#    # from within R install library("plyr")
#    install.packages("plyr")
#    # select 1 for 0-Cloud
#    # say yes to install in local package directory
#    sessionInfo()
#    #R version 3.4.0 (2017-04-21)
#    #Platform: x86_64-redhat-linux-gnu (64-bit)
#    #Running under: CentOS release 6.7 (Final)
#    #Matrix products: default
#    #BLAS: /usr/lib64/R/lib/libRblas.so
#    #LAPACK: /usr/lib64/R/lib/libRlapack.so
#    #locale:
#    # [1] LC_CTYPE=en_US.UTF-8 LC_NUMERIC=C
#    # [3] LC_TIME=en_US.UTF-8 LC_COLLATE=en_US.UTF-8
#    # [5] LC_MONETARY=en_US.UTF-8 LC_MESSAGES=en_US.UTF-8
#    # [7] LC_PAPER=en_US.UTF-8  LC_NAME=C
#    # [9] LC_ADDRESS=C LC_TELEPHONE=C
#    # [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
#    #attached base packages:
#    # [1] stats graphics grDevices utils datasets methods base
#    #
#    #other attached packages:
#    # [1] plyr_1.8.4
#    #loaded via a namespace (and not attached):
#    # [1] compiler_3.4.0 Rcpp_0.12.11
# F.Create coding and non-coding RNA references for USEARCH
#    # non-coding RNA reference
#    # create non-coding RNA reference from (Y for QNAP) Y:\REFERENCES\Nonencode references\ncRNA_mm10\VirtualReference\ncRNAmm-mRNA-CATG-28\NONCODE2016_mouse.fa\28.unique.all.virtual.tags.rec
#    # 28.unique.all.virtual.tags.rec renamed to ncRNA_28.unique.all.virtual.tags
#    ~/bin/usearch -makeudb_usearch ncRNA_28.unique.all.virtual.tags -output SAGE28_ncRNA_ref.udb
#    # coding RNA reference
#    # R:\solid\GCF\solidsageReferences\mm10\GRCm38.p1\VirtualReference\GRCm38.p1-mRNA-CATG-28\28.unique.all.virtual.tags.rec
#    ~/bin/usearch -makeudb_usearch 28.unique.all.virtual.tags.rec -output SAGE28_ref.udb
#
#    # note that you can update the virtual reference, following these steps (the old virtual reference was 
#    # mm10/GRCm38.p1, newest version is mm10/GRCm38.p5
#
#    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.25_GRCm38.p5/GCF_000001635.25_GRCm38.p5_rna.fna.gz
#    gunzip GCF_000001635.25_GRCm38.p5_rna.fna.gz
#
#    # run these 3 lines together
#    perl -pne "s/>(\w+_\d+\.\d+) (.+)/print'>gi|'; print int(rand 999999999); print'|ref|';print $1;/ge" GCF_000001635.25_GRCm38.p5_rna.fna | \
#    grep -Pv "^1" | \
#    perl -pe "s/(\w+_\d+\.\d+)/\1|/g" > GCF_000001635.25_GRCm38.p5_rna.fa
#
#    # note that creat-virtref.pl is in this github repo github.com/jelber2/sage/
#    ./create-virtref.pl
#    # question 1's answer = paste name of directory housing GCF_000001635.25_GRCm38.p5_rna.fa
#    # ex: /mnt/c/Users/ElbersJP/Desktop/old_sage_analysis_scripts/solid_sage_dianna
#    # question 2's answer = 1
#    # question 3's answer = ENTER
#    # question 4's answer = ENTER
#    # Proceed? [Y] = y
#    # this creates a new 28.unique.all.virtual.tags.rec file in the folder:
#    # ./solid_sage_diana/solid_sage_dianna-mRNA-28/GCF_000001635.25_GRCm38.p5_rna.fa/
#    # the desired input for usearch -makeudb_usearch is: 28.unique.all.virtual.tags.rec
## G.Barcode file format
#    #note salbaum-barcodes.txt is in this format (minus the # and four leading spaces)
#    CA N61M36BP70.fq
#    CT N61M36CP70.fq
#    CG N61M36DP70.fq
#    AC N58M65AP70.fq
#    AA N58M65BP70.fq
#    TC D58M04AP70.fq
#    TA D58M04BP70.fq
#    TT D58M04CP70.fq
#    GC D58M04DP70.fq
#    GA D61M47AP70.fq
#    GT D61M47BP70.fq
#    GG D61M47CP70.fq
###############################################################################
# Analysis portion of script
# Print variables
if [ $# -eq 0 ]
then
    echo ''
    echo 'You need to provide a config file'
    echo ''
    echo 'For example:'
    echo './convert_SAGE_FASTQ_to_RNA_count_table.sh config.txt'
    echo ''
    exit 1
else
    echo ''
    echo "$1 is the config file you specified"
    echo ''
    echo "It's variables are as follows:"
    echo ''
    sed -n '2p' $1
    sed -n '3p' $1
    sed -n '4p' $1
    sed -n '5p' $1
    sed -n '6p' $1
    sed -n '7p' $1
    sed -n '8p' $1
    sed -n '9p' $1
    sed -n '10p' $1
    sed -n '11p' $1
    sed -n '12p' $1
    sed -n '13p' $1
    date_config=$(sed -n '2p' $1|perl -pe 's/\t.+//g')
    analyst_config=$(sed -n '3p' $1|perl -pe 's/\t.+//g')
    fq_gz_config=$(sed -n '4p' $1|perl -pe 's/\t.+//g')
    cat_fq_config=$(sed -n '5p' $1|perl -pe 's/\t.+//g')
    barcodes_config=$(sed -n '6p' $1|perl -pe 's/\t.+//g')
    demux_config=$(sed -n '7p' $1|perl -pe 's/\t.+//g')
    fq_ext_config=$(sed -n '8p' $1|perl -pe 's/\t.+//g')
    rna_ref_config=$(sed -n '9p' $1|perl -pe 's/\t.+//g')
    ncrna_ref_config=$(sed -n '10p' $1|perl -pe 's/\t.+//g')
    rna_out_config=$(sed -n '11p' $1|perl -pe 's/\t.+//g')
    ncrna_out_config=$(sed -n '12p' $1|perl -pe 's/\t.+//g')
    usearch_config=$(sed -n '13p' $1|perl -pe 's/\t.+//g')
    if [ ! -f $barcodes_config ]
    then
        echo ''
        echo ''
        echo "$barcodes_config does not exist. Fix config file..."
        echo ''
        exit 1
    fi
    if [ ! -f $rna_ref_config ]
    then
        echo ''
        echo ''
        echo "$rna_ref_config does not exist. Fix config file..."
        echo ''
        exit 1
    fi
    if [ ! -f $ncrna_ref_config ]
    then
        echo ''
        echo ''
        echo "$ncrna_ref_config does not exist. Fix config file..."
        echo ''
        exit 1
    fi
    if [ ! -f $usearch_config ]
    then
        echo ''
        echo ''
        echo "$usearch_config not in specified folder. Fix config file..."
        echo ''
        exit 1
    fi
    if [ ! -f $barcodes_config ]
    then
        echo ''
        echo ''
        echo "$barcodes_config does not exist. Fix config file..."
        echo ''
        exit 1
    fi
    echo ''
    echo ''
    echo 'Are you satisfied with the config file settings?'
    echo 'y or n'
    echo ''
    while read -r question1
    do
        if [ $question1 = "y" ]
        then
            echo ''
            echo ''
            echo 'Do you wish to demultiplex the data?'
            echo 'y or n'
            echo ''
            while read -r question2
            do
                if [ $question2 = "y" ]
                then
                    echo ''
                    echo ''
                    #Step 1
                    if [ -f $cat_fq_config ]
                    then
                        rm $cat_fq_config
                    fi
                    echo 'Starting Step 1: Concatenate reads'
                    cat *$fq_gz_config > $cat_fq_config
                    echo "Done with Step 1"
                    echo ''
                    echo ''
                    #Step 2
                    echo 'Starting Step 2: Demultiplex reads'
                    ~/bin/sabre/sabre se -m 0 -f $cat_fq_config -b $barcodes_config -u /dev/null > $demux_config
                    cat $demux_config
                    echo "Done with Step 2"
                    echo ''
                    echo ''
                    break
                elif [ $question2 = "n" ]
                then
                    echo ''
                    echo ''
                    break
                else
                    echo ''
                    echo ''
                    echo "You didn't enter y or n"
                    echo ''
                    echo "Please enter y or n"
                fi
            done
            #Step 3
            echo 'Starting Step 3: FASTQ -> FASTA -> reverse complement'
            ls *.fq | perl -pe "s/.fq//g" > samples
            ~/bin/parallel-20130922/src/parallel '~/bin/seqtk/seqtk seq -a {} | ~/bin/seqtk/seqtk seq -r > {.}.rc' ::: *$fq_ext_config
            echo "Done with Step 3"
            echo ''
            echo ''
            #Step 4
            echo 'Starting Step 4: Get tags starting with CATG followed by 21 to 23 bases'
           ~/bin/parallel-20130922/src/parallel 'grep -Po "CATG\w{21,23}" {} > {.}.fa.txt && rm {}' ::: *.rc
            echo "Done with Step 4"
            echo ''
            echo ''
            #Step 5
            echo 'Starting Step 5: Convert file of tags into FASTA format'
            while read i;do
                export LENGTH="$(wc -l < $i.fa.txt)"
                seq 1 $LENGTH | perl -pe "s/(\d+)/>\1/" > bar
                awk '{print $0 "\n"}' $i.fa.txt |sed '1s/^/\n/' > test
                paste -d '\n' bar <(awk 'NR%2==0' test) > tmp && mv tmp $i.fa.txt.fa && rm bar test && rm $i.fa.txt
            done < samples
            echo "Done with Step 5"
            echo ''
            echo ''
            #Step 6
            echo 'Starting Step 6: Dereplicates FASTA file and performs usearch_global on coding and non-coding RNA refereneces'
            while read i;do
                $usearch_config -fastx_uniques $i.fa.txt.fa -fastaout $i.derep.fa -sizeout
                $usearch_config -usearch_global $i.derep.fa -db $rna_ref_config -id 1.0 -strand both -threads 32 -blast6out $i.usearch.out
                $usearch_config -usearch_global $i.derep.fa -db $ncrna_ref_config -id 1.0 -strand both -threads 32 -blast6out $i.usearch.nc.out
                echo "Done with sample $i"
                echo ''
                echo ''
                echo ''
            done < samples
            echo "Done with Step 6"
            echo ''
            echo ''
            #Step 7
            echo 'Starting Step 7:Process usearch output'
            echo ''
            echo 'Processing coding RNA'
            while read i;do
                awk -F"\t" '$4 == "25"||$4 == "26"||$4 == "27" {print $1"\t"$2}' $i.usearch.out | \
                perl -pe "s/\w+ (.+)/\1/" | \
                perl -pe "s/( )/\/\/\//g" | \
                awk '{n = split($2, t, "///"); _2 = x
                split(x, _)
                for (i = 0; ++i <= n;)
                _[t[i]]++ || _2 = _2 ? _2 "," t[i] : t[i]
                $2 = _2
                }-2' | grep -v "," | perl -pe "s/.+;size=(\d+);/\1/g" |sort -k 2,2 | sed "1i $i\tGene" > $i.usearch.out.count
            done < samples
            echo 'Processing non-coding RNA'
            while read i;do
                awk -F"\t" '$6 == "25"||$6 == "26"||$6 == "27" {print $1"\t"$4}' $i.usearch.nc.out | \
                perl -pe "s/( )/\/\/\//g" | \
                awk '{n = split($2, t, "///"); _2 = x
                split(x, _)
                for (i = 0; ++i <= n;)
                _[t[i]]++ || _2 = _2 ? _2 "," t[i] : t[i]
                $2 = _2
                }-2' | grep -v "," | perl -pe "s/.+;size=(\d+);/\1/g" |sort -k 2,2 | sed "1i $i\tGene" > $i.usearch.nc.out.count
            done < samples
            echo "Done with Step 7"
            echo ''
            #Step 8
            echo 'Starting Step 8a:Make coding RNA count table'
            echo "library(plyr)" > make.count.table.R
            echo "files <-list.files(pattern='usearch.out.count')" >> make.count.table.R
            echo "rawcounts <- list()" >> make.count.table.R
            echo "counts <- list()" >> make.count.table.R
            echo "for (i in 1:length(files)) {rawcounts[[i]] <- read.table(files[i],header=T)}" >> make.count.table.R
            echo "for (i in 1:length(files)) {counts[[i]] <- ddply(data.frame(rawcounts[i]),'Gene',numcolwise(sum))}" >> make.count.table.R
            echo "merged.data.frame = Reduce(function(...) merge(..., all=T), counts)" >> make.count.table.R
            echo "write.table(merged.data.frame, file='$rna_out_config',quote=F,sep='\t',eol='\r\n',na='0',row.names=F,col.names=T)" >> make.count.table.R
            Rscript make.count.table.R
            echo ''
            echo ''
            echo 'Starting Step 8b:Make non-coding RNA count table'
            echo "library(plyr)" > make.nc.count.table.R
            echo "files <-list.files(pattern='.usearch.nc.out.count')" >> make.nc.count.table.R
            echo "rawcounts <- list()" >> make.nc.count.table.R
            echo "counts <- list()" >> make.nc.count.table.R
            echo "for (i in 1:length(files)) {rawcounts[[i]] <- read.table(files[i],header=T)}" >> make.nc.count.table.R
            echo "for (i in 1:length(files)) {counts[[i]] <- ddply(data.frame(rawcounts[i]),'Gene',numcolwise(sum))}" >> make.nc.count.table.R
            echo "merged.data.frame = Reduce(function(...) merge(..., all=T), counts)" >> make.nc.count.table.R
            echo "write.table(merged.data.frame, file='$ncrna_out_config',quote=F,sep='\t',eol='\r\n',na='0',row.names=F,col.names=T)" >> make.nc.count.table.R
            Rscript make.nc.count.table.R
            echo "Done with Step 8"
            echo ''
            echo ''
            #Step 9
            echo "Step 9: Remove intermediate files"
            rm $cat_fq_config
            while read i
            do
                rm $i$fq_ext_config
                rm $i.fa.txt.fa
                rm $i.derep.fa
                rm $i.usearch.out
                rm $i.usearch.out.count
                rm $i.usearch.nc.out
                rm $i.usearch.nc.out.count
            done < samples
            rm make.count.table.R
            rm make.nc.count.table.R
            rm samples
            echo ''
            echo 'Done with all steps'
            echo ''
            exit 1
        elif [ $question1 = "n" ]
        then
            echo ''
            echo 'Please edit config file'
            echo ''
            exit 1
        else
            echo ''
            echo "You didn't enter y or n"
            echo ''
            echo "Please enter y or n"
        fi
    done
fi
