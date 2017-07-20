#!/usr/bin/perl -w

#############################################################################
# Create a solid-sage virtual reference 
#
# There can be multiple reference files that are processed together.
# For example, each chromosome could be in its own separate file.
#
#############################################################################
# 
# Original code base from http://solidsoftwaretools.com/gf/project/solid-sage/
#Data analysis GUI for SOLiD SAGE
#Version 1.1.0
#03/2010
#Developed by Xiequn (Tony) Xu
#
# Modified by Diana Holmes 10/2011
# 
# Strip out the creation of the virtual reference file from the solid-sage
# program so that it can be run separately.
#############################################################################
use strict;
use warnings;

use File::Basename;
use Getopt::Std;

# This helps us find things based on where the script is run from
# It sets the lib path to include the place where we were called
# from and from there, we can construct the path to the auxiliary
# programs we use.
use FindBin;
use lib "$FindBin::Bin";

#############################################################################
# globals are stored in a common file so that all the auxiliary programs
# can use the same values.
#############################################################################
use solidsage;
use vars qw(%c);
*c=\%solidsage::c;

#############################################################################
# debugging leaves intermediate files around
# 0 to delete all intermediate files
# 1 to leave them around
#############################################################################
my $debugging = 1;


#############################################################################
#                              
# MAIN                         
#                              
#############################################################################

my $version = 'v0.1';

# Get the name of this script as it was run
my $progname = fileparse($0);

# declare the perl command line flags/options we want to allow
# 
# by default, we assume we are in interactive mode and we will prompt for 
# the information we need
#
# When in interactive mode, all commandline options are ignored.
#
# -s: silent mode, uses built-in defaults plus commandline options
#
# -r ReferenceDirectory
# 
# -t ReferenceType (mRNA, aRNA, gDNA)
my $ReferenceDirInput;
my $ReferenceTypeInput;
my $InteractiveMode;

my %options=();
getopts("r:t:sh", \%options);

if (defined $options{h}) {
	# Print usage and exit
	Usage();
	exit(0);
}

# interactive mode
# 0 = batch mode
# 1 = interactive mode
if (defined $options{s}) {

	# batch mode - everything is specified with commandline options
	$InteractiveMode = 0;

	# Directory containing the reference file to use for
	# creating the virtual tag reference file.
	if (defined $options{r}) {
		$ReferenceDirInput = $options{r};

		if (! -d $ReferenceDirInput) {
			print STDERR "$ReferenceDirInput does not exist\n";
			Usage();
			exit 1;
		}
	}
	else {
		print STDERR "Must include -r to specify reference directory\n";
		Usage();
		exit 1;
	}

	# Reference type
	if (defined $options{t}) {
		$ReferenceTypeInput = $options{t};

		# TBD - need to do case insensitive comparisons here
		if (($ReferenceTypeInput ne $c{mRNA}) && ($ReferenceTypeInput ne $c{aRNA}) && ($ReferenceTypeInput ne $c{gDNA})) {
			print STDERR "Invalid reference type $ReferenceTypeInput\nValid values are $c{mRNA}, $c{aRNA}, $c{gDNA}\n";
			Usage();
			exit 1;
		}
	}
	else {
		print STDERR "Must include -t to specify reference type\n";
		Usage();
		exit 1;
	}
}
else {
	$InteractiveMode = 1;
}

# Set up the log file
my $runLog = "$progname.$c{runLogDefault}";
open RUNLOG, ">$runLog" or die "Can't open $runLog: $!\n";

# Capture all the stderr output - this is mainly used to get the
# output from the auxiliary programs we use.
my $errorLog = "$progname.$c{errorLogDefault}";
open STDERR, ">$errorLog" or die "Can't redirect errors to $errorLog $!\n";

my $defaultSite = $c{defaultSite};
my $defaultVirtualTagLength = $c{defaultVirtualTagLength};

my $setReferenceDir;
my $setReferenceType;
my $setVirtualTagLength;
my $setSite;
my $setOutputDir;

if ($InteractiveMode) {

	# Interactive Mode - get all arguments from the user

	if (defined $ReferenceDirInput) {
		$setReferenceDir = $ReferenceDirInput;
	}
	else {
		$setReferenceDir = $c{defaultReferenceDir};
	}

	# Get the name of the directory with the reference sequences
	$setReferenceDir = promptUser("Select directory containing the reference sequences", $setReferenceDir);
	if (!defined $setReferenceDir) {
		exit(1);
	}

	# and the reference sequence type
	my $prompt = "";
	my @choices = ($c{mRNAChoice}, $c{aRNAChoice}, $c{gDNAChoice});
	foreach my $i (0 .. $#choices) {
		$prompt = $prompt . ($i + 1) . ": " . $choices[$i] . "\n";
	}
	my $rtv = promptUser("\n" . $prompt . "Choose reference type", "1");
	if ($rtv == 3) {
		$setReferenceType = $c{gDNA};
	} elsif ($rtv == 2) {
		$setReferenceType = $c{aRNA};
	} else {
		$setReferenceType = $c{mRNA};
	}
	if (!defined $setReferenceType) {
		exit(1);
	}

	# get virtual tag length
	$setVirtualTagLength = promptUser("Enter virtual tag length", $defaultVirtualTagLength);
	if (!defined $setVirtualTagLength) {
		exit(1);
	}

	# get digestion site
	$setSite = promptUser("Enter digestion site", $defaultSite);
	if (!defined $setSite) {
		exit(1);
	}
}
else {
	# Use defaults for the parameters

	if (defined $ReferenceDirInput) {
		$setReferenceDir = $ReferenceDirInput;
	}
	else {
		$setReferenceDir = $c{defaultReferenceDir};
	}
	$setReferenceType = $ReferenceTypeInput;
	$setVirtualTagLength = $defaultVirtualTagLength;
	$setSite = $defaultSite;
}

# validate the reference file(s)
if (!ValidateReferenceFiles($setReferenceDir, $setReferenceType)) {
	# error was already displayed
	exit(1);
}

# for the output directory, use the last part of the reference directory name
# and add in the sequence type and tag length, all separated by dashes
$setOutputDir = fileparse($setReferenceDir) . "-" . $setReferenceType . "-" . $setVirtualTagLength;

if (-d $setOutputDir) {
	ErrorMsg("$setOutputDir already exists. Delete or rename and try again\n");		
	exit(1);
}

if ( ! mkdir($setOutputDir)) {
	ErrorMsg("Unable to create $setOutputDir - $!\n");
	exit(1);
}

# Display the parameters we're going to use
my $mapMsg = "\n\nParameters:\n";
$mapMsg .= "\tReference folder\t\t$setReferenceDir\n";
$mapMsg .= "\tReference type\t\t\t$setReferenceType\n";
$mapMsg .= "\tVirtual reference size\t\t$setVirtualTagLength\n";
$mapMsg .= "\tDigestion site\t\t\t$setSite\n";
$mapMsg .= "\tOutput Directory\t\t$setOutputDir\n";

if ($InteractiveMode) {
	# Confirm that the user wants to continue
	print $mapMsg;
	my $proc = promptUser("Proceed?", "Y");
	if (($proc ne "y") && ($proc ne "Y")) {
		exit;
	}
}

LogMsg($mapMsg);

# Do it!
CreateVirtualReference();


#############################################################################
#
# SUBROUTINES
#
#############################################################################

#############################################################
#############################################################
sub Usage {
	print "Usage: \n";
	print "	$progname\n";
	print "		You will be prompted for all parameters\n";
	print "\n";
	print "	$progname -s -r ReferenceDirectory -t ReferenceType\n";
	print "		Uses built-in defaults plus commandline options\n";
	print "\n";
	print "			ReferenceDirectory contains the .fa files\n";
	print "			ReferenceType is mRNA, aRNA, or gDNA\n";
	print "			Output directory will use the last part of the reference directory name\n";
	print "				and add in the sequence type and tag length, all separated by dashes.\n";
	print "			OutputDirectory will contain all the virtual reference files.\n";
}

#############################################################
#############################################################
sub CreateVirtualReference {

	# Get the list of all the .fa files from the directory
	my @refFileList = readInRefDir($setReferenceDir);
	if (scalar @refFileList == 0) {
		die "Error: No fasta format file found in the designated directory\n$setReferenceDir\nPlease double check and be sure the names of all reference fasta files end with .fa\n";
	}

	# set the reference type for all of the files we found
	my %refFiles;
	foreach my $f (@refFileList) {
		$refFiles{$f} = $setReferenceType;
	}

	#############################################################
	# INTERNALFILE: store all the giID information
	# This is a tab delimited file that contains the gi ID's from all of the reference files
	# Column 1: gi ID
	# Column 2: full description line from the reference file for this gi ID
	# Column 3: NM (accession) number
	# Column 4: gene description
	# Column 5: gene name from description (multiple genes will be separated by space)
	#
	# getVirtualTagsOfAFile writes to the giDef file as it looks for the virtual tags in the reference file
	#
	# NOTE: There is no checking for duplicates in the GIDEF file so hopefully
	# all the reference files use unique giID's 
	#############################################################
	my $giDefFile = $setOutputDir . "/" . $c{giDefFileDefault};
	LogMsg( "Opening $giDefFile\n");
	open (GIDEF, ">$giDefFile") || die "Can't open $giDefFile $!\n";

	# go through each reference file and create the virtual reference
	foreach my $f (sort {$b cmp $a} keys %refFiles) {

		LogMsg( "Processing reference file $f\n");

		my $referenceType = $refFiles{$f};

		# the virtual reference output files for this reference file will go 
		# in a directory that is the same as the last part of the reference
		# file path
		my $refBasename = $setOutputDir . "/" . fileparse($f);
		chomp($refBasename);

		# Create the directory that will hold the results. 
		# The name of the directory matches the name of the reference file
		if (-d $refBasename) {
			ErrorMsg("$refBasename already exists\n");
			exit(1);
		}
		my $refBasenameMakeDir = system("mkdir $refBasename");
		if ($refBasenameMakeDir ne 0) {
			ErrorMsg("mkdir $refBasename failed $!\n");
			exit(1);
		}

		my $VirtTagTabFileDefault = "$setVirtualTagLength.all.virtual.tags.tab";
		my $VirtTagFAFileDefault = "$setVirtualTagLength.unique.all.virtual.tags.fa";
		my $VirtTagRECFileDefault = "$setVirtualTagLength.unique.all.virtual.tags.rec";

		my $VirtTagTabFile;
		my $VirtTagFAFile;
		my $VirtTagRECFile;

		# Get all the virtual tags from the sequence file
		# side-effect: writes to GIDEF file
		# returns array indexed by virtual tag that contains the list of places this tag occurs (giID_index)
		# returns array indexed by gi id that contains the corresponding gene for that gi
		my ($giIDLocationsByTag, $geneMappingsByTag) = getVirtualTagsOfAFile($f, $referenceType, $setSite, $setVirtualTagLength);
		my %giIDLocationsByTag = %$giIDLocationsByTag;
		my %geneMappingsByTag = %$geneMappingsByTag;

		if ($debugging == 1) {
			my $numkeys = keys %giIDLocationsByTag;
			LogMsg( "numkeys for giIDLocationsByTag=$numkeys\n");

			$numkeys = keys %geneMappingsByTag;
			LogMsg( "numkeys for geneMappingsByTag=$numkeys\n");
		}

		#############################################################
		# INTERNALFILE: all.virtual.tags.tab file is tab-delimited
		# Column 1: tag sequence
		# Column 2: a list of places where this tag occurs, each separated by a space. 
		#           The place is specified as giID_indexInReferenceSequence
		# Column 3: a list of genes that corresponds to the giIDs in column 2
		#############################################################
		$VirtTagTabFile = "$refBasename/$VirtTagTabFileDefault";
		LogMsg( "Creating $VirtTagTabFile\n");

		open (AVT, ">$VirtTagTabFile") || die "Can't open $VirtTagTabFile: $!\n";
		foreach my $tag (keys %giIDLocationsByTag) {
			print AVT "$tag\t@{ $giIDLocationsByTag{$tag} }\t@{ $geneMappingsByTag{$tag} }\n";
		}
		close (AVT) || die "Can't close $VirtTagTabFile: $!\n";
		%giIDLocationsByTag = ();
		%geneMappingsByTag = ();

		# read from all.virtual.tags.tab
		# returns array indexed by tag sequence that contains the places the tag occurs
		# NOTE: Not sure why this is done here since we just created the file and
		# already had these arrays but it was in the original code for some reason....
		($giIDLocationsByTag, $geneMappingsByTag) = getUniqueVirtualTags($VirtTagTabFile);
		%giIDLocationsByTag = %$giIDLocationsByTag;
		%geneMappingsByTag = %$geneMappingsByTag;

		if ($debugging == 1) {
			my $numkeys2 = keys %giIDLocationsByTag;
			LogMsg( "numkeys for giIDLocationsByTag=$numkeys2\n");

			$numkeys2 = keys %geneMappingsByTag;
			LogMsg( "numkeys for geneMappingsByTag=$numkeys2\n");
		}

		#############################################################
		# INTERNALFILE: unique.all.virtual.tags.fa contains two lines
		# first line: >conc_seq
		# second line: All the tags concatenated together with 'NNN' spacer between tags
		#############################################################
		my $spacer = $c{spacer};
		$VirtTagFAFile = "$refBasename/$VirtTagFAFileDefault";

		#############################################################
		# INTERNALFILE: unique.all.virtual.tags.rec contains two lines for each tag
		# first line: >sequence_number TAB locations TAB genes
		# 	where locations is places where tag occurs, each separated by a space 
		# 			(the place is specified as giID_indexInReferenceSequence
		#         genes is genes that correspond to the giID
		# second line: tag sequence
		#############################################################
		$VirtTagRECFile = "$refBasename/$VirtTagRECFileDefault";

		# Create the contents for unique.all.virtual.tags.fa and unique.all.virtual.tags.rec
		my $seqCt = 1;
		LogMsg( "Creating $VirtTagFAFile and $VirtTagRECFile\n");
		open (FAFILE, ">$VirtTagFAFile") || die "Can't open $VirtTagFAFile: $!\n";
		open (RECFILE, ">$VirtTagRECFile") || die "Can't open $VirtTagRECFile: $!\n";
		print FAFILE ">conc_seq\n";
		foreach my $tag (keys %giIDLocationsByTag) {
			# write out each tag sequence with the spacer in between
			print FAFILE $tag . $spacer;

			# if all the gene names for this tag are not the same, we mark
			# this with an "I" at the end of the sequence number
			my $currentGene;
			my $Ignore = ""; # assume we don't ignore this
			my @genes = split(/ /, $geneMappingsByTag{$tag}[0]);
			if (scalar @genes > 1) {
				my $current = $genes[0];
				foreach my $i (1 .. $#genes) {
					if ($current ne $genes[$i]) {
						$Ignore = "I";
						#LogMsg("Ignoring $seqCt\n");
						last;
					}
				}
			}

			# write out sequence number followed by tag locations
			# next line is the tag sequence
			print RECFILE ">$Ignore$seqCt\t@{ $giIDLocationsByTag{$tag} }\t@{ $geneMappingsByTag{$tag} }\n$tag\n";

			$seqCt++;
		}
		print FAFILE "\n";
		%giIDLocationsByTag = ();
		%geneMappingsByTag = ();
		close (FAFILE) || die "Can't close $VirtTagFAFile: $!\n";
		close (RECFILE) || die "Can't close $VirtTagRECFile: $!\n";

		LogMsg( "Done Processing reference file $f\n");
	}

	close (GIDEF) || die "Can't close $giDefFile $!\n";

	my $donemsg="\n************ALL DONE!!!************\n\nResults are stored in current directory.\nLog file is $runLog and error file is $errorLog\n";
	LogMsg( $donemsg);

	close (RUNLOG);
}

####################################
####################################
sub readInRefDir {
	my $rrd = shift;
	my @rf;
	LogMsg( "Reading $rrd\n");
	if (opendir (RRD, $rrd)) {
		@rf = grep { /\.fa$/ } readdir RRD;
		if (closedir (RRD)) {
		} else {
			ErrorMsg(" Can't close directory $rrd: $!");
		}
	} else {
		ErrorMsg(" Can't open directory $rrd: $!");
	}

	foreach my $f (@rf) {
		$f = $rrd . '/' . $f;
	}
	return @rf;
}

####################################
####################################
sub revcom{
	my $original = shift;
	my $reverse = reverse($original);
	$reverse =~ tr/ATCGatcg/TAGCtagc/;
	return $reverse;
}

#############################################################
# Input:
# Name of all.virtual.tags.tab file
#
# File format:
# Column 1: tag sequence
# Column 2: a list of places where this tag occurs, each separated by a space. 
#           The place is specified as giID_indexInReferenceSequence
# Column 3: a list of genes that corresponds to the giIDs in column 2
#
# Output:
# Array indexed by tag sequence that contains the places the tag occurs
# Array indexed by tag sequence that contains the genes that correspond to the giIDs
#############################################################
sub getUniqueVirtualTags {
	my ($VirtTagTabFile) = @_;

	my %LocationsByTag;
	my %GenesByTag;

	LogMsg( "Reading in $VirtTagTabFile\n");

	if (open (AF, $VirtTagTabFile)) {
		while (<AF>) {
			chomp;

			# split at the TAB
			# temp[0] is tag sequence
			# temp[1] is list of places
			# temp[2] is list of genes
			my @temp = split(/\t/, $_);

			push @{ $LocationsByTag{$temp[0]} }, $temp[1];
			push @{ $GenesByTag{$temp[0]} }, $temp[2];
		}
		if (close (AF)) {
		} else {
			ErrorMsg(" Can't close $VirtTagTabFile $!\n");
		}
	} else {
		ErrorMsg(" Can't open $VirtTagTabFile $!\n");
	}
	return (\%LocationsByTag, \%GenesByTag);
}

#############################################################
# Input parameters:
# reference file name
# reference type (gDNA, mRNA or aRNA)
# digestion site
# virtual tag length
#
# Output:
# Array indexed by virtual tag that contains the list of places this tag occurs (giID_index)
# Array indexed by gi ID that contains the gene that corresponds to this gi ID
#
# Side effect:
# Writes to the tags definition file (GIDEF)
#############################################################

sub getVirtualTagsOfAFile {
	my ($refFile, $referenceType, $site, $virtualTagLength) = @_;

	LogMsg( "Creating Virtual Tags file\n");

	my %giIDLocationsByTag;	# gi ID and position indexed by virtual tag
	my %geneMappingsByTag;	# genes that correspond to the gi ID indexed by virtual tag

	my %geneByGIid;		# array of genes indexed by gi ID
	my %refFileSeqByGIid; 	# sequences from the reference file indexed by gi ID

	my $fileType;
	my $fileTypeRefSeq = 'r';
	my $fileTypeGenomic = 'g';

	# Open the reference file
	LogMsg( "Reading in $refFile and writing to GIDEF file\n");
	if (open (F, $refFile)) {
		my $giID;

		# read each line of the file
		while (<F>) {
			chomp;

			# If we haven't determined the file type (RefSeq or Genomic), then do that now
			# TBD turns out that the NCBI genomic reference uses >gi in the description line
			# TBD so that isn't really a good way to be able to tell the difference. May have to
			# TBD ask the user and do it that way to ensure it's handled correctly.
			if (!defined ($fileType)) {
				if (/^>gi/) {
					# RefSeq file
					$fileType = $fileTypeRefSeq;
					LogMsg("Determined $refFile is RefSeq reference\n");
				}
				else {
					if (/^>chr/) {
						# Genomic file
						$fileType = $fileTypeGenomic;
						LogMsg("Determined $refFile is genomic reference\n");
					}
					else {
						die "Unrecognized reference file type\n";
					}
				}
			}

			# Process the file according to what type it is
			if ($fileType eq $fileTypeRefSeq) {
				# This is a RefSeq file

				# Is this a description line or a sequence line?
				if (/^>/) {
					# It's a description line that is delimited by the bar character
					# The description is formatted as follows:
					# temp[0] = "gi"
					# temp[1] = gi number
					# temp[2] = source (usually "ref")
					# temp[3] = NM (accession) number
					# temp[4] = description (which may contain commas)


					# Split up the description line at the bar character
					my @temp = split(/\|/, $_);

					# The description field might contain a gene name in parentheses.
					# There may be more than one set of parentheses.
					# Assume the gene name is the last set.
					my $desc = $temp[4];
					my $gene;
					my $lastopenparen = rindex($desc, "(");
					my $lastcloseparen = rindex($desc, ")");
					if ($lastcloseparen < 0) {
						# no gene name found
						$gene = "unknown";
					}
					else {
						$gene = substr($desc, $lastopenparen+1, $lastcloseparen - $lastopenparen - 1);
					}

					# Write a line to the tags definition file which is tab delimited
					# Column 1: gi ID
					# Column 2: full line from reference file
					# Column 3: NM (accession) number
					# Column 4: description
					# Column 5: gene name from description
					#
					# NOTE: This contains repeated information because the original code
					# didn't break up the description any further. New functionality which
					# needs the separate pieces was added but to make sure the old
					# code doesn't break, we won't change the order of things and we'll
					# have some duplicate information
					$giID = $temp[1];
					print GIDEF "$giID\t$_\t$temp[3]\t$temp[4]\t$gene\n";

					# Save the gene in the geneByGIid (indexed by gi ID)
					$geneByGIid{$giID} = $gene;
				} 
				else {
					# Save the sequence in the refFileSeqByGIid array (indexed by the gi ID)
					$refFileSeqByGIid{$giID} .= $_;
				}
			} # End RefSeq handling

			if ($fileType eq $fileTypeGenomic) {
				# This is a Genomic reference
				# That means that there is one line with the description and
				# then the rest of the file is the sequence.
				# The description may or may not contain useful information
				# (we assume it only contains the chromosome name).

				# Is this a description line or a sequence line?
				if (/^>/) {

					# The GIDEF files needs to be created as we find the virtual tags
					# so that each tag is considered a unique entity
					# For now, we'll just add a line to the GIDEF file and read
					# in the sequence but we'll add the real GIDEF contents later

					# It's a description line that specifies what chromosome this is
					# We use the chromosome name as the giID
					$giID = substr $_, 1;

					# Write a line to the tags definition file which is tab delimited
					print GIDEF "$giID\t$_\t$giID\t$giID\t$giID\n";

					# Save the chr in the geneByGIid (indexed by gi ID)
					# Saving just the chr isn't useful for genomic. We really
					# want to know the position in the chr sequence.
					# We'll add that further below but for now, just store the
					# chr here.
					$geneByGIid{$giID} = $giID;
				} 
				else {
					# Save the sequence in the refFileSeqByGIid array (indexed by the gi ID)
					$refFileSeqByGIid{$giID} .= $_;
				}
			}
		}
		close (F) || die "Can't close $refFile: $!\n";
	} else {
		ErrorMsg(" Can't open $refFile at getVirtualTagsOfAFile: $!\n");
	}

	if ($debugging) {
		my $numkeys = keys %geneByGIid;
		LogMsg("GetVirtualTagsOfAFile: geneByGIid numkeys=$numkeys\n");

		$numkeys = keys %refFileSeqByGIid;
		LogMsg("GetVirtualTagsOfAFile: refFileSeqByGIid numkeys=$numkeys\n");
	}

	# Go through the array of sequences looking for the digestion site in each one
	REF:foreach my $giID (keys %refFileSeqByGIid) {
		my $seq = $refFileSeqByGIid{$giID};
	  	my $index = 0;
		my $last_index = 0;

		INDEX:while ($index != -1) {
			# Look for the next instance of the digestion site in the sequence
			$index = index($seq, $site, $last_index + 1);
			if ($index == -1) {
				# Site not found in this sequence. Go to the next sequence.
				next REF;
			}

			# The site was found. Pull out the virtual tag starting at 
			# the digestion site and continuing for the length of the virtual tag
			my $vtag;
			$vtag = substr($seq, $index, $virtualTagLength);

			# Make sure there aren't any N's in the tag and make sure there aren't
			# any repeat sections (designated with lowercase letters)
			if ($vtag =~ /[Ncatg]/ ) {
				$last_index = $index;
				next INDEX;
			}

			# Add the tag to the array according to what type of reference we're using
			# NOTE: There can be multiple giID.position entries for a given tag
			if ($referenceType eq $c{mRNA} || $referenceType eq $c{gDNA}) {
				# sense strand
				if (length($vtag) == $virtualTagLength) {

					# RefSeq
					if ($fileType eq $fileTypeRefSeq) {
						# Save the location of this site for this tag
						push @{ $giIDLocationsByTag{$vtag} }, $giID . '_' . $index ;

						# Add the gene to the list for this tag
						push @{ $geneMappingsByTag{$vtag} }, $geneByGIid{$giID};
			      		}

					# Genomic
					if ($fileType eq $fileTypeGenomic) {
						# The giID will be the chromosome name (which is currently in $giID)
						# followed by a period and then the position of this site
						# i.e. chr1.position
						my $giIDgenomic = $giID . '.' . $index;

						# We need to add a line to GIDEF for this site
						# We just set all fields to the same value
						print GIDEF "$giIDgenomic\t$giIDgenomic\t$giIDgenomic\t$giIDgenomic\t$giIDgenomic\n";
						# Save the location of this site for this tag	
						push @{ $giIDLocationsByTag{$vtag} }, $giIDgenomic . '_' . $index ;

						# Push the same thing that we pushed for the location
						push @{ $geneMappingsByTag{$vtag} }, $giIDgenomic . '_' . $index ;
					}
				}
			}

			if ($referenceType eq $c{gDNA} || $referenceType eq $c{aRNA}) {
				# antisense strand
				$vtag = revcom(substr($seq, $index - $virtualTagLength + length($site), $virtualTagLength));
				if (length($vtag) == $virtualTagLength) {

					# RefSeq
					if ($fileType eq $fileTypeRefSeq) {
						push @{ $giIDLocationsByTag{$vtag} }, $giID . '_' . '-' . $index;

						# Add the gene to the list for this tag
						push @{ $geneMappingsByTag{$vtag} }, $geneByGIid{$giID};
					}

					# Genomic
					if ($fileType eq $fileTypeGenomic) {
						# The giID will be the chromosome name (which is currently in $giID)
						# followed by a period and then the position of this site
						# i.e. chr1.position
						my $giIDgenomic = $giID . '.' . '-' . $index;

						# We need to add a line to GIDEF for this site
						# We just set all fields to the same value
						print GIDEF "$giIDgenomic\t$giIDgenomic\t$giIDgenomic\t$giIDgenomic\t$giIDgenomic\n";
						# Save the location of this site for this tag	
						push @{ $giIDLocationsByTag{$vtag} }, $giIDgenomic . '_' . '-' . $index;

						# Push the same thing that we pushed for the location
						push @{ $geneMappingsByTag{$vtag} }, $giIDgenomic . '_' . '-' . $index;
					}
				}
			}

			$last_index = $index;
		}
	}
	LogMsg( "Done creating Virtual Tags file\n");
	return (\%giIDLocationsByTag, \%geneMappingsByTag);
}

####################################
# Validate reference file
####################################
sub ValidateReferenceFiles {
	my ($referenceDir, $referenceType) = @_;

	if (defined $referenceDir) {
		my @refFile = readInRefDir($referenceDir);
		if (scalar @refFile == 0) {
			ErrorMsg("No .fa file found\n");
			return 1;
		}

		my $allFCheck = 1;
		foreach my $f (@refFile) {
			LogMsg("Checking format of $f\n");
			my $fFC = faFormatCheck($f);
			if ($fFC) {

			} else {
				LogMsg("Incorrect fasta format. Reference sequences must follow the standard NCBI GenBank fasta format, with definition lines in the following format: >gi|GI_number|gene_name|......\n$f");
				$allFCheck = 0;
			}
		}

		if ($allFCheck) {
			return 1;
		} else {

		}
	}
	return 0;
}

####################################
#File format check, but only check the first lineToCheck lines
####################################
sub faFormatCheck {
	my $sfa = shift @_;

	my $lineToCheck = $c{lineToCheck};

	my $res = 1;

	my $cmd = "grep \">\" $sfa | wc -l";
	my $defLineNumber = `$cmd`;
	chomp($defLineNumber);
	$defLineNumber =~ s/\D//g;
	if ($defLineNumber == 0) {
		return 0;
	} else {
		my $clCt = 0;
		if (open (SFA, "head -$lineToCheck $sfa |")) {
			my $fline = <SFA>;
			if ($fline =~ /^>/) {
			} else {
				ErrorMsg(" File $sfa does not begin with \> .");
				$res = 0;
				if (close (SFA)) {
				} else {
					ErrorMsg(" Can't close $sfa: $!");
				}
				return $res;
			}
			while (<SFA>) {
				chomp;
				if (/^>/) {
				} elsif (/[^ACGTURYSWKMBDHVNX\-\.]+/i) {
					ErrorMsg(" Illegal base(s) in line\n$_\nfrom $sfa");
					$res = 0;
					last;
				} else {
				}
			}
			if (close (SFA)) {
			} else {
				ErrorMsg(" Can't close $sfa: $!");
			}
		} else {
			ErrorMsg(" Can't open $sfa: $!");
		}
	}
	return $res;
}

####################################
####################################
sub promptUser {
	my($prompt, $default) = @_;
	my $defaultValue = $default ? "[$default]" : "";
	print "\n--> $prompt $defaultValue: ";
	chomp(my $input = <STDIN>);
	return $input ? $input : $default;
}

####################################
####################################
sub LogMsg {
	my($Msg) = @_;

	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();

	my $myTime = sprintf("%02d:%02d:%02d", $hour, $minute, $second);

	print RUNLOG "$myTime: $Msg";

	if ($InteractiveMode) {
		print "$myTime: $Msg";
	}
}

####################################
####################################
sub ErrorMsg {
	my($Msg) = @_;

	LogMsg("**** ERROR: $Msg");
}
