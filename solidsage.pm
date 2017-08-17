#----------------------------------------------------------
# This package holds the common stuff needed by the various
# solid-sage related programs.
#----------------------------------------------------------
package solidsage;

use strict;
use warnings;

use vars qw(%c);

%c = (
	defaultcsfastaDir => 'final',
	defaultReferenceDir => '/mnt/bulksolid/GCF/VirtualTagReference/RefSeq-mRNA-28',
	giDefFileDefault => 'GIDefinition.tab',
	SequenceCountsFileDefault => 'sequencecounts.tab',
	NMCountsFileDefault => 'nmcounts.tab',
	GeneCountsFileDefault => 'genecounts.tab',
	ReadIDToSequenceOneGeneFileDefault => 'readidtosequenceOneGene.tab',
	AllMappedReadsFileDefault =>'mappedreads-all.txt',
	UniqueMappedReadsFileDefault =>'mappedreads-uniq.txt',
	MultiMappedReadsFileDefault =>'mappedreads-multi.txt',
	IgnoredMappedReadsFileDefault =>'mappedreads-ignored.txt',
	UnmappedReadsFileDefault =>'unmappedreads.csfasta',
	ReadCountsFileDefault =>'readcounts.csv',
	defaultReferenceType => 'mRNA',
	defaultTagLength => 22,
	defaultMismatches => 1,
	defaultSite => 'CATG',
	alternateSite => 'GATC',
	defaultReadSplit => 10000000,
	defaultVirtualTagLength => 28,
	spacer => 'N' x 3,
	lineToCheck => 100,
	errorLogDefault => 'error.log',
	runLogDefault => 'run.log',
	mRNAChoice => 'Map to sense strand only (mRNA), i.e., the provided sequence.',
	aRNAChoice => 'Map to antisense strand only (aRNA), i.e., the reverse complementary of provided sequence.',
	gDNAChoice => 'Map to both strands (gDNA), i.e., the provided and its reverse complementary sequences.',
	mRNA => 'mRNA',
	aRNA => 'aRNA',
	gDNA => 'gDNA',
);
 

1; # need to return a true value from the file

