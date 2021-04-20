#!/usr/local/bin/perl


use strict;
use Getopt::Long;
use Cwd qw(cwd);
use POSIX qw(strftime);


=head1 NAME

preprocessData.pl


=head2 SYNOPSIS

perl preprocessData.pl -i input-dir -o output-dir


=head3 DESCRIPTION

This PERL script integrates the WDL-coded workflows paired-fastq-to-unmapped-bam.wdl and  
processing-for-fusion-detection.wdl from GATK Best Practices to generate mapped    
BAM files from a batch of fastq raw sequence files for variant dection at a local 
computer. It takes a batch of fastq files and output a set of BAM files.


=head4 AUTHOR

Xiaokang Pan (Xiaokang.Pan@osumc.edu)


=head5 LAST UPDATE

04/19/2021

=cut


my $usage = <<EOS;
   Usage: perl $0 [-options]

   -i|inpath	[string] (full path of input dir storing fastq files in '.gz' format)
   -o|outpath	[string] (path of output dir storing bam files)
   -h|help (help information)

EOS

my ($inpath, $outpath, $help);
GetOptions (
  "inpath=s"  => \$inpath,        # full path of input dir storing fastq files
  "outpath:s" => \$outpath,       # path of  output fir storing BAM files
  "help:s"    => \$help		  # help information
);
die "\n$usage\n" if ($help or
                     !defined($inpath));

my $current_path = cwd;
if (!$outpath) {

    system("mkdir mbam");
    $outpath = $current_path."/mbam";

} elsif ($outpath !~ /^\/.+\/.+$/) {
    $outpath = $current_path."/".$outpath;
}

my @in_files = <$inpath/*.fastq.gz>;
my $nsamples = @in_files /2;

my %in_files = map {$_=>1} @in_files if @in_files;

#generate paired-fastq-to-unmapped-bam_inputs.json
open (OUTPUT1, ">paired-fastq-to-ubam_inputs.json") || die "cannot open paired-fastq-to-ubam_inputs.json";

my $date = strftime("%y-%m-%d", localtime(time));
print OUTPUT1 "{\n";
print OUTPUT1 "   \"ConvertPairedFastQsToUnmappedBamWf.readgroup_list\": [\n";

my $file_count = 0;
my %dup;
my $raw_seq_path;
foreach my $file (sort keys %in_files) {
    $file_count++;

    my $file_prefix;
    if ($file =~ /^(.+)\/(.+).R\d.fastq.gz$/) {
        $raw_seq_path = $1;
        $file_prefix  = $2;
    }

    next if $dup{$file_prefix};

    if ($file_count == 1) {
        print OUTPUT1 "    \"$file_prefix\"";
    } else {
        print OUTPUT1 ", \"$file_prefix\"";
    }
    $dup{$file_prefix} = 1;
} 
print OUTPUT1 "\n   ],\n";

print "\n\nprocessing from raw reads to ubam ..........\n\n";

system("mkdir ubam");
print OUTPUT1 "   \"ConvertPairedFastQsToUnmappedBamWf.output_loc\": ";
print OUTPUT1 "   \"$current_path/ubam\",\n";
print OUTPUT1 "   \"ConvertPairedFastQsToUnmappedBamWf.PairedFastQsToUnmappedBAM.disk_size\": 200,\n";
print OUTPUT1 "   \"ConvertPairedFastQsToUnmappedBamWf.PairedFastQsToUnmappedBAM.mem_size\": \"2 GB\",\n";
print OUTPUT1 "   \"ConvertPairedFastQsToUnmappedBamWf.metadata\":  {\n";

my $file_count = 0;
my %dup2;
foreach my $file (sort keys %in_files) {
    $file_count++;
     
    my $file_prefix;
    if ($file =~ /^(.+)\/(.+).R\d.fastq.gz$/) {
        $raw_seq_path = $1;
        $file_prefix  = $2;
    }

    next if $dup2{$file_prefix};
    
    if ($file_count > 1) {
        print OUTPUT1 "      ],\n";
    }
    print OUTPUT1 "    \"$file_prefix\": [\n";
    print OUTPUT1 "       \"$file_prefix\", \"TruSeqLT\", \"NextSeq\", \"$date\", \"Illumina\", \"JML\"\n";

    $dup2{$file_prefix} = 1;
}
print OUTPUT1 "    ]\n";
print OUTPUT1 "  },\n";

$file_count = 0;
my %dup3;
print OUTPUT1 "  \"ConvertPairedFastQsToUnmappedBamWf.fastq_pairs\": {\n";
foreach my $file (sort keys %in_files) {
    $file_count++;

    my $file_prefix;
    if ($file =~ /^(.+)\/(.+).R\d.fastq.gz$/) {
        $raw_seq_path = $1;
        $file_prefix  = $2;
    }

    next if $dup3{$file_prefix};

    if ($file_count > 1) {
        print OUTPUT1 "     ],\n";
    }

    print OUTPUT1 "     \"$file_prefix\": [\n";

    my $c = 0;
 
    my ($file_prefix2) = $file =~ /^.+\/(.+_R)\d.fastq.gz$/; 
    foreach my $file (keys %in_files) {
         if ($file =~ /$file_prefix2/) {
             $c++;
             if ($c == 1) {
                  print OUTPUT1 "        \"$file\",\n";
                  
             } else {
                  print OUTPUT1 "        \"$file\"\n";
             }
         }
        
    }
    $dup3{$file_prefix} = 1;
}
print OUTPUT1 "     ]\n";
print OUTPUT1 "  }\n";
print OUTPUT1 "}\n";
close(OUTPUT1);

# run paired-fastq-to-unmapped-bam.wdl
system("java -jar \$CROMWELL run paired-fastq-to-unmapped-bam.wdl --inputs paired-fastq-to-ubam_inputs.json");

# start to map sequences to human genome ..... to generate bam files
# get ubam files
print "\n\nStart sequence mapping from ubam to mbam ...... \n\n";
my @ubam_files = <$current_path/ubam/*.unmapped.bam>;
my %ubam_files = map {$_=>1} @ubam_files if @ubam_files;
my $file_count = 0;
foreach my $file (sort {$a cmp $b} keys %ubam_files) {
    $file_count++;
    my ($file_prefix) = $file =~ /^.+\/(.+)\.unmapped.bam$/;

    # generate processing-for-variant-discovery-gatk4.hg19.LNGS.inputs.json for each sample dataset
    my $json_file = "processing-for-variant-discovery-gatk4.hg19.LNGS.inputs".$file_count.".json";
    open(OUTPUT2, ">$json_file");
    print OUTPUT2 "{\n";
    print OUTPUT2 "   \"##_COMMENT1\": \"SAMPLE NAME AND UNMAPPED BAMS\",\n";
    print OUTPUT2 "   \"PreProcessingForVariantDiscovery_GATK4.sample_name\": \"$file_prefix\",\n";
    print OUTPUT2 "   \"PreProcessingForVariantDiscovery_GATK4.ref_name\": \"hg19\",\n";
    print OUTPUT2 "   \"PreProcessingForVariantDiscovery_GATK4.flowcell_unmapped_bams_list\": [\"$file\"],\n";
    print OUTPUT2 "   \"PreProcessingForVariantDiscovery_GATK4.unmapped_bam_suffix\": \".bam\",\n";
    print OUTPUT2 "   \"PreProcessingForVariantDiscovery_GATK4.output_loc\": \"$outpath\",\n";   
 
my $rest_text = <<EOF;

    "##_COMMENT2": "REFERENCE FILES",
    "PreProcessingForVariantDiscovery_GATK4.ref_dict": "/data/reference-genomes/Homo_sapiens/UCSC/hg19/Sequence/GATKBundle/hg19.dict",
    "PreProcessingForVariantDiscovery_GATK4.ref_fasta": "/data/reference-genomes/Homo_sapiens/UCSC/hg19/Sequence/GATKBundle/hg19.fasta",
    "PreProcessingForVariantDiscovery_GATK4.ref_fasta_index": "/data/reference-genomes/Homo_sapiens/UCSC/hg19/Sequence/GATKBundle/hg19.fasta.fai",
    "PreProcessingForVariantDiscovery_GATK4.SamToFastqAndBwaMem.ref_sa": "/data/reference-genomes/Homo_sapiens/UCSC/hg19/Sequence/GATKBundle/hg19.fasta.sa",
    "PreProcessingForVariantDiscovery_GATK4.SamToFastqAndBwaMem.ref_amb": "/data/reference-genomes/Homo_sapiens/UCSC/hg19/Sequence/GATKBundle/hg19.fasta.amb",
    "PreProcessingForVariantDiscovery_GATK4.SamToFastqAndBwaMem.ref_bwt": "/data/reference-genomes/Homo_sapiens/UCSC/hg19/Sequence/GATKBundle/hg19.fasta.bwt",
    "PreProcessingForVariantDiscovery_GATK4.SamToFastqAndBwaMem.ref_ann": "/data/reference-genomes/Homo_sapiens/UCSC/hg19/Sequence/GATKBundle/hg19.fasta.ann",
    "PreProcessingForVariantDiscovery_GATK4.SamToFastqAndBwaMem.ref_pac": "/data/reference-genomes/Homo_sapiens/UCSC/hg19/Sequence/GATKBundle/hg19.fasta.pac",


    "##_COMMENT3": "KNOWN SITES RESOURCES",
    "PreProcessingForVariantDiscovery_GATK4.dbSNP_vcf": "/data/reference-genomes/Homo_sapiens/UCSC/hg19/Annotation/Variation/dbsnp_138.hg19.vcf",
    "PreProcessingForVariantDiscovery_GATK4.dbSNP_vcf_index": "/data/reference-genomes/Homo_sapiens/UCSC/hg19/Annotation/Variation/dbsnp_138.hg19.vcf.idx",
    "PreProcessingForVariantDiscovery_GATK4.known_indels_sites_VCFs": [
    "/data/reference-genomes/Homo_sapiens/UCSC/hg19/Annotation/Variation/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
    ],
    "PreProcessingForVariantDiscovery_GATK4.known_indels_sites_indices": [
      "/data/reference-genomes/Homo_sapiens/UCSC/hg19/Annotation/Variation/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx"
    ],

    "##_COMMENT4": "MISC PROGRAM PARAMETERS",
    "PreProcessingForVariantDiscovery_GATK4.bwa_commandline": "bwa mem -K 100000000 -p -v 3 -t 4 -Y \$bash_ref_fasta",
    "PreProcessingForVariantDiscovery_GATK4.compression_level": 5,
    "PreProcessingForVariantDiscovery_GATK4.SamToFastqAndBwaMem.num_cpu": "4",

    "##_COMMENT5": "DOCKERS",
    "## PreProcessingForVariantDiscovery_GATK4.gotc_docker": "String? (optional)",
    "## PreProcessingForVariantDiscovery_GATK4.gatk_docker": "String? (optional)",
    "## PreProcessingForVariantDiscovery_GATK4.picard_docker": "String? (optional)",
    "## PreProcessingForVariantDiscovery_GATK4.python_docker": "python:2.7",

    "##_COMMENT6": "PATHS",
    "PreProcessingForVariantDiscovery_GATK4.gotc_path": "/usr/bin/",
    "PreProcessingForVariantDiscovery_GATK4.picard_path": "/usr/local/lib/",
    "PreProcessingForVariantDiscovery_GATK4.gatk_path": "/usr/local/lib/gatk-4.1.4.1/",

    "##_COMMENT7": "JAVA OPTIONS",
    "PreProcessingForVariantDiscovery_GATK4.SamToFastqAndBwaMem.java_opt": "-Xms3000m",
    "PreProcessingForVariantDiscovery_GATK4.MergeBamAlignment.java_opt": "-Xms3000m",
    "PreProcessingForVariantDiscovery_GATK4.MarkDuplicates.java_opt": "-Xms4000m",
    "PreProcessingForVariantDiscovery_GATK4.SortAndFixTags.java_opt_sort": "-Xms4000m",
    "PreProcessingForVariantDiscovery_GATK4.SortAndFixTags.java_opt_fix": "-Xms500m",
    "PreProcessingForVariantDiscovery_GATK4.BaseRecalibrator.java_opt": "-Xms4000m",
    "PreProcessingForVariantDiscovery_GATK4.GatherBqsrReports.java_opt": "-Xms3000m",
    "PreProcessingForVariantDiscovery_GATK4.ApplyBQSR.java_opt": "-Xms3000m",
    "PreProcessingForVariantDiscovery_GATK4.GatherBamFiles.java_opt": "-Xms2000m",

    "##_COMMENT8": "MEMORY ALLOCATION",
    "PreProcessingForVariantDiscovery_GATK4.GetBwaVersion.mem_size": "1 GB",
    "PreProcessingForVariantDiscovery_GATK4.SamToFastqAndBwaMem.mem_size": "14 GB",
    "PreProcessingForVariantDiscovery_GATK4.MergeBamAlignment.mem_size": "3500 MB",
    "PreProcessingForVariantDiscovery_GATK4.MarkDuplicates.mem_size": "7 GB",
    "PreProcessingForVariantDiscovery_GATK4.SortAndFixTags.mem_size": "5000 MB",
    "PreProcessingForVariantDiscovery_GATK4.CreateSequenceGroupingTSV.mem_size": "2 GB",
    "PreProcessingForVariantDiscovery_GATK4.BaseRecalibrator.mem_size": "6 GB",
    "PreProcessingForVariantDiscovery_GATK4.GatherBqsrReports.mem_size": "3500 MB",
    "PreProcessingForVariantDiscovery_GATK4.ApplyBQSR.mem_size": "3500 MB",
    "PreProcessingForVariantDiscovery_GATK4.GatherBamFiles.mem_size": "3 GB",

    "##_COMMENT9": "DISK SIZE ALLOCATION",
    "PreProcessingForVariantDiscovery_GATK4.agg_small_disk": 200,
    "PreProcessingForVariantDiscovery_GATK4.agg_medium_disk": 300,
    "PreProcessingForVariantDiscovery_GATK4.agg_large_disk": 400,
    "PreProcessingForVariantDiscovery_GATK4.flowcell_small_disk": 100,
    "PreProcessingForVariantDiscovery_GATK4.flowcell_medium_disk": 200,

    "##_COMMENT10": "PREEMPTIBLES",
    "PreProcessingForVariantDiscovery_GATK4.preemptible_tries": 4,
    "PreProcessingForVariantDiscovery_GATK4.agg_preemptible_tries": 4
EOF
 
    print OUTPUT2 "$rest_text\n";
    print OUTPUT2 "}\n";
    close(OUTPUT2);

    # run processing-for-variant-discovery-gatk4.wdl
    system("nohup java -jar \$CROMWELL run processing-for-variant-discovery-gatk4.wdl --inputs $json_file > log$file_count &");

}

system("sleep 120m");

my @bam_files = <$outpath/*.bam>;
my $bam_file_num = scalar @bam_files;

while ($bam_file_num < $nsamples) {
   system("sleep 10m");
   my @bam_files2 = <$outpath/*.bam>;
   $bam_file_num = @bam_files2;
}
system("sleep 20m");

system("rename .hg19.bam .bam $outpath/*.hg19.bam");
system("rename .hg19.bai .bam.bai $outpath/*.hg19.bai");
#system("rename .bai .bam.bai $outpath/*.bai");

# remove all the json files  
 foreach (my $i=1; $i<= $bam_file_num; $i++) {
     my $json_file = "processing-for-variant-discovery-gatk4.hg19.LNGS.inputs".$i.".json";
     system("rm $json_file");
 }

system("rm -r ubam");
system("rm -r cromwell-executions");
system("rm -r cromwell-workflow-logs");
system("rm *log*");
print "\n\n Done sequence mapping ...... \n\n";
