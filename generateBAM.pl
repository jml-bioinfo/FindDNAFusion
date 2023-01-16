#!/usr/local/bin/perl


use strict;
use Getopt::Long;
use Cwd qw(cwd);
use POSIX qw(strftime);


=head1 NAME

runGeneFuse.pl


=head2 SYNOPSIS

perl generateBAM.pl -i input-dir -o output-dir


=head3 DESCRIPTION

This PERL script to generate BAM files for multiple samples.


=head4 AUTHOR

Xiaokang Pan (Xiaokang.Pan@osumc.edu)


=head5 LAST UPDATE

01/02/2023

=cut


my $usage = <<EOS;
   Usage: perl $0 [-options]

   -i|inpath	[string] (input dir storing fastq files at directory)
   -o|outpath	[string] (output dir storng outputs from GeneFuse)
   -c|cpus      [string] CPU number
   -t|time      [string] (expected waiting time in minutes to start actual run)
   -h|help (help information)

EOS

my ($inpath, $outpath, $cpus, $time, $help);
GetOptions (
  "inpath=s"  => \$inpath,        # input dir 
  "outpath:s" => \$outpath,       # output dir
  "cpus:s"    => \$cpus,          # number of CPUs
  "time:s"    => \$time,          # expected waiting time in minutes
  "help:s"    => \$help		  # help information
);
die "\n$usage\n" if ($help or
                     !defined($inpath));

my $current_path = cwd;

if ($inpath !~ /^\/.+\/.+$/) {
    $inpath =  $current_path."/".$inpath;
}

if (!$outpath) {
    system("mkdir mbam");
    $outpath = $current_path."/mbam";

} elsif ($outpath !~ /^\/.+\/.+$/) {
    $outpath = $current_path."/".$outpath;
}

if (!$cpus) {
    $cpus = 16;
}

if ($time) {
   $time  = $time."m";
   system("sleep $time");
}

my @files = <$inpath/*.fastq.gz>;

my %samples;
foreach my $file (@files) {

    my $sample_name;
        
    if ($file =~ /^.+\/(.+)_R\d\.fastq\.gz$/) {
        $sample_name = $1;
        $samples{$sample_name} = 1;
    }
}

my $ref_file = "database/hg19.fasta";

# run BWA MEM to map on Human Genome
foreach my $sample (sort keys %samples) {

   next if ($sample =~ /-RNA$/ || $sample =~ /-R$/);

   my $r1_fastq_file = "$inpath/$sample"."_R1.fastq.gz";
   my $r2_fastq_file = "$inpath/$sample"."_R2.fastq.gz";
   system("bwa mem -M -t $cpus $ref_file $r1_fastq_file $r2_fastq_file | samtools sort -@ $cpus -o $outpath/$sample.bam - &");
}
system("sleep 90m");

# check if mapping completed.
# if not, wait for a few minutes
my @tmp_files = <$outpath/*.tmp.*>;
while (@tmp_files) {
    system("sleep 5m");
    @tmp_files = <$outpath/*.tmp.*>;
}

# make index
my @bam_files = <$outpath/*.bam>;
foreach my $bam_file (@bam_files) {

    system("samtools index $bam_file &");
}
system("sleep 10m");
print "BAM generation is completed.\n\n";

