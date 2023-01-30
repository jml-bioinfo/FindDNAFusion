#!/usr/local/bin/perl


use strict;
use Getopt::Long;
use Cwd qw(cwd);
use POSIX qw(strftime);


=head1 NAME

runGeneFuse.pl


=head2 SYNOPSIS

perl runGeneFuse.pl -i input-dir -o output-dir


=head3 DESCRIPTION

This PERL script to run multiple jobs using GeneFuse.


=head4 AUTHOR

Xiaokang Pan (Xiaokang.Pan@osumc.edu)


=head5 LAST UPDATE

01/29/2023

=cut


my $usage = <<EOS;
   Usage: perl $0 -i inpath -r ref_genome [-options]

   -i|inpath	  [string] (input sequence directory storing fastq files)
   -r|ref_genome  [string] (reference genome file in FASTA format)
   -o|outpath	  [string] (output directory)
   -c|cpus        [string] CPU number per sample to be used
   -h|help        (help information)

EOS

my ($inpath, $outpath, $ref_genome, $cpus, $help);
GetOptions (
  "inpath=s"     => \$inpath,        # input dir 
  "ref_genome=s" => \$ref_genome,    # reference genome
  "outpath:s"    => \$outpath,       # output dir
  "cpus:s"       => \$cpus,          # number of CPUs
  "help:s"       => \$help	     # help information
);
die "\n$usage\n" if ($help or
                     !defined($inpath) or
                     !defined($ref_genome));

my $current_path = cwd;

if ($inpath !~ /^\/.+\/.+$/) {
    $inpath =  $current_path."/".$inpath;
}

if (!$outpath) {
    system("mkdir genefuse-output");
    $outpath = $current_path."/genefuse-output";

} elsif ($outpath !~ /^\/.+\/.+$/) {
    $outpath = $current_path."/".$outpath;
}

if (!$cpus) {
    $cpus = 12;
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

#my $ref_file = "database/hg19.fasta";
my $fusionGenes = "database/PossibleFusionGenes.csv";

foreach my $sample (sort keys %samples) {

   next if ($sample =~ /-RNA$/ || $sample =~ /-R$/);

   my $r1_fastq_file = "$inpath/$sample"."_R1.fastq.gz";
   my $r2_fastq_file = "$inpath/$sample"."_R2.fastq.gz";
   system("nohup genefuse -r $ref_genome -f $fusionGenes -1 $r1_fastq_file -2 $r2_fastq_file -t $cpus -u 2 -h $outpath/$sample.html > $outpath/$sample-fusion.txt &");

}


