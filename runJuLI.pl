#!/usr/bin/perl


use strict;
use Getopt::Long;
use Cwd qw(cwd);
use POSIX qw(strftime);


=head1 NAME

runJuLI.pl


=head2 SYNOPSIS

perl runJuLI.pl -i input-dir -o output-dir


=head3 DESCRIPTION

This PERL script to run multiple jobs simultaneously using JuLI to detect gene fusions.


=head4 AUTHOR

Xiaokang Pan (Xiaokang.Pan@osumc.edu)


=head5 LAST UPDATE

01/25/2023

=cut


my $usage = <<EOS;
   Usage: perl $0 -i inpath [-options]

   -i|inpath	    [string] (input directory storing BAM and their index files)
   -r|ref_genome    [file] (reference genome file in FASTA format)
   -o|outpath	    [string] (output directory)
   -h|help          (help information)

EOS

my ($inpath, $ref_genome, $outpath, $help);
GetOptions (
  "inpath=s"     => \$inpath,        # input dir 
  "ref_genome:s" => \$ref_genome,    # reference genome
  "outpath:s"    => \$outpath,       # output dir
  "help:s"       => \$help	     # help information
);
die "\n$usage\n" if ($help or
                     !defined($inpath));

my $current_path = cwd;

if ($inpath !~ /^.+\/.+$/) {
    $inpath =  $current_path."/".$inpath;
}

if (!$ref_genome) {
    $ref_genome = "$current_path/database/hg19.fasta";
}

if (!$outpath) {

    system("mkdir JuLI-output");
    $outpath = $current_path."/JuLI-output";

} elsif ($outpath !~ /^\/.+\/.+$/) {
    $outpath = $current_path."/".$outpath;
}

my @files = <$inpath/*.bam>;

foreach my $file (@files) {

    my $sample_name;
        
    if ($file =~ /^.+\/(.+)\.bam$/) {
        $sample_name = $1;
    }

    system("./run-JuLI-single-sample.R $file $sample_name $ref_genome $outpath &");
}

