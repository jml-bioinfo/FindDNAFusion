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

   -i|inpath	    [string] (input dir storing BAM files)
   -t|time	    [string] (expected waiting time in minutes to start actual run)
   -r|ref_genome    [string] (reference genome file with full path)
   -o|outpath	    [string] (output dir storng outputs from JuLI)
   -h|help          (help information)

EOS

my ($inpath,  $time, $ref_genome, $outpath, $help);
GetOptions (
  "inpath=s"     => \$inpath,        # input dir 
  "time:s"       => \$time,          # expected waiting time in minutes
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

if ($time) {
   $time  = $time."m";
   system("sleep $time");
}

my @files = <$inpath/*.bam>;

foreach my $file (@files) {

    my $sample_name;
        
    if ($file =~ /^.+\/(.+)\.bam$/) {
        $sample_name = $1;
    }

    system("./run-JuLI-single-sample.R $file $sample_name $ref_genome $outpath &");
}

