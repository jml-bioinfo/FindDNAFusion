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

02/21/2021

=cut


my $usage = <<EOS;
   Usage: perl $0 [-options]

   -i|inpath	    [string] (input dir storing BAM files)
   -t|time	    [string] (expected waiting time in minutes to start actual run)
   -o|outpath	    [string] (output dir storng outputs from JuLI)
   -h|help          (help information)

EOS

my ($inpath, $time, $outpath, $help);
GetOptions (
  "inpath=s"  => \$inpath,        # input dir 
  "time:s"    => \$time,          # expected waiting time in minutes
  "outpath:s" => \$outpath,       # output dir
  "help:s"    => \$help		  # help information
);
die "\n$usage\n" if ($help or
                     !defined($inpath));

my $current_path = cwd;
my $this_path = "/data/lab-genomics/LNGS/FUSION";
if (!$outpath) {

    system("mkdir JuLI-output");
    $outpath = $current_path."/JuLI-output";

} elsif ($outpath !~ /^.+\/.+$/) {
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

    system("./run-JuLI-single-sample.R $file $sample_name $outpath $current_path &");
}

