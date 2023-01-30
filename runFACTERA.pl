#!/usr/local/bin/perl


use strict;
use Getopt::Long;
use Cwd qw(cwd);
use POSIX qw(strftime);


=head1 NAME

runFACTERA.pl


=head2 SYNOPSIS

perl runFACTERA.pl -i input-dir -o output-dir


=head3 DESCRIPTION

This PERL script to run multiple jobs using FACTERA via an algoritm of one fusion gene 
by one fusion gene.


=head4 AUTHOR

Xiaokang Pan (Xiaokang.Pan@osumc.edu)


=head5 LAST UPDATE

01/29/2023

=cut


my $usage = <<EOS;
   Usage: perl $0 -i inpath [-options]

   -i|inpath	[string] (input directory storing alignment BAM and their index files )
   -p|pos_file	[file]   (BED file conain position info of targeted intron regions)
   -o|outpath	[string] (output directory)
   -h|help (help information)

EOS

my ($inpath, $pos_file, $outpath, $help);
GetOptions (
  "inpath=s"  => \$inpath,        # input directory contains alignment BAM and their index files 
  "pos_file:s" => \$pos_file,	  # BED file contains position info of targeted intron regions in bed
  "outpath:s" => \$outpath,       # output directory
  "help:s"    => \$help		  # help information
);
die "\n$usage\n" if ($help or
                     !defined($inpath));

my $current_path = cwd;

if ($inpath !~ /^.+\/.+$/) {
    $inpath =  $current_path."/".$inpath;
}

if (!$pos_file) {
     $pos_file = "database/CTP-targeted-intron-regions.bed";
}

if (!$outpath) {

    system("mkdir FACTERA-output");
    $outpath = $current_path."/factera-output";

} elsif ($outpath !~ /^\/.+\/.+$/) {
    $outpath = $current_path."/".$outpath;
}

open (INPUT, "$pos_file") || die $!;
system("mkdir tmp_beds");
my $gene0;
while (<INPUT>) {
    my $line = $_;
    chomp($line);
    next if !$line;

    my @items = split(/\t/, $line);
    my $gene = $items[3];
    if (!$gene0 || $gene0 ne $gene) {
        if ($gene0 && $gene0 ne $gene) {
            close(OUTPUT);
        }
        open (OUTPUT, ">>tmp_beds/$gene.bed") || die $!; 
    } 
    print OUTPUT "$line\n";       
 
    $gene0 = $gene;
}

my @files = <$inpath/*.bam>;
my @bed_files = <tmp_beds/*.bed>;

my %outputpath;

foreach my $file (@files) {

   my $sample_name;
   if ($file =~ /^.+\/(.+)\.bam$/) {
       $sample_name = $1;
   
       my $outpath2 = "$outpath/$sample_name";
       system("mkdir $outpath2");

       foreach my $bed_file (@bed_files) {
           my ($fgene) = $bed_file =~ /^.+\/(.+)\.bed$/;
     
           my $outpath3 = "$outpath2/$fgene";
         
           system("mkdir $outpath3"); 

           $outputpath{$sample_name}{$outpath3} = 1;      

           system("nohup factera.pl -p 4 -s 2 -F -o $outpath3 $file database/exons.bed database/hg19.2bit $bed_file &");
      }
      system("sleep 3m");
   } 
}

system("sleep 15m");

foreach my $sample (keys %outputpath) {

    open (OUTPUT2, ">>$outpath/$sample.factera.fusions.txt") || die $!;
    my $title = "Est_Type\tRegion1\tRegion2 Break1\tBreak2\tBreak_support1\tBreak_support2\tBreak_offset\tOrientation\tOrder1\tOrder2\tBreak_depth\tProper_pair_support\tUnmapped_support\tImproper_pair_support\tPaired_end_depth\tTotal_depth\tFusion_seq\tNon-templated_seq";
    print OUTPUT2 "$title\n";
    my $gene_num = 0;
    foreach my $path (keys %{$outputpath{$sample}}) {
        $gene_num++;

        if (open (INPUT2, "$path/$sample.factera.fusions.txt")) {
            my $line_num = 0;       
            while(<INPUT2>) {
            
                my $line = $_;
                chomp($line);
                next if !$line;

                $line_num++;
            
                if ($line_num == 1) {
                    next;
                } else {
                    print OUTPUT2 "$line\n";
                }
            }
            close(INPUT2);
        }
        system("rm -r $path");
    }
    close(OUTPUT2);
    system("rm -r $outpath/$sample");
}
system("rm -r tmp_beds");
#system("rm nohup.out");
