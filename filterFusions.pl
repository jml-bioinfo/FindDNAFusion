#!/usr/local/bin/perl


use strict;
use Getopt::Long;
use Cwd qw(cwd);


=head1 NAME

filterFusions.pl


=head2 SYNOPSIS

perl filterFusions.pl -i input-directory -o output-file


=head3 DESCRIPTION

This script filters artifacts in the outputs of FACTERA, JuLI and GeneFuse in a run based on 
our testing datasets of over 50 runs. It outputs a file containing potentially reportable 
fusions in an order by sample by sample. 


=head4 AUTHOR

Xiaokang Pan (Xiaokang.Pan@osumc.edu)


=head5 LAST UPDATE

04/06/2021

=cut


my $usage = <<EOS;
   Usage: perl $0 [-options]

   -i|inpath	[string] (path of the directory storing fusion outputs from all the three methods
   -o|outfile	[string] output file
   -h|help (help information)

EOS

my ($inpath, $outfile, $help);
GetOptions (
  "inpath=s"  => \$inpath,        # path of the directory storing input fusion files 
  "outfile:s" => \$outfile,       # output file
  "help:s"    => \$help		  # help information
);
die "\n$usage\n" if ($help or
                     !defined($inpath));
$outfile = "outfile" if !$outfile;

my $current_path = cwd;
if (!$outfile) {

    system("touch fusion-output-file");
    $outfile = $current_path."/fusion-output-file";

} elsif ($outfile !~ /^.+\/.+$/) {
    $outfile = $current_path."/".$outfile;
}

open(OUTPUT, ">$outfile") || die "cannot open $outfile\n";

# read in artifacts
open(INPUT, "input_data/artifacts.txt") || die "cannot open artifact file";

my %artifacts;
while (<INPUT>) {
    my $line = $_;
    chomp($line);
    next if !$line;

    my @items = split(/\t/, $line);
    my $gene1 = $items[0];
    my $gene2 = $items[1];

    $artifacts{$gene1.$gene2} = 1;
    $artifacts{$gene2.$gene1} = 1;
}
close(INPUT);

#read targeted fusion genes
open(TARGET, "input_data/targeted_fusion_genes.txt") || die "cannot open targeted fusion genes file";
my %target_genes;
while (<TARGET>) {
    my $line = $_;
    chomp($line);
    next if !$line;

    $target_genes{$line} = 1;
}
close(TARGET);

#read all fusion output files 
my @f_files = <$inpath/factera/*fusions.txt>;
my @j_files = <$inpath/juli/*annotated.txt>;
my @g_files = <$inpath/genefuse/*fusion.txt>;

#filter one file by one file
my %f_files = map {$_=>1} @f_files if @f_files;
my %j_files = map {$_=>1} @j_files if @j_files;
my %g_files = map {$_=>1} @g_files if @g_files;

my %samples;
my %factera_fusions;
foreach my $f_file (sort {$a cmp $b} keys %f_files) {

    my ($sample_id) = $f_file =~ /^.+\/(.+)\.factera.fusions.txt$/;
    $samples{$sample_id} = 1; 

    open (INPUT1, "$f_file") || die "cannot open $f_file file\n";
    my $count = 0;
    while (<INPUT1>) {

        $count++;
        next if $count == 1;

        my $line = $_;

        chomp($line);

        next if !$line;        

        my @list = split(/\t/, $line);
        my $type = $list[0];
        my $reg1 = $list[1];
        my $reg2 = $list[2];
        my $brk1 = $list[3];
        my $brk2 = $list[4]; 
        my $sup1 = $list[5];
        my $sup2 = $list[6];      
        my $total_dp = $list[16];        

        my $fusion_line;
        $fusion_line = "$type\t$reg1\t$reg2\t$brk1\t$brk2\t$sup1\t$sup2\t$total_dp";

        if (($sup1 < 10 || $sup2 < 10) && !$target_genes{$reg1} && !$target_genes{$reg2}) { 
            next;

        } elsif ($artifacts{$reg1.$reg2} || $artifacts{$reg2.$reg1}) {

            $factera_fusions{$sample_id}{$reg1.$reg2} = "(".$fusion_line.")";

        } else {
            if ($target_genes{$reg1} || $target_genes{$reg2}) {

                $factera_fusions{$sample_id}{$reg1.$reg2} = "*".$fusion_line;
            } else {

                $factera_fusions{$sample_id}{$reg1.$reg2} = $fusion_line;
            }
        }
    }
    close(INPUT1);
}

my %juli_fusions;
foreach my $j_file (sort {$a cmp $b} keys %j_files) {

    my ($sample_id) = $j_file =~ /^.+\/(.+)\.annotated.txt$/;

    open (INPUT2, "$j_file") || die "cannot open $j_file file\n";
    my $count = 0;
    while (<INPUT2>) {

        $count++;
        next if $count == 1;

        my $line = $_;

        chomp($line);

        next if !$line;

        my @list = split(/\t/, $line);
        my $ChrA = $list[0];
        my $BreakA = $list[1];
        my $brkA = $ChrA.":".$BreakA;

        my $OriA = $list[2];
        my $DisA = $list[3];
        my $SplitA = $list[4];
        my $ChrB = $list[5];
        my $BreakB = $list[6];
        my $brkB = $ChrB.":".$BreakB;

        my $OriB = $list[7];
        my $DisB = $list[8];
        my $SplitB = $list[9]; 
        my $Event = $list[10];
        my $GeneA = $list[11];
        my $StrGeneA = $list[12];
        my $GeneB = $list[13];
        my $StrGeneB = $list[14];
        my $InfoA = $list[15];
        my $InfoB = $list[16];
        my $Direction = $list[17];
        my $Frame = $list[18];
        my $Cosmic = $list[19];
        $Cosmic = "No match" if !$Cosmic;

        my $fusion_line;
        if ($target_genes{$GeneA} || $target_genes{$GeneB}) {

            $fusion_line = "*"."$Event\t$GeneA\t$GeneB\t$brkA\t$brkB\t$DisA\t$SplitA\t$Direction\t$DisB\t$SplitB\t$Cosmic";
        } else {

            $fusion_line = "$Event\t$GeneA\t$GeneB\t$brkA\t$brkB\t$DisA\t$SplitA\t$Direction\t$DisB\t$SplitB\t$Cosmic";
        }    

        if (($DisA < 10 || $DisB < 10) && !$target_genes{$GeneA} && !$target_genes{$GeneB}) {
             next;

        } elsif ($artifacts{$GeneA.$GeneB} || $artifacts{$GeneB.$GeneA}) {

            $juli_fusions{$sample_id}{$GeneA.$GeneB} = "(".$fusion_line.")";

        } elsif ($Event eq "Deletion" && $GeneA eq $GeneB) {
             next;

        } elsif ($Event eq "Tandem" && $GeneA eq $GeneB) {
             next;

        } elsif ($Event eq "Inversion" && $GeneA eq $GeneB) {
             next;

        } elsif ($GeneA =~ /Flanking/ || $GeneB =~ /Flanking/) {
             $juli_fusions{$sample_id}{$GeneA.$GeneB} = "<".$fusion_line.">";

        } elsif ($GeneA =~ /$GeneB/ || $GeneB =~ /$GeneA/) {
             next;

        } elsif ($GeneA =~ /Compl_/ && $GeneB =~ /Compl_/) {
             my ($geneA) =  $GeneA =~ /^Compl_(.+)$/;
             my ($geneB) =  $GeneB =~ /^Compl_(.+)$/;

             if ($artifacts{$geneA.$geneB} || $artifacts{$geneB.$geneA}) {
                 next;
             } else {
                 $juli_fusions{$sample_id}{$GeneA.$GeneB} = $fusion_line;
             }
        } else {
             $juli_fusions{$sample_id}{$GeneA.$GeneB} = $fusion_line;
        }
    }
    close(INPUT2);
}

my %genefuse_fusions;
foreach my $g_file (sort {$a cmp $b} keys %g_files) {

    my ($sample_id) = $g_file =~ /^.+\/(.+)\-fusion.txt$/;

    open (INPUT3, "$g_file") || die "cannot open $g_file file\n";
    my $count = 0;
    while (<INPUT3>) {

        $count++;
        next if $count == 1;

        my $line = $_;

        chomp($line);

        next if !$line;

        my ($geneA, $geneB, $total_reads, $unique_reads);
        if ($line =~ /^#Fusion:\s(.+)_ENST.+___(.+)_ENST.+\s\s\(total:\s(\d*),\sunique:(\d*)\)$/) {
            $geneA = $1;
            $geneB = $2;
            $$total_reads = $3;
            $unique_reads = $4;

            if ($$total_reads < 10 && !$target_genes{$geneA} && !$target_genes{$geneB}) {
                next;
            } elsif ($artifacts{$geneA.$geneB} || $artifacts{$geneB.$geneA}) {
                $line = "(".$line.")";
            
            } else {
                if ($target_genes{$geneA} || $target_genes{$geneB}) {
                    $line = "*".$line;
                } else {
                    $genefuse_fusions{$sample_id}{$geneA.$geneB} =  "$line\n";
                }
            }
        }
    }
    close(INPUT3);
}

my $sample_num = 0;
print OUTPUT "Note: A row included in () shows an artifact fusion. A row contianed in <> indicates that the fusion gene or partner is flanking of an element. A row with * marks containg a targeted fusion gene.\n\n";
foreach my $sample (sort {$a cmp $b} keys %samples) {
    $sample_num++;

    if ($sample_num == 1) {
        print OUTPUT "------ Sample $sample ------\n\n";
    } else {
        print OUTPUT "\n\n------ Sample $sample ------\n\n";
    } 

    print OUTPUT "Method FACTERA:\n\n";

    my $f_fusions = 0;
    foreach my $fusion (sort {$a cmp $b} keys %{$factera_fusions{$sample}}) {
        $f_fusions++;
        if ($f_fusions == 1) {
            print OUTPUT "Est_Type\tGene A\tGene B\tBreak A\tBreak B\tBreak_support A\tBreak_support B\tTotal_Depth\n";
        }
        my $fusion_line = $factera_fusions{$sample}{$fusion};    
        print OUTPUT "$fusion_line\n";
    }
    print OUTPUT "No fusion detected\n" if $f_fusions == 0;

    print OUTPUT "\n\nMethod JULI:\n\n";
    
    my $j_fusions = 0;
    foreach my $fusion (sort {$a cmp $b} keys %{$juli_fusions{$sample}}) {
        $j_fusions++;
        if ($j_fusions == 1) {
            print OUTPUT "Event\tGene A\tGene B\tBreak A\tBreak B\tDiscord A\tSplit A\tDirection\tDiscord B\tSplit B\tCosmic\n";
        }
        my $fusion_line2 = $juli_fusions{$sample}{$fusion};
        print OUTPUT "$fusion_line2\n";
    }
    print OUTPUT "No fusion detected\n" if $j_fusions == 0;

    print OUTPUT "\n\nMethod GeneFuse:\n\n";
    
    my $g_fusions = 0;
    foreach my $fusion (sort {$a cmp $b} keys %{$genefuse_fusions{$sample}}) {
        $g_fusions++;
        my $fusion_line3 = $genefuse_fusions{$sample}{$fusion};
        print OUTPUT "$fusion_line3\n";
    }
    print OUTPUT "No fusion detected\n" if $g_fusions == 0;

}


# to trim the white spaces at ends of a line or string
sub trimEndWhiteSpace {
    my $str = shift;
    $str =~ s/^\s+//;  # strip leading whitespace
    $str =~ s/\s+$//;  # strip trailing whitespace
    return $str;
}

