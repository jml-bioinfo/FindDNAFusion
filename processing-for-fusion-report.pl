#!/usr/local/bin/perl


use strict;
use Getopt::Long;
use Cwd qw(cwd);


=head1 NAME

processing-for-fusion-report.pl


=head2 SYNOPSIS

perl processing-for-fusion-report.pl -i input-directory -o output-directory


=head3 DESCRIPTION

    This script parses outputs of JuLI, FACTERA and GeneFuse fusion callers, filters out common 
tool-specific artifactual calls based on a blacklist of artifacts we established, selects 
reportable fusions according to critera we developed and annotates the selected fusions by importing
a visualization application. It generates files "Possible-fusions.txt", "Reportable-fusions.txt" 
and a set of annotated fusion files per file per sample in PDF format.

Selection criteria for reportable fusions in DNA:

1) 2 or 3 methods agreed with 20 or more supporting reads at two breakpoints or agreed 
with 10 or more at each breakpoints;
2) 1 agreed with strong supporting reads (>=50) in total and with reciprocal calls and 
found in the “known-fusion” database (or expected in the tumor being tested to be 
implemented)


=head4 AUTHOR

Xiaokang Pan (Xiaokang.Pan@osumc.edu)


=head5 LAST UPDATE

12/28/2022

=cut


my $usage = <<EOS;
   Usage: perl $0 [-options]

   -i|inpath	[string] (path of the directory storing fusion outputs from all the three methods)
   -b|bampath   [string] (path of the directory storing BAM files
   -o|outpath	[string] (path of the directory storing  results)
   -h|help (help information)

EOS

my ($inpath, $outpath, $bampath, $help);
GetOptions (
  "inpath=s"  => \$inpath,        # path of the directory storing input fusion files 
  "bampath=s" => \$bampath,       # path of the directory storing BAM files
  "outpath:s" => \$outpath,       # output path
  "help:s"    => \$help		  # help information
);
die "\n$usage\n" if ($help or
                     !defined($inpath) or
                     !defined($bampath));

my $current_path = cwd;
if (!$outpath) {
    $outpath = $current_path;
}

my $outfile  = "Possible-fusions.txt";
my $outfile2 = "Reportable-fusions.txt";

open(OUTPUT, ">$outpath/$outfile") || die "cannot open $outfile\n";
open(OUTPUT2, ">$outpath/$outfile2") || die "cannot open $outfile2\n";

# read in artifacts
open(INPUT, "database/artifacts.txt") || die "cannot open artifact file";

my %artifacts;
while (<INPUT>) {
    my $line = $_;
    chomp($line);
    next if !$line;

    my @items = split(/\t/, $line);
    my $gene1 = $items[0];
    my $gene2 = $items[1];

    $artifacts{"$gene1-$gene2"} = 1;
    $artifacts{"$gene2-$gene1"} = 1;
}
close(INPUT);

#read targeted fusion genes
open(TARGET, "database/targeted_fusion_genes.txt") || die "cannot open targeted fusion genes file";
my %target_genes;
while (<TARGET>) {
    my $line = $_;
    chomp($line);
    next if !$line;

    $target_genes{$line} = 1;
}
close(TARGET);

#read known fusions
open(FUSION, "database/known-fusions.txt") || die "cannot open known fusion file";
my %known_fusions;
while (<FUSION>) {
    my $line = $_;
    chomp($line);
    next if !$line;

    my @list = split(/\t/, $line);
    my ($geneA, $geneB) = ($list[0], $list[1]);

    $known_fusions{"$geneA-$geneB"} = 1;
}
close(FUSION);

#read all fusion output files 
my @f_files = <$inpath/factera/*fusions.txt>;
my @j_files = <$inpath/juli/*annotated.txt>;
my @g_files = <$inpath/genefuse/*fusion.txt>;

#filter one file by one file
my %f_files = map {$_=>1} @f_files if @f_files;
my %j_files = map {$_=>1} @j_files if @j_files;
my %g_files = map {$_=>1} @g_files if @g_files;

# parse FACTERA output
my %samples;
my (%factera_fusions, %factera_fusions2, %factera_fusions3);
my (%factera_fusions_org, %factera_fusions2_org, %factera_fusions3_org);
my $factera_fusions_title;
foreach my $f_file (sort {$a cmp $b} keys %f_files) {

    my ($sample_id) = $f_file =~ /^.+\/(.+)\.factera.fusions.txt$/;
    $samples{$sample_id} = 1; 

    open (INPUT1, "$f_file") || die "cannot open $f_file file\n";
    my $count = 0;
    while (<INPUT1>) {

        $count++;

        my $line = $_;
        if ($count == 1) {
            $factera_fusions_title = $line;
        }
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
        my $discordant_mates = $list[14];      
        my $total_dp = $list[16];        

        my $fusion_line;
        $fusion_line = "$type\t$reg1\t$reg2\t$brk1\t$brk2\t$sup1\t$sup2\t$discordant_mates\t$total_dp";

        if (($sup1 < 5 || $sup2 < 5) && !$target_genes{$reg1} && !$target_genes{$reg2}) { 
            next;

        } elsif ($artifacts{"$reg1-$reg2"} || $artifacts{"$reg2-$reg1"}) {

            #$factera_fusions{$sample_id}{"$reg1-$reg2"}{"$brk1-$brk2"} = "(".$fusion_line.")";

        } else {
            if ($target_genes{$reg1} || $target_genes{$reg2}) {

                $factera_fusions{$sample_id}{"$reg1-$reg2"}{"$brk1-$brk2"} = "*".$fusion_line;
                $factera_fusions_org{$sample_id}{"$reg1-$reg2"}{"$brk1-$brk2"} = $line;
                
                if ($sup1 >= 20 || $sup2 >= 20 || ($sup1 >= 10 && $sup2 > 10)) {
                    $factera_fusions2{$sample_id}{"$reg1-$reg2"}{"$brk1-$brk2"} = "*".$fusion_line;
                    $factera_fusions2_org{$sample_id}{"$reg1-$reg2"}{"$brk1-$brk2"} = $line;
                }

                if ($sup1 >= 50 || $sup2 >= 50) {
                    $factera_fusions3{$sample_id}{"$reg1-$reg2"}{"$brk1-$brk2"} = "*".$fusion_line;
                    $factera_fusions3_org{$sample_id}{"$reg1-$reg2"}{"$brk1-$brk2"} = $line;
                }
            } else {
                $factera_fusions{$sample_id}{"$reg1-$reg2"}{"$brk1-$brk2"} = $fusion_line;
                $factera_fusions_org{$sample_id}{"$reg1-$reg2"}{"$brk1-$brk2"} = $line;
                
                if ($sup1 >= 20 || $sup2 >= 20 || ($sup1 >= 10 && $sup2 > 10)) {
                    $factera_fusions2{$sample_id}{"$reg1-$reg2"}{"$brk1-$brk2"} = $fusion_line;
                    $factera_fusions2_org{$sample_id}{"$reg1-$reg2"}{"$brk1-$brk2"} = $line;
                } 

                if ($sup1 >= 50 || $sup2 >= 50) {
                    $factera_fusions3{$sample_id}{"$reg1-$reg2"}{"$brk1-$brk2"} = "*".$fusion_line;
                    $factera_fusions3_org{$sample_id}{"$reg1-$reg2"}{"$brk1-$brk2"} = $line;
                }
            }
        }
    }
    close(INPUT1);
}

#parse JuLI outputs
my (%juli_fusions, %juli_fusions2, %juli_fusions3);
my (%juli_fusions_org, %juli_fusions2_org, %juli_fusions3_org);
my $juli_fusions_title;
foreach my $j_file (sort {$a cmp $b} keys %j_files) {

    my ($sample_id) = $j_file =~ /^.+\/(.+)\.annotated.txt$/;

    open (INPUT2, "$j_file") || die "cannot open $j_file file\n";
    my $count = 0;
    while (<INPUT2>) {

        $count++;
        my $line = $_;
        if ($count == 1) {
            $juli_fusions_title = $line;
        }

        chomp($line);

        next if !$line;

        my @list = split(/\t/, $line);
        my $ChrA = $list[0];
        my $BreakA = $list[1];
        my $brkA = "$ChrA:$BreakA";

        my $OriA = $list[2];
        my $DisA = $list[3];
        my $SplitA = $list[4];
        my $ChrB = $list[5];
        my $BreakB = $list[6];
        my $brkB = "$ChrB:$BreakB";

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

        if (($SplitA < 10 || $SplitB < 10) && !$target_genes{$GeneA} && !$target_genes{$GeneB}) {
             next;
        } 

        if ($artifacts{"$GeneA-$GeneB"} || $artifacts{"$GeneB-$GeneA"}) {
             next;
        } 

        if ($GeneA =~ /$GeneB/ || $GeneB =~ /$GeneA/) {
             next;
        }
       
        if ($GeneA =~ /Flanking/ || $GeneB =~ /Flanking/) {
            $juli_fusions{$sample_id}{"$GeneA-$GeneB"}{"$brkA-$brkB"} = "<".$fusion_line.">";

        } else {
            $juli_fusions{$sample_id}{"$GeneA-$GeneB"}{"$brkA-$brkB"} = $fusion_line;
            $juli_fusions_org{$sample_id}{"$GeneA-$GeneB"}{"$brkA-$brkB"} = $line;
            if ($SplitA >= 20 || $SplitB >= 20 || ($SplitA >= 10 && $SplitB >= 10)) {
                $juli_fusions2{$sample_id}{"$GeneA-$GeneB"}{"$brkA-$brkB"} = $fusion_line;
                $juli_fusions2_org{$sample_id}{"$GeneA-$GeneB"}{"$brkA-$brkB"} = $line;
            }

            if ($SplitA >= 50 || $SplitB >= 50) {
                $juli_fusions3{$sample_id}{"$GeneA-$GeneB"}{"$brkA-$brkB"} = $fusion_line;
                $juli_fusions3_org{$sample_id}{"$GeneA-$GeneB"}{"$brkA-$brkB"} = $line;
            }
        }
    }
    close(INPUT2);
}

#parse GeneFuse outputs
my (%genefuse_fusions, %genefuse_fusions2, %genefuse_fusions3);
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

        my ($geneA, $geneB, $brkA, $brkB, $total_reads, $unique_reads);
        if ($line =~ /#Fusion:\s(.+)_ENST.+(chr.+)___(.+)_ENST.+(chr.+)\s\s\(total:\s(\d*),\sunique:(\d*)\)$/) {
            $geneA = $1;
            $geneB = $3;
            $brkA  = $2;
            $brkB  = $4;
            $total_reads = $5;
            $unique_reads= $6;

            if ($total_reads < 10 && !$target_genes{$geneA} && !$target_genes{$geneB}) {
                next;
            } elsif ($total_reads < 5) {
                next;
            } elsif ($artifacts{"$geneA-$geneB"} || $artifacts{"$geneB-$geneA"}) {
                #$line = "(".$line.")";
            
            } else {
                if ($target_genes{$geneA} || $target_genes{$geneB}) {
                 
                    $genefuse_fusions{$sample_id}{"$geneA-$geneB"}{"$brkA-$brkB"} = "*".$line;

                    if ($total_reads >= 10) {
                        $genefuse_fusions2{$sample_id}{"$geneA-$geneB"}{"$brkA-$brkB"} = "*".$line;
                    } 

                    if ($total_reads >= 50) {
                        $genefuse_fusions3{$sample_id}{"$geneA-$geneB"}{"$brkA-$brkB"} = "*".$line;
                    }
                } else {
                    $genefuse_fusions{$sample_id}{"$geneA-$geneB"}{"$brkA-$brkB"} = $line;

                    if ($total_reads >= 10) {
                        $genefuse_fusions2{$sample_id}{"$geneA-$geneB"}{"$brkA-$brkB"} = $line;
                    }

                    if ($total_reads >= 50) {
                        $genefuse_fusions3{$sample_id}{"$geneA-$geneB"}{"$brkA-$brkB"} = "*".$line;
                    }
                }
            }
        }
    }
    close(INPUT3);
}

# output possible fusions in file Possible-fusions.txt
my $sample_num = 0;
my $note;
$note = "Note: A row contianed in <> indicates that the fusion gene or partner is flanking of an element. ";
$note.= "A row with * marks containg a targeted fusion gene.";
print OUTPUT "$note\n\n";
foreach my $sample (sort {$a cmp $b} keys %samples) {
    $sample_num++;

    if ($sample_num == 1) {
        print OUTPUT "\n------ Sample $sample ------\n\n";
    } else {
        print OUTPUT "\n\n\n------ Sample $sample ------\n\n";
    } 

    print OUTPUT "Method FACTERA:\n\n";

    my $f_fusions = 0;
    foreach my $fusion (sort {$a cmp $b} keys %{$factera_fusions{$sample}}) {
        $f_fusions++;
        if ($f_fusions == 1) {
            print OUTPUT "Est_Type\tGene A\tGene B\tBreak A\tBreak B\tSupport A\tSupport B\tDiscord_mates\tTotal_Depth\n";
        }

        foreach my $brk (sort {$a cmp $b} keys %{$factera_fusions{$sample}{$fusion}}) {

           my $fusion_line = $factera_fusions{$sample}{$fusion}{$brk};    
           print OUTPUT "$fusion_line\n";

        } 
    }
    print OUTPUT "No fusion detected\n" if $f_fusions == 0;
    print OUTPUT "\n\nMethod JULI:\n\n";
    
    my $j_fusions = 0;
    foreach my $fusion (sort {$a cmp $b} keys %{$juli_fusions{$sample}}) {
        $j_fusions++;
        if ($j_fusions == 1) {

            print OUTPUT "Event\tGene A\tGene B\tBreak A\tBreak B\tDiscord A\tSplit A\tDirection\tDiscord B\tSplit B\tCosmic\n";
        }

        foreach my $brk (sort {$a cmp $b} keys %{$juli_fusions{$sample}{$fusion}}) {
           my $fusion_line2 = $juli_fusions{$sample}{$fusion}{$brk};
           print OUTPUT "$fusion_line2\n";
        }
    }
    print OUTPUT "No fusion detected\n" if $j_fusions == 0;
    print OUTPUT "\n\nMethod GeneFuse:\n\n";
    
    my $gfusions = 0;
    foreach my $fusion (sort {$a cmp $b} keys %{$genefuse_fusions{$sample}}) {
        $gfusions++;

        foreach my $brk (sort {$a cmp $b} keys %{$genefuse_fusions{$sample}{$fusion}}) {
           my $fusion_line3 = $genefuse_fusions{$sample}{$fusion}{$brk};
           print OUTPUT "$fusion_line3\n";
        }
    }
    print OUTPUT "No fusion detected\n" if $gfusions == 0;
}

# select reportable fusions and store them in file Reportable-fusions.txt
$sample_num = 0;
my (%selected_sample, %selected_f_fusion, %selected_j_fusion);

my $fusion_num = 0;
print OUTPUT2 "$note\n\n";
foreach my $sample (sort {$a cmp $b} keys %samples) {
    $sample_num++;

    my ($juli_fusions, $genefuse_fusions); 

    foreach my $fusion (sort {$a cmp $b} keys %{$juli_fusions2{$sample}}) {
        my ($GeneA, $GeneB) = split(/-/, $fusion);

        if ($GeneA =~ /Compl_/ && $GeneB =~ /Compl_/) {

            my ($geneA) = $GeneA =~ /^Compl_(.+)$/;
            my ($geneB) = $GeneB =~ /^Compl_(.+)$/;
            $juli_fusions .= "$geneA-$geneB ";
        } else {
            $juli_fusions .= $fusion." ";
        }
    }
    foreach my $fusion (sort {$a cmp $b} keys %{$genefuse_fusions2{$sample}}) {
        
        $genefuse_fusions .= $fusion." ";
    }

    # first, if FACTERA identified, to see if JuLI or GeneFuse also identified
    my %f_fusions;
    my $f_fusion_num = 0;
    my $pre_fusion;
    foreach my $fusion (sort {$a cmp $b} keys %{$factera_fusions2{$sample}}) {
        $f_fusions{$fusion} = 1;
        my ($gene_A, $gene_B) = split(/-/, $fusion);
        my $fusion2 = "$gene_B-$gene_A";

        next if $fusion2 eq $pre_fusion;

        if ($juli_fusions =~ /$fusion/ || $juli_fusions =~ /$fusion2/ || 
            $genefuse_fusions =~ /$fusion/ || $genefuse_fusions =~ /$fusion2/) {
            $f_fusion_num++;
            $fusion_num++;

            if ($f_fusion_num == 1) {

                if ($sample_num == 1) {
                    print OUTPUT2 "------ Sample $sample ------\n";
                } else {
                    print OUTPUT2 "\n\n------ Sample $sample ------\n";
                }
            }

            foreach my $brk (sort {$a cmp $b} keys %{$factera_fusions2{$sample}{$fusion}}) {

                if ($sample_num == 1) {
                    print OUTPUT2 "Method FACTERA:\n";
                } else {
                    print OUTPUT2 "\nMethod FACTERA:\n";
                }
                print OUTPUT2 "Est_Type\tGene A\tGene B\tBreak A\tBreak B\tSupport A\tSupport B\tDiscord_Mates\tTotal_Depth\n";
                #$selected_f_fusion{$sample}{"title"} = $factera_fusions_title;           
            
                my $fusion_line = $factera_fusions2{$sample}{$fusion}{$brk};
                print OUTPUT2 "$fusion_line\n";
                my $fusion_line_org = $factera_fusions2_org{$sample}{$fusion}{$brk};
                $selected_f_fusion{$sample}{$fusion}{$brk} = $fusion_line_org; 

                foreach my $brk2 (sort {$a cmp $b} keys %{$factera_fusions2{$sample}{$fusion2}}) {

                   my $fusion_line2 = $factera_fusions2{$sample}{$fusion2}{$brk2};
                   $fusion_line2 = $factera_fusions{$sample}{$fusion2}{$brk2} if !$fusion_line2;
                   print OUTPUT2 "$fusion_line2\n" if $fusion_line2;
                
                   my $fusion_line2_org = $factera_fusions2_org{$sample}{$fusion2}{$brk2};
                   $fusion_line2_org = $factera_fusions_org{$sample}{$fusion2}{$brk2} if !$fusion_line2_org;
                   $selected_f_fusion{$sample}{$fusion2}{$brk} = $fusion_line2_org if $fusion_line2_org;
               }
            }

            # parse JuLI output
            my $fus;
            if ($juli_fusions2{$sample}{$fusion}) {
                $fus = $fusion;
            } elsif ($juli_fusions2{$sample}{$fusion2}) {
                $fus = $fusion2;
            } elsif ($juli_fusions2{$sample}{"Compl_$gene_A-Compl_$gene_B"}) {
                $fus = "Compl_$gene_A-Compl_$gene_B";
            } elsif ($juli_fusions2{$sample}{"$gene_A-Compl_$gene_B"}) {
                $fus = "$gene_A-Compl_$gene_B";
            } elsif ($juli_fusions2{$sample}{"Compl_$gene_A-$gene_B"}) {
                $fus = "Compl_$gene_A-$gene_B";
            }

            my %direction;
            foreach my $brk (sort {$a cmp $b} keys %{$juli_fusions2{$sample}{$fus}}) {
                my $fusion_line_j = $juli_fusions2{$sample}{$fusion}{$brk};
                my $fusion_line_j_org = $juli_fusions2_org{$sample}{$fusion}{$brk};
                my ($g_A, $g_B) = split(/-/, $fusion);
                #next if $geneAs{$g_A};
       
                my @list = split(/\t/, $fusion_line_j);
                my $dir = $list[7];
                next if $direction{$dir}; 

                my $fusion_line_j_compl = $juli_fusions2{$sample}{"Compl_$g_A-Compl_$g_B"}{$brk};
                my $fusion_line_j_compl2 = $juli_fusions2{$sample}{"$g_A-Compl_$g_B"}{$brk};
                my $fusion_line_j_compl3 = $juli_fusions2{$sample}{"Compl_$g_A-$g_B"}{$brk};

                my ($fusion_line_j2, $fusion_line_j2_compl, $fusion_line_j2_compl2, $fusion_line_j2_compl3);
                my $fusion_line_j2_org;
                
                foreach my $brk2 (sort {$a cmp $b} keys %{$juli_fusions2{$sample}{$fusion2}}) {
                   $fusion_line_j2 = $juli_fusions2{$sample}{$fusion2}{$brk2};
                   $fusion_line_j2_org = $juli_fusions2_org{$sample}{$fusion2}{$brk2};
                   $fusion_line_j2 = $juli_fusions{$sample}{$fusion2}{$brk2} if !$fusion_line_j2;
                   $fusion_line_j2_org = $juli_fusions_org{$sample}{$fusion2}{$brk2} if !$fusion_line_j2_org;
                   $fusion_line_j2_compl = $juli_fusions2{$sample}{"Compl_$g_B-Compl_$g_A"}{$brk2};
                   $fusion_line_j2_compl = $juli_fusions{$sample}{"Compl_$g_B-Compl_$g_A"}{$brk2} if !$fusion_line_j2_compl;
                   $fusion_line_j2_compl2 = $juli_fusions2{$sample}{"$g_B-Compl_$g_A"}{$brk2} if !$fusion_line_j2_compl;
                   $fusion_line_j2_compl2 = $juli_fusions{$sample}{"$g_B-Compl_$g_A"}{$brk2} if !$fusion_line_j2_compl2;
                   $fusion_line_j2_compl3 = $juli_fusions2{$sample}{"Compl_$g_B-$g_A"}{$brk2} if !$fusion_line_j2_compl;
                   $fusion_line_j2_compl3 = $juli_fusions{$sample}{"Compl_$g_B-$g_A"}{$brk2} if !$fusion_line_j2_compl3;
                }

                if ($fusion_line_j || $fusion_line_j2 || 
                    $fusion_line_j_compl || $fusion_line_j_compl2 || $fusion_line_j_compl3 || 
                    $fusion_line_j2_compl || $fusion_line_j2_compl2 || $fusion_line_j2_compl3) {
                    print OUTPUT2 "\nMethod JULI:\n";
                    print OUTPUT2 "Event\tGene A\tGene B\tBreak A\tBreak B\tDiscord A\tSplit A\tDirection\tDiscord B\tSplit B\tCosmic\n";
                    #$selected_j_fusion{$sample}{"title"} = $juli_fusions_title;
                    print OUTPUT2 "$fusion_line_j\n" if $fusion_line_j;
                    $selected_j_fusion{$sample}{$fusion}{$brk} = $fusion_line_j_org if $fusion_line_j_org;
                    print OUTPUT2 "$fusion_line_j_compl\n" if $fusion_line_j_compl;

                    print OUTPUT2 "$fusion_line_j_compl2\n" if $fusion_line_j_compl2;
                    print OUTPUT2 "$fusion_line_j_compl3\n" if $fusion_line_j_compl3;
                    print OUTPUT2 "$fusion_line_j2\n" if $fusion_line_j2;
                    $selected_j_fusion{$sample}{$fusion2}{$brk} = $fusion_line_j2_org if ($fusion_line_j2_org && !$fusion_line_j_org);
                    print OUTPUT2 "$fusion_line_j2_compl\n" if $fusion_line_j2_compl;
                    print OUTPUT2 "$fusion_line_j2_compl2\n" if $fusion_line_j2_compl2;
                    print OUTPUT2 "$fusion_line_j2_compl3\n" if $fusion_line_j2_compl3;
                    $direction{$dir} = 1;
                }
            }

            my $breakpoint_count = 0;
            foreach my $brk (sort {$a cmp $b} keys %{$genefuse_fusions2{$sample}{$fusion}}) {

                my $fusion_line_g = $genefuse_fusions2{$sample}{$fusion}{$brk};

                my $fusion_line_g2;

                foreach my $brk2 (sort {$a cmp $b} keys %{$genefuse_fusions2{$sample}{$fusion2}}) {
                   $fusion_line_g2 = $genefuse_fusions2{$sample}{$fusion2}{$brk2};
                   $fusion_line_g2 = $genefuse_fusions{$sample}{$fusion2}{$brk2} if !$fusion_line_g2;
                }

                if ($fusion_line_g || $fusion_line_g2) {
                    $breakpoint_count++;
                    print OUTPUT2 "\nMethod GeneFuse:\n" if $breakpoint_count == 1;
                    print OUTPUT2 "$fusion_line_g\n" if $fusion_line_g;
                    print OUTPUT2 "$fusion_line_g2\n" if $fusion_line_g2;
                }
            } 
            $selected_sample{$sample} = 1;            
            print OUTPUT2 "................ End of last fusion ................\n";
        } else {

            foreach my $brk (sort {$a cmp $b} keys %{$factera_fusions3{$sample}{$fusion}}) {
                if ($factera_fusions3{$sample}{$fusion}{$brk} && $factera_fusions3{$sample}{$fusion2}{$brk})  { # reciprocal
                    my $fusion_line = $factera_fusions3{$sample}{$fusion}{$brk};
                    my $fusion_line2 = $factera_fusions3{$sample}{$fusion2}{$brk};
                    $fusion_line2 = $factera_fusions2{$sample}{$fusion2}{$brk} if !$fusion_line2;
                    $fusion_line2 = $factera_fusions{$sample}{$fusion2}{$brk} if !$fusion_line2;
                    if ($known_fusions{$fusion} || $known_fusions{$fusion2}) {
                        $f_fusion_num++;
                        $fusion_num++;
                        print OUTPUT2 "Method FACTERA:\n";
                        print OUTPUT2 "Est_Type\tGene A\tGene B\tBreak A\tBreak B\tSupport A\tSupport B\tDiscord_mates\tTotal_Depth\n";
                        my $fusion_line = $factera_fusions3{$sample}{$fusion}{$brk};
                        print OUTPUT2 "$fusion_line\n";

                        foreach my $brk2 (sort {$a cmp $b} keys %{$factera_fusions3{$sample}}){
                           $selected_f_fusion{$sample}{$fusion}{$brk} = $fusion_line;
                           my $fusion_line2 = $factera_fusions3{$sample}{$fusion2}{$brk2};
                           print OUTPUT2 "$fusion_line2\n" if $fusion_line2;
                           my $fusion_line2_org = $factera_fusions3_org{$sample}{$fusion2}{$brk2};
                           $selected_f_fusion{$sample}{$fusion2}{$brk} = $fusion_line2_org;
                       }
                    }
                    $selected_sample{$sample} = 1;
                    print OUTPUT2 "................ End of last fusion ................\n";
                }
            }
        }
        $pre_fusion = $fusion;      
    }     
    # second, if GeneFuse identified, to see if JuLI identified, too.
    my %g_fusions;
    my $gfusion_num = 0;
    my $pre_gfusion;
    foreach my $gfusion (sort {$a cmp $b} keys %{$genefuse_fusions2{$sample}}) {
        $g_fusions{$gfusion} = 1;
        my ($gene_A, $gene_B) = split(/-/, $gfusion);
        my $gfusion2 = "$gene_B-$gene_A";

        next if $gfusion2 eq $pre_gfusion;

        if (!$f_fusions{$gfusion} && !$f_fusions{$gfusion2}) {

            my $fus;
            if ($juli_fusions2{$sample}{$gfusion}) {
                $fus = $gfusion;
            } elsif ($juli_fusions2{$sample}{$gfusion2}) {
                $fus = $gfusion2;
            } elsif ($juli_fusions2{$sample}{"Compl_$gene_A-Compl_$gene_B"}) {
                $fus = "Compl_$gene_A-Compl_$gene_B";
            } elsif ($juli_fusions2{$sample}{"$gene_A-Compl_$gene_B"}) {
                $fus = "$gene_A-Compl_$gene_B";
            } elsif ($juli_fusions2{$sample}{"Compl_$gene_A-$gene_B"}) {
                $fus = "Compl_$gene_A-$gene_B";
            }
            
            my %direction;
            foreach my $brk (sort {$a cmp $b} keys %{$juli_fusions2{$sample}{$fus}}) {
               my $fusion_line_j = $juli_fusions2{$sample}{$gfusion}{$brk};

               my $fusion_line_j_org = $juli_fusions2_org{$sample}{$gfusion}{$brk};
               my ($g_A, $g_B) = split(/-/, $fus);

               my @list = split(/\t/, $fusion_line_j);
               my $dir = $list[7];
               next if $direction{$dir};

               my $fusion_line_j_compl = $juli_fusions2{$sample}{"Compl_$g_A-Compl_$g_B"}{$brk};
               my $fusion_line_j_compl2 = $juli_fusions2{$sample}{"$g_A-Compl_$g_B"}{$brk};
               $fusion_line_j_compl2 = $juli_fusions{$sample}{"$g_A-Compl_$g_B"}{$brk} if !$fusion_line_j_compl2;
               my $fusion_line_j_compl3 = $juli_fusions2{$sample}{"Compl_$g_A-$g_B"}{$brk};
               $fusion_line_j_compl3 = $juli_fusions{$sample}{"Compl_$g_A-$g_B"}{$brk} if !$fusion_line_j_compl3;
            
               my $fusion_line_j2 = $juli_fusions2{$sample}{$gfusion2}{$brk};
               my $fusion_line_j2_org = $juli_fusions2_org{$sample}{$gfusion2}{$brk};
               $fusion_line_j2 = $juli_fusions{$sample}{$gfusion2}{$brk} if !$fusion_line_j2;
               $fusion_line_j2_org = $juli_fusions_org{$sample}{$gfusion2}{$brk} if !$fusion_line_j2_org;
               my $fusion_line_j2_compl = $juli_fusions2{$sample}{"Compl_$g_B-Compl_$g_A"}{$brk};
               $fusion_line_j2_compl = $juli_fusions{$sample}{"Compl_$g_B-Compl_$g_A"}{$brk} if !$fusion_line_j2_compl;
               my $fusion_line_j2_compl2 = $juli_fusions2{$sample}{"$g_B-Compl_$g_A"}{$brk} if !$fusion_line_j2_compl;
               $fusion_line_j2_compl2 = $juli_fusions{$sample}{"$g_B-Compl_$g_A"}{$brk} if !$fusion_line_j2_compl2;
               my $fusion_line_j2_compl3 = $juli_fusions2{$sample}{"Compl_$g_B-$g_A"}{$brk} if !$fusion_line_j2_compl;
               $fusion_line_j2_compl3 = $juli_fusions{$sample}{"Compl_$g_B-$g_A"}{$brk} if !$fusion_line_j2_compl3;
            
               if ($fusion_line_j || $fusion_line_j2 ||
                   $fusion_line_j_compl || $fusion_line_j_compl2 || $fusion_line_j_compl3 ||
                   $fusion_line_j2_compl || $fusion_line_j2_compl2 || $fusion_line_j2_compl3) {
                   $gfusion_num++;
                   $fusion_num++;

                   if (!$selected_sample{$sample}) {
                       if ($gfusion_num == 1) {

                           if ($sample_num == 1) {
                               print OUTPUT2 "------ Sample $sample ------\n";
                           } else {
                               print OUTPUT2 "\n------ Sample $sample ------\n";
                           }
                       }
                   }
                   if ($sample_num == 1) {
                        print OUTPUT2 "Method JULI:\n";
                   } else {
                        print OUTPUT2 "\nMethod JULI:\n";
                   }
                   print OUTPUT2 "Event\tGene A\tGene B\tBreak A\tBreak B\tDiscord A\tSplit A\tDirection\tDiscord B\tSplit B\tCosmic\n";
                   print OUTPUT2 "$fusion_line_j\n";
                   $selected_j_fusion{$sample}{$gfusion}{$brk} = $fusion_line_j_org if $fusion_line_j_org;

                   print OUTPUT2 "$fusion_line_j_compl\n" if $fusion_line_j_compl;
                   print OUTPUT2 "$fusion_line_j_compl2\n" if $fusion_line_j_compl2;
                   print OUTPUT2 "$fusion_line_j_compl3\n" if $fusion_line_j_compl3;

                   print OUTPUT2 "$fusion_line_j2\n" if $fusion_line_j2;
                   $selected_j_fusion{$sample}{$gfusion2}{$brk} = $fusion_line_j2_org if $fusion_line_j2_org;            
                
                   print OUTPUT2 "$fusion_line_j2_compl\n" if $fusion_line_j2_compl;
                   print OUTPUT2 "$fusion_line_j2_compl2\n" if $fusion_line_j2_compl2;
                   print OUTPUT2 "$fusion_line_j2_compl3\n" if $fusion_line_j2_compl3;
                   $direction{$dir} = 1;

                   my $breakpoint_count = 0;
                   foreach my $brk (sort {$a cmp $b} keys %{$genefuse_fusions2{$sample}{$gfusion}}) {
                       $breakpoint_count++;
                       my $fusion_line_g = $genefuse_fusions2{$sample}{$gfusion}{$brk};

                       my $fusion_line_g2;

                       foreach my $brk2 (sort {$a cmp $b} keys %{$genefuse_fusions2{$sample}{$gfusion2}}) {
                          $fusion_line_g2 = $genefuse_fusions2{$sample}{$gfusion2}{$brk2};
                          $fusion_line_g2 = $genefuse_fusions{$sample}{$gfusion2}{$brk2} if !$fusion_line_g2;
                       }
                       print OUTPUT2 "\nMethod GeneFuse:\n" if $breakpoint_count == 1;
                       print OUTPUT2 "$fusion_line_g\n";
                       print OUTPUT2 "$fusion_line_g2\n" if $fusion_line_g2;
                       $selected_sample{$sample} = 1;
                   }
                   print OUTPUT2 "................ End of last fusion ................\n";
               }
            }
        } else {
            my $breakpoint_count = 0;
            foreach my $brk (sort {$a cmp $b} keys %{$genefuse_fusions3{$sample}{$gfusion}}) {
               if ($genefuse_fusions3{$sample}{$gfusion}{$brk} && $genefuse_fusions3{$sample}{$gfusion2}{$brk}) { # reciprocal
  
                   if ($known_fusions{$gfusion} || $known_fusions{$gfusion2}) {
                       $fusion_num++;
                       $breakpoint_count++;
                       print OUTPUT2 "\nMethod GeneFuse:\n" if $breakpoint_count == 1;
                       print OUTPUT2 "\n$genefuse_fusions3{$sample}{$gfusion}{$brk}\n";
                       foreach my $brk2 (sort {$a cmp $b} keys %{$genefuse_fusions3{$sample}{$gfusion2}}) {
                          print OUTPUT2 "$genefuse_fusions3{$sample}{$gfusion2}{$brk2}\n"
                          if $genefuse_fusions3{$sample}{$gfusion2}{$brk2};
                       }                       
                   }
                   $selected_sample{$sample} = 1;
                   print OUTPUT2 "................ End of last fusion ................\n";
               }
            }
        }
        $pre_gfusion = $gfusion;
    }

    # last, check if JULI identified solely
    my $jfusion_num = 0;
    my $pre_jfusion;
    foreach my $jfusion (sort {$a cmp $b} keys %{$juli_fusions3{$sample}}) {
        
        my ($gene_A, $gene_B) = split(/-/, $jfusion);
        my $jfusion2 = "$gene_B-$gene_A";

        next if $jfusion2 eq $pre_jfusion;
   
        my ($geneA, $geneB) = ($gene_A, $gene_B);
        if ($gene_A =~ /Compl_(.+)$/) {
            $geneA = $1;
        }

        if ($gene_B =~ /Compl_(.+)$/) {
            $geneB = $1;
        }     

        foreach my $brk (sort {$a cmp $b} keys %{$juli_fusions3{$sample}{$jfusion}}) { 
            my $fusion_line_j = $juli_fusions3{$sample}{$jfusion}{$brk};
            my $fusion_line_j_org = $juli_fusions3_org{$sample}{$jfusion}{$brk};

            my ($fusion_line_j2, $fusion_line_j2_org);
            foreach my $brk2 (sort {$a cmp $b} keys %{$juli_fusions3{$sample}{$jfusion2}}) {
               $fusion_line_j2 = $juli_fusions3{$sample}{$jfusion2}{$brk};
               $fusion_line_j2_org = $juli_fusions3_org{$sample}{$jfusion2}{$brk};
               $fusion_line_j2 = $juli_fusions3{$sample}{$jfusion2}{$brk} if !$fusion_line_j2;
               $fusion_line_j2_org = $juli_fusions3_org{$sample}{$jfusion2}{$brk} if !$fusion_line_j2_org;
            }

            # if it is reciprocal
            if ($fusion_line_j2) { 
                if (!$f_fusions{"$geneA-$geneB"} && !$f_fusions{"$geneB-$geneA"} && 
                    !$g_fusions{"$geneA-$geneB"} && !$g_fusions{"$geneB-$geneA"} ) {

                    # if it can be found in the "known-fusions" database
                    if ($known_fusions{$jfusion} || $known_fusions{$jfusion2}) {
                        $jfusion_num++;
                        $fusion_num++;
                        if (!$selected_sample{$sample}) {
                            if ($jfusion_num == 1) {

                                if ($sample_num == 1) {
                                    print OUTPUT2 "------ Sample $sample ------\n";
                                } else {
                                    print OUTPUT2 "\n------ Sample $sample ------\n";
                                }
                            }
                        }
                        print OUTPUT2 "\nMethod JuLI:\n";
                        print OUTPUT2 "Event\tGene A\tGene B\tBreak A\tBreak B\tDiscord A\tSplit A\tDirection\tDiscord B\tSplit B\tCosmic\n";
                        print OUTPUT2 "\n$fusion_line_j\n";
                        $selected_j_fusion{$sample}{$jfusion}{$brk} = $fusion_line_j_org if $fusion_line_j_org;

                        print OUTPUT2 "\n$fusion_line_j2\n";
                        $selected_j_fusion{$sample}{$jfusion2}{$brk} = $fusion_line_j2_org if $fusion_line_j2_org;
                      
                        $selected_sample{$sample} = 1;
                        print OUTPUT2 "................ End of last fusion ................\n";
                    }
                }
            }
        }
    }
}
if ($fusion_num == 0) {
    print OUTPUT2 "\nNo reportable fusion detected!\n";
}

# generate TSV files in Arriba output format for draw_fusion.R program as input
system("mkdir $outpath/tmp-dir");
my %selected_samples; 
foreach my $sample (sort {$a cmp $b} keys %selected_f_fusion) {
   $selected_samples{$sample} = 1;
   open(OUTPUT22, ">$outpath/tmp-dir/$sample.tsv") || die "cannot open $sample.tsv";
 
   my $arriba_tsv_title = "#gene1\tgene2\tstrand1(gene/fusion)\tstrand2(gene/fusion)\tbreakpoint1\tbreakpoint2\tsite1\tsite2\ttype\tsplit_reads1\tsplit_reads2\tdiscordant_mates\tcoverage1\tcoverage2\tconfidence\treading_frame\ttags\tretained_protein_domains\tclosest_genomic_breakpoint1\tclosest_genomic_breakpoint2\tgene_id1]tgene_id2\ttranscript_id1\ttranscript_id2\tdirection1\tdirection2\tfilters\tfusion_transcript\tpeptide_sequence\tread_identifiers"; 
   print OUTPUT22 "$arriba_tsv_title\n";
   my %selected_fusions; 
   foreach my $fusion (sort {$a cmp $b} keys %{$selected_f_fusion{$sample}}) {
       $selected_fusions{$fusion} = 1;   

       foreach my $brk (sort {$a cmp $b} keys %{$selected_f_fusion{$sample}{$fusion}}) {
           my $f_funsion_line = $selected_f_fusion{$sample}{$fusion}{$brk};

           my @f_items=split(/\t/, $f_funsion_line); 
           my $Est_Type = $f_items[0];
           if ($Est_Type eq "TRA") {
               $Est_Type = "translocation";
           } elsif ($Est_Type eq "INV") {
               $Est_Type = "inversion";
           } elsif ($Est_Type eq "DEL") {
               $Est_Type = "deletion";
           } elsif ($Est_Type eq "DUP") {
               $Est_Type = "duplication";
           } elsif ($Est_Type eq "-") {
               $Est_Type = "not determined";
           }

           my $Region1  = $f_items[1];
           my $Region2  = $f_items[2];
           my $Break1   = $f_items[3];
           my $Break2   = $f_items[4];
           my $Support1 = $f_items[5];
           my $Support2 = $f_items[6];
           my $Offset   = $f_items[7];
           my $Orient   = $f_items[8];
           my ($str1, $str2) = split(/\s/, $Orient);
           my ($S1) = $str1 =~ /^\d(.)$/;
           my $Strand1 = "$S1/$S1";
           my ($S2) = $str2 =~ /^\d(.)$/;
           my $Strand2 = "$S2/$S2";
           my $Order1   = $f_items[9];
           if ($Order1 eq "NC") {
              $Order1 = "downstream";
           } else {
              $Order1 = "upstream";
           }
           my $Order2   = $f_items[10];
           if ($Order2 eq "NC") {
               $Order2 = "downstream";
           } else {
               $Order2 = "upstream";
           }
           my $Break_depth = $f_items[11];

           my $confidence;
           if ($Break_depth >= 50) {
              $confidence = "high";
           } elsif ($Break_depth >= 20) {
              $confidence = "medium";
           } else {
              $confidence = "low";
           }
           my $discordant_mates = $f_items[14];
           my $Total_depth = $f_items[16];
           my $Fusion_seq  = $f_items[17];

           my $arriba_line = "$Region1\t$Region2\t$Strand1\t$Strand2\t$Break1\t$Break2\t";
           $arriba_line .= ".\t.\t$Est_Type\t$Support1\t$Support2\t$discordant_mates\t";
           $arriba_line .= ".\t.\t$confidence\t.\t.\t.\t.\t.\t.\t.\t.\t$Order1\t$Order2\t.\t$Fusion_seq\t.\t.";
           print OUTPUT22 "$arriba_line\n";
       }
   }

   foreach my $fusion (sort {$a cmp $b} keys %{$selected_j_fusion{$sample}}) {

       next if ($selected_fusions{$fusion});
       Juli_to_Arriba_fusion_line($sample, $fusion);
   }
   close(OUTPUT22);
}

foreach my $sample (sort {$a cmp $b} keys %selected_j_fusion) {

   next if $selected_samples{$sample};
   open(OUTPUT22, ">$outpath/tmp-dir/$sample.tsv") || die "cannot open $sample.tsv";

   my $arriba_tsv_title = "#gene1\tgene2\tstrand1(gene/fusion)\tstrand2(gene/fusion)\tbreakpoint1\tbreakpoint2\tsite1\tsite2\ttype\tsplit_reads1\tsplit_reads2\tdiscordant_mates\tcoverage1\tcoverage2\tconfidence\treading_frame\ttags\tretained_protein_domains\tclosest_genomic_breakpoint1\tclosest_genomic_breakpoint2\tgene_id1]tgene_id2\ttranscript_id1\ttranscript_id2\tdirection1\tdirection2\tfilters\tfusion_transcript\tpeptide_sequence\tread_identifiers";
   print OUTPUT22 "$arriba_tsv_title\n";
   foreach my $fusion (sort {$a cmp $b} keys %{$selected_j_fusion{$sample}}) {
      Juli_to_Arriba_fusion_line($sample, $fusion);
   }
   close(OUTPUT22);
}

#draw annotation of each fusion by calling draw_fusions.R
my @files = <$outpath/tmp-dir/*.tsv>;
my %files = map {$_=>1} @files if @files;
foreach my $file (sort {$a cmp $b} keys %files) {
   my ($sample_id) = $file =~ /^.+\/(.+)\.tsv$/;

   system("nohup ./draw_fusions.R --fusions='$outpath/tmp-dir/$sample_id.tsv' --alignments='$bampath/$sample_id.bam' --output='$outpath/$sample_id-fusions.pdf' --annotation=database/GENCODE19.gtf --cytobands=database/cytobands_hg19_hs37d5_GRCh37_v2.1.0.tsv --proteinDomains=database/protein_domains_hg19_hs37d5_GRCh37_v2.1.0.gff3 &"); 

}

sub Juli_to_Arriba_fusion_line {

   my ($sample, $fusion) = @_;
   foreach my $brk (sort {$a cmp $b} keys %{$selected_j_fusion{$sample}{$fusion}}) {
      my $j_funsion_line = $selected_j_fusion{$sample}{$fusion}{$brk};

      my @items=split(/\t/, $j_funsion_line);
      my $ChrA = $items[0];
      my $BreakA = "$ChrA:$items[1]";
      my $OriA = $items[2];
      my $direction1;
      if ($OriA == 1) {
          $OriA = "-";
          $direction1 = "upstream";
      } elsif ($OriA == 0) {
          $OriA = "+";
          $direction1 = "downstream";
      }
      my $DisA = $items[3];
      my $SplitA = $items[4];
      my $ChrB = $items[5];
      my $BreakB = "$ChrB:$items[6]";
      my $OriB = $items[7];
      my $direction2;
      if ($OriB == 1) {
          $OriB = "-";
          $direction2 = "upstream";
      } elsif ($OriB == 0) {
          $OriB = "+";
          $direction2 = "downstream";
      }
      my $DisB = $items[8];
      my $SplitB = $items[9];
      my $Event= $items[10];
      my $GeneA= $items[11];
      my $StrGeneA = $items[12];
      my $Strand1 = "$StrGeneA/$OriA";
      my $GeneB=$items[13];
      my $StrGeneB = $items[14];
      my $Strand2 = "$StrGeneB/$OriB";
      my $InfoA=$items[15];
      my $InfoB=$items[16];
      my $Direction= $items[17];
      my $Frame=$items[18];

      my $discordant_mates = $DisA + $DisB;
      my $Break_depth = $SplitA + $SplitB;
      my $confidence;
      if ($Break_depth >= 50) {
          $confidence = "high";
      } elsif ($Break_depth >= 20) {
          $confidence = "medium";
      } else {
          $confidence = "low";
      }
      my $arriba_line = "$GeneA\t$GeneB\t$Strand1\t$Strand2\t$BreakA\t$BreakB\t";
      $arriba_line .= "$InfoA\t$InfoA\t$Event\t$SplitA\t$SplitB\t$discordant_mates\t";
      $arriba_line .= ".\t.\t$confidence\t$Frame\t.\t.\t.\t.\t.\t.\t.\t$direction1\t$direction2\t.\t.\t.\t.";
      print OUTPUT22 "$arriba_line\n";
   }
}
