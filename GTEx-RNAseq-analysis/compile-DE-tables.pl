#!/usr/bin/perl
use strict;
use warnings;
my $fh;
my $ofh;
my @comps = ('3-2','4-3','5-4','6-5','7-6');
my %hash;
my %genes;

open ($ofh ,">GTEx-lung-DEseq2-results.LRT.timeBin.combined.SexSmokingIschemia.all.txt") || die "c";
foreach my $f (@comps){
    open ($fh,"<GTEx-lung-DEseq2-results.LRT.timeBin.$f.SexSmokingIschemia.txt") || die "c";
    <$fh>;
    my %temphash;
    while (<$fh>){
        s/[\r\n]//g;
        s/"//g;
        my @lines = split ("\t",$_);
        my $gene = $lines[0];
        my $basemean = $lines[1];
        my $lfc = $lines[2];
        my $se = $lines[3];
        my $stat = $lines[4];
        my $pval = $lines[5];
        my $adjp = $lines[6];
        if ($adjp ne "NA"){
            #if ($adjp < 0.05){
                $temphash{$gene} = [$lfc,$se];
                if ($f eq '3-2'){
                    $genes{$gene} = [$basemean,$stat,$pval,$adjp];
                }
            #}
        }
    }
    close $fh;

    $hash{$f} = \%temphash;
}

print $ofh "gene\tbaseMean\t3.2_log2FC\t3.2_lfcSE\t4.3_log2FC\t4.3_lfcSE\t5.4_log2FC\t5.4_lfcSE\t6.5_log2FC\t6.5_lfcSE\t7.6_log2FC\t7.6_lfcSE\tstat\tpvalue\tpadj\n";
foreach my $g (sort {$genes{$a}[3] <=> $genes{$b}[3] } keys %genes){
    print $ofh "$g\t$genes{$g}[0]";
    foreach my $f (@comps){
        my @data = @{$hash{$f}{$g}};
        print $ofh "\t$data[0]\t$data[1]";
    }
    print $ofh "\t$genes{$g}[1]\t$genes{$g}[2]\t$genes{$g}[3]\n";
}
