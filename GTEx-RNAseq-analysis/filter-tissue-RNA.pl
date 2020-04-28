#!/usr/bin/perl
use strict;
use warnings;
my $fh;

my @wanttissues = ("Lung","Colon","Kidney","Liver","Small Intestine");
my %tissuehash;
foreach my $t (@wanttissues){
    $tissuehash{$t} = 1;
}

my %clinhash;
open ($fh,"<sample-info.txt") || die "no";
<$fh>;
while (<$fh>){
    s/[\r\n]//g;
    s/"//g;
    my $line = $_;
    my @lines = split("\t",$line);
    my $tissue = $lines[1];

    if (defined $tissuehash{$tissue}){
        $clinhash{$lines[0]} = 1;
    }
}
close $fh;


open ($fh,"<GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct") || die "c";

#change this to the TPM file and rerun
#open ($fh,"<GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct") || die "c";
<$fh>;
<$fh>;
my $head = <$fh>;
$head =~ s/[\r\n]//g;
my @headers = split ("\t",$head);

print "Name";
my @wantcols;
for (my $i = 2; $i < scalar @headers; $i++){
    if (defined $clinhash{$headers[$i]}){
        print "\t$headers[$i]";
        push (@wantcols, $i);
    }
}
print "\n";

while (<$fh>){
    s/[\r\n]//g;
    s/"//g;
    my $line = $_;
    my @lines = split("\t",$line);
    print "$lines[1]";
    foreach my $i (@wantcols){
        print "\t$lines[$i]";
    }
    print "\n";
}
close $fh;
