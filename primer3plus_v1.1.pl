#!/usr/bin/perl
use strict;
use warnings;
use Bio::Tools::Run::Primer3;
use Bio::SeqIO;

#AUTHOR: Rebecca Allan
#DATE: Mar 3/15
#PURPOSE: Produce primers from an inputted DNA sequence (fasta file) using primer3 tool
#NOTES: Must install BioPerl first following the instructions at http://www.bioperl.org/wiki/Installing_BioPerl_on_Windows
        #Must download Bio/Tools/Run and Bio/Roots/Roots modules and place in Perl>lib for primer3 tool to be accessed
        #Must copy primer3_core.exe in same location as THIS perl file to avoid "primer3 cannot be found" error
        #Use primer3 version 1.1.4 to avoid "missing SEQUENCE tag" error
        #Copy primer3_config directory in same location as THIS perl file to avoid "thermodynamic approach" error
        #Check output file (temp.out) to see actual left/right primer sequences

#Input a fasta file (Primer3Test.fa)
my $seqio = Bio::SeqIO->new(-file=>'Primer3Test.fa');
my $seq = $seqio->next_seq;
my $primer3 = Bio::Tools::Run::Primer3->new(-seq => $seq,
                                            -outfile => "temp.out",
                                            -path => "/usr/bin/primer3_core");

unless ($primer3->executable) {
    print STDERR "primer3 can not be found. Is it installed?\n";
    exit(-1)
}

#displays all arguments within the primer3 tool (big amount of text)
#my $args = $primer3->arguments;
#print "ARGUMENT\tMEANING\n";
#foreach my $key (keys %{$args}) {print "$key\t", $$args{$key}, "\n"}

#set the maximum and minimum Tm of the primer
$primer3->add_targets('PRIMER_MIN_TM'=>56, 'PRIMER_MAX_TM'=>70);

#design the primers. This runs primer3 and returns a Bio::Tools::Run::Primer3 object with the results
my $results = $primer3->run;

#see the Bio::Tools::Run::Primer3 pod for things that you can get from this
#displays number of left/right primers produced (1 primer = 1 left and 1 right primer)
print "There were ", $results->number_of_results, " primers\n";
