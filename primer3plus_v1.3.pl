#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Tools::Run::Primer3;
use Bio::Tools::Run::StandAloneBlastPlus;

#AUTHOR: Team Dynamic
#DATE STARTED: Mar 3/15
#DATE UPDATED: Mar 4/15
#PURPOSE: Produce primers from an inputted DNA sequence (fasta file) using primer3 tool
#NOTES: Must install BioPerl first following the instructions at http://www.bioperl.org/wiki/Installing_BioPerl_on_Windows
        #Must download Bio/Tools/Run and Bio/Roots/Roots modules and place in Perl>lib for primer3 tool to be accessed
        #Must copy primer3_core.exe in same location as THIS perl file to avoid "primer3 cannot be found" error
        #Use primer3 version 1.1.4 to avoid "missing SEQUENCE tag" error
        #Alter -path => in Primer3->new() to avoid "SH: Command not found" error
        #Copy primer3_config directory in same location as THIS perl file to avoid "thermodynamic approach" error
        #Check output file (temp.out) to see actual left/right primer sequences

#*******************************************************************************************************#
#ATTACH JAVA GUI TO THIS PROGRAM FILE



#*******************************************************************************************************#
#USER-SELECTED ORGANISM FROM GUI USED AS INPUT



#*******************************************************************************************************#
#ACCESS ORGANISM'S FASTA FILE FROM DATABASE



#*******************************************************************************************************#
#PRIMER3 TOOL TO PRODUCE PRIMERS FROM FASTA FILE (Rebecca Allan, sourced from Chad Matsalla)
#Input a fasta file and create a new primer3 list
my $seqio = Bio::SeqIO->new(-file=>'Primer3Test.fa');
my $seq = $seqio->next_seq;
my $primer3 = Bio::Tools::Run::Primer3->new(-seq => $seq,
                                            -outfile => "temp.out",
                                            -path => "C:/Users/Rebecca/Desktop/IGP/PERL_PROGRAM/primer3_core");

#Test to see if primer3_core.exe is within the directory
unless ($primer3->executable) {
    print STDERR "Primer3 can not be found. Is it installed?\n";
    exit(-1)
}

#Display all arguments within the primer3 tool
#my $args = $primer3->arguments;
#print "ARGUMENT\tMEANING\n";
#foreach my $key (keys %{$args}) {
#    print "$key\t", $$args{$key}, "\n";
#}

#Adjust default values of specific arguments in primer3
$primer3->add_targets('PRIMER_MIN_TM'=>56, 'PRIMER_MAX_TM'=>65);
$primer3->add_targets('PRIMER_MIN_SIZE'=>20, 'PRIMER_MAX_SIZE'=>27, 'PRIMER_OPT_SIZE'=>20, 'PRIMER_DEFAULT_SIZE'=>20);
#Next line on primer product size range does not work, must be fixed
#$primer3->add_targets('PRIMER_PRODUCT_SIZE_RANGE'=>100..500);
#Number of primers produced
$primer3->add_targets('PRIMER_NUM_RETURN'=>5);
my $results = $primer3->run;

#Print the number of results (proves successful run)
print "There were ", $results->number_of_results, " primers\n";

#*******************************************************************************************************#
#RETRIEVE PRIMERS FROM TEMP.OUT FILE (Rebecca Allan and Phuong Ma)
#Array to contain temp.out contents
my @file_contents;
#Arrays to contain the primers and their positions
my @left_primers = ();
my @left_primers_pos = ();
my @right_primers = ();
my @right_primers_pos = ();
#Array to contain amplimers for each primer pair
my @amplimers = ();

#Read the entire temp.out file into an array and then close file
open (my $filehandle, "<", "temp.out") or die "Error reading file.\n";
foreach (<$filehandle>) {
    push @file_contents, $_;
}
close $filehandle;
chomp @file_contents;

#Scalar value $temp_out for array @file_contents for regex
my $temp_out = join '', @file_contents;

for (my $i = 0; $i < $results->number_of_results; $i++) {
    #Perform different regex statements depending on the index value of each primer pair element
    if ($i > 0) {
        $temp_out =~ /(PRIMER_LEFT_$i\_SEQUENCE=)(.+)(PRIMER_RIGHT_$i\_SEQUENCE=)(.+)(PRIMER_LEFT_$i=)(.+)(\,.+)(PRIMER_RIGHT_$i=)(.+?)(\,.+)/;
        #Find primers and their positions, then push them into the correct array
        my $temp_left_seq = $2;
        my $temp_left_pos = $6;
        my $temp_right_seq = $4;
        my $temp_right_pos = $9;
        push @left_primers, $temp_left_seq;
        push @left_primers_pos, $temp_left_pos;
        push @right_primers, $temp_right_seq;
        push @right_primers_pos, $temp_right_pos;
        
        #Produce the amplimers and place in an array
        my $right = scalar reverse $temp_right_seq;
        $right =~ tr/AGCT/TCGA/;
        $temp_out =~ /($temp_left_seq)(.+)($right)/;
        my $amplimer = $1.$2.$3;
        push @amplimers, $amplimer;
    } else {
        $temp_out =~ /(PRIMER_LEFT_SEQUENCE=)(.+)(PRIMER_RIGHT_SEQUENCE=)(.+)(PRIMER_LEFT=)(.+)(\,.+)(PRIMER_RIGHT=)(.+?)(\,.+)/;
        #Find primers and their positions, then push them into the correct arra
        my $temp_left_seq = $2;
        my $temp_left_pos = $6;
        my $temp_right_seq = $4;
        my $temp_right_pos = $9;
        push @left_primers, $temp_left_seq;
        push @left_primers_pos, $temp_left_pos;
        push @right_primers, $temp_right_seq;
        push @right_primers_pos, $temp_right_pos;
        
        #Produce the amplimers and place in an array
        my $right = scalar reverse $temp_right_seq;
        $right =~ tr/AGCT/TCGA/;
        $temp_out =~ /($temp_left_seq)(.+)($right)/;
        my $amplimer = $1.$2.$3;
        push @amplimers, $amplimer;
    }
}

#Display the primers, their position and the amplimers produced by them for the user
for (my $i = 0; $i < $results->number_of_results; $i++) {
    print "Primer set $i: \t $left_primers[$i] \t $right_primers[$i] \n";
    print "Position: \t $left_primers_pos[$i] \t\t\t $right_primers_pos[$i]\n";
    print "Amplimer $i: \t $amplimers[$i]\t\n";
}
print "\n";

#*******************************************************************************************************#
#PERFORM REMOTE BLAST TESTING ON AMPLIMERS (HUMAN FILTER) (Phuong Ma, sourced from Mark A. Jensen)
#WARNING: This segment takes MINUTES to run (~5 mins per BLAST search)
#Create a new list for standaloneblast to perform a search
my $blast_directory = './blast/blast/bin';
my $blast_database = 'nr';
my $blast_result = 1;
my @remote_blast_params = (-prog_dir => $blast_directory, #flag to specify directory where blast programs are
                           -db_name => $blast_database,
                           -remote => $blast_result );    #get your results
my $remote_blast = Bio::Tools::Run::StandAloneBlastPlus->new(@remote_blast_params);

#Take a fasta file to perform the BLAST search on and output the results to file query.bls
my $query_file = 'data/Primer3Test.fa';
my $output_file = 'query.bls';
my $method = [ '-num_alignments' => 10 ];
my @result_params = (-query => $query_file,        
                     -outfile => $output_file,     
                     -method_args => $method);
my $result = $remote_blast->blastn(@result_params);

$remote_blast->cleanup;

#*******************************************************************************************************#
#PERFORM LOCAL BLAST TESTING ON PRIMERS (PRIMER-DIMER FILTER)
#Bio::Tools::Run::StandAloneBlast?



#*******************************************************************************************************#
#CREATE DYNAMIC VIEWER
#GBrowse?



#*******************************************************************************************************#
#DISPLAY OUTPUT OF SUITABLE PRIMERS -> SEQUENCE, POSITION (DYNAMIC VIEWER?), TM
#DISPLAY OUTPUT OF THEIR AMPLIMERS -> SEQUENCE, BLAST RESULTS


