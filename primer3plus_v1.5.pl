#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Bio::Tools::Run::Primer3;
use Bio::Tools::Run::StandAloneBlastPlus;

#AUTHOR: Team Dynamic
#DATE STARTED: Mar 3/15
#DATE UPDATED: Mar 6/15
#PURPOSE: Produce primers from an inputted DNA sequence (fasta file) using primer3 tool
#NOTES: Must install BioPerl first following the instructions at http://www.bioperl.org/wiki/Installing_BioPerl_on_Windows
        #Must download Bio/Tools/Run and Bio/Roots/Roots modules and place in Perl>lib for primer3 tool to be accessed
        #Must copy primer3_core.exe in same location as THIS perl file to avoid "primer3 cannot be found" error
        #Use primer3 version 1.1.4 to avoid "missing SEQUENCE tag" error
        #Alter -path => in Primer3->new() to avoid "SH: Command not found" error
        #Copy primer3_config directory in same location as THIS perl file to avoid "thermodynamic approach" error
        #To avoid "can't locate object method 'new'" error, install correct version of Primer3.pm

#*******************************************************************************************************#
#ATTACH JAVA GUI TO THIS PROGRAM FILE (SINGLEPLEX AND MULTIPLEX)



#*******************************************************************************************************#
#USER-SELECTED ORGANISM FROM GUI USED AS INPUT (SINGLEPLEX)



#*******************************************************************************************************#
#ACCESS ORGANISM'S FASTA FILE FROM DATABASE (SINGLEPLEX)



#*******************************************************************************************************#
#ALLOW USER TO UPLOAD CONCATENATED FASTA FILE WITH MULTIPLE SEQUENCES(MULTIPLEX)



#*******************************************************************************************************#
#PRIMER3 TOOL TO PRODUCE PRIMERS FROM FASTA FILE (SINGLEPLEX) (Rebecca Allan, sourced from Chad Matsalla)
#Input a fasta file and create a new primer3 list
my $fasta_file = 'data/Primer3Test.fa';
my $primer3_output_file = 'data/temp.out';
my $primer3_path = 'C:/Users/Rebecca/Desktop/IGP/PERL_PROGRAM/primer3_core'; #USE YOUR OWN PATH!

#Put the fasta file in the correct SeqIO format and declare the variables for primer3
my $sequence_io = Bio::SeqIO->new(-file=>$fasta_file);
my $sequence = $sequence_io->next_seq;
my $primer3 = Bio::Tools::Run::Primer3->new(-seq => $sequence,
                                            -outfile => $primer3_output_file,
                                            -path => $primer3_path); 

#Test to see if primer3_core.exe is within the directory
unless ($primer3->executable) {
    print STDERR "Primer3 can not be found. Is it installed?\n";
    exit(-1)
}

#Message for the user to let them know the program is running
print "Please wait while primer candidates are being found in the sequence. This may take up to 5 minutes. \n\n";

#Adjust default values of specific arguments in primer3
$primer3->add_targets('PRIMER_MIN_TM'=>56, 'PRIMER_MAX_TM'=>65, 'PRIMER_OPT_TM'=>60);
$primer3->add_targets('PRIMER_MIN_SIZE'=>20, 'PRIMER_MAX_SIZE'=>27, 'PRIMER_OPT_SIZE'=>20, 'PRIMER_DEFAULT_SIZE'=>20);
$primer3->add_targets('PRIMER_MIN_GC'=>30, 'PRIMER_MAX_GC'=>80, 'PRIMER_OPT_GC_PERCENT'=>50);
#Next line on primer product size range does not work, must be fixed
#$primer3->add_targets('PRIMER_PRODUCT_SIZE_RANGE'=>100..500);
#Number of primer candidates returned
$primer3->add_targets('PRIMER_NUM_RETURN'=>5);

#Run the actual primer3 application
my $primer3_results = $primer3->run;

#Print the number of results for the user to prove a successful run
print "There were ", $primer3_results->number_of_results, " primer candidates found. \n";

#*******************************************************************************************************#
#RETRIEVE PRIMERS FROM TEMP.OUT FILE (Rebecca Allan and Phuong Ma)
#Array to contain temp.out contents
my @primer3_file_contents;
#Arrays to contain the primers and their positions
my @left_primers = ();
my @left_primers_pos = ();
my @right_primers = ();
my @right_primers_pos = ();
#Array to contain amplimers for each primer pair
my @amplimers = ();

#Read the entire temp.out file into an array and then close file
open (my $primer3_fh, "<", $primer3_output_file) or die "Error reading $primer3_output_file file.\n";
foreach (<$primer3_fh>) {
    push @primer3_file_contents, $_;
}
close $primer3_fh;
chomp @primer3_file_contents;

#Scalar value $temp_out for array @file_contents for regex
my $primer3_temp_out = join '', @primer3_file_contents;

#For each primer candidate, collect the results
for (my $i = 0; $i < $primer3_results->number_of_results; $i++) {
    #Perform different regex statements depending on the index value of each primer pair element
    if ($i > 0) {
        $primer3_temp_out =~ /(PRIMER_LEFT_$i\_SEQUENCE=)(.+)(PRIMER_RIGHT_$i\_SEQUENCE=)(.+)(PRIMER_LEFT_$i=)(.+)(\,.+)(PRIMER_RIGHT_$i=)(.+?)(\,.+)/;
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
        $primer3_temp_out =~ /($temp_left_seq)(.+)($right)/;
        my $amplimer = $1.$2.$3;
        push @amplimers, $amplimer;
    } else {
        $primer3_temp_out =~ /(PRIMER_LEFT_SEQUENCE=)(.+)(PRIMER_RIGHT_SEQUENCE=)(.+)(PRIMER_LEFT=)(.+)(\,.+)(PRIMER_RIGHT=)(.+?)(\,.+)/;
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
        $primer3_temp_out =~ /($temp_left_seq)(.+)($right)/;
        my $amplimer = $1.$2.$3;
        push @amplimers, $amplimer;
    }
}

#*******************************************************************************************************#
#QC TO REMOVE PRIMERS THAT ARE NOT IDEAL CANDIDATES/REDUNDANT (Rebecca Allan)
#Delete primers that have the same left/left or right/right primer positions
for(my $a = 0; $a < @left_primers; $a++) {
    for(my $b = 0; $b < @left_primers; $b++) {
        #If elements do not have the same index value (not the primer being compared to itself)
        if ($a != $b) {
            #If the position numbers are the same (thus the primer will be the same)
            if ($left_primers_pos[$a] == $left_primers_pos[$b]) {
                #Remove the primer pair and its amplimer from the arrays
                splice @left_primers, $b, 1;
                splice @left_primers_pos, $b, 1;
                splice @right_primers, $b, 1;
                splice @right_primers_pos, $b, 1;
                splice @amplimers, $b, 1;
            #Repeat comparison for right primers if the left primer positions are not equal
            } elsif ($right_primers_pos[$a] == $right_primers_pos[$b]) {
                splice @left_primers, $b, 1;
                splice @left_primers_pos, $b, 1;
                splice @right_primers, $b, 1;
                splice @right_primers_pos, $b, 1;
                splice @amplimers, $b, 1;
            }
        }
    }
}

print "After quality control ", (($primer3_results->number_of_results) - @left_primers), " primer candidates were eliminated. \n\n";

#*******************************************************************************************************#
#VIEW PRIMER/AMPLIMER RESULTS (Rebecca Allan)
#Display the primers, their position and the amplimers produced by them for the user
for (my $i = 0; $i < @left_primers; $i++) {
    #Use printf for formatting purposes
    printf "%10s\t%25s\t%25s", "Primer Set $i:", "$left_primers[$i]", "$right_primers[$i]";
    print "\n";
    printf "%10s\t%25s\t%25s", "Position $i:", "$left_primers_pos[$i]", "$right_primers_pos[$i]";
    print "\n";
    #printf "%10s\t%25s", "Amplimer $i:", "$amplimers[$i]";
    #print "\n";
}
print "\n";

#*******************************************************************************************************#
#CREATE CONCATENATED FASTA FILES THAT CONTAIN PRIMERS AND AMPLIMERS (Phuong Ma and Rebecca Allan)
#Open filehandle, $primers_out, and create concatenated fasta file for all primers
my $primer_filename = 'data/primers.fa';
open (my $primers_out, '>', "$primer_filename") or die "Could not create $primer_filename.";

#Append every left primer to concatenated fasta file
for (my $i = 0; $i < @left_primers; $i++) {
    print $primers_out ">Left Primer $i\n";
    print $primers_out "$left_primers[$i]\n";
}
#Append every right primer to concatenated fasta file
for (my $j = 0; $j < @right_primers; $j++) {
    print $primers_out ">Right Primer $j\n";
    print $primers_out "$right_primers[$j]\n";
}

#Close filehandle
print "Primers written to $primer_filename\n";
close $primers_out;

#Open filehandle, $amplimers_out, and create concatenated fasta file for all amplimers
my $amplimer_filename = 'data/amplimers.fa';
open (my $amplimers_out, ">", "$amplimer_filename") or die "Could not create $amplimer_filename.";

#Append every amplimer to concatenated fasta file
for (my $k = 0; $k < @left_primers; $k++) {
    my $amplimers_info = ">Amplimer $k from $fasta_file \n$amplimers[$k]";
    $amplimers_info =~ s/(.{1,80})/$1\n/g;
    print $amplimers_out $amplimers_info;
}

#Close filehandle
print "Amplimers written to $amplimer_filename\n\n";
close $amplimers_out;

#*******************************************************************************************************#
#PERFORM REMOTE BLAST TESTING ON AMPLIMERS (HUMAN FILTER) (Phuong Ma and Rebecca Allan, sourced from Mark A. Jensen)
#WARNING: This segment takes MINUTES to run (~5 mins per BLAST search)
#Create sample fasta file for amplimer to be read by StandAloneBlast
open (my $amplimer_out, ">", "data/amplimer0.fa") or die "Could not create file.";
my $amplimer_info = ">Amplimer 0 from $fasta_file \n$amplimers[0]";
$amplimer_info =~ s/(.{1,80})/$1\n/g;
print $amplimer_out $amplimer_info;
close $amplimer_out;

my $user_result = 0;

my $blast_directory = './blast/blast/bin';
my $blast_database = 'nr';
my $blast_result = 1;

my $remote_query_file = "data/amplimer0.fa";
my $remote_blast_output_file = 'data/query.bls';
my $remote_method = [ '-num_alignments' => 10 ];

#Start a loop that requires the user to input the correct input ("y" or "n") to escape it
while ($user_result == 0) {
    print "Would you like to perform a remote BLAST search on the first amplimer? (y/n) ";
    my $user_prompt = <STDIN>;
    chomp $user_prompt;
    
    if ($user_prompt eq "y") {
        print "Please wait while your remote BLAST search is being performed. This may take up to 10 minutes.\n";
        
        #Create a new list for standaloneblast to perform a search
        my @remote_blast_params = (-prog_dir => $blast_directory, #flag to specify directory where blast programs are
                                   -db_name => $blast_database,
                                   -remote => $blast_result );    #get your results
    
        #Take a fasta file to perform the BLAST search on and output the results to file query.bls
        my @result_params = (-query => $remote_query_file,        
                             -outfile => $remote_blast_output_file,     
                             -method_args => $remote_method);

        #Run the remote blast search
        my $remote_blast = Bio::Tools::Run::StandAloneBlastPlus->new(@remote_blast_params);
        my $result = $remote_blast->blastn(@result_params);
        $remote_blast->cleanup;
        
        #Send a message to the user to know the remote blast search has been completed
        print "Your results are in $remote_blast_output_file. \n\n";
        $user_result = 1;
    } elsif ($user_prompt eq "n") {
        $user_result = 1;
    #Continue the loop if the user enters an invalid input
    } else {
        print "Please enter a valid answer!\n";
        $user_result = 0;
    }
}

print "\n";

#*******************************************************************************************************#
#PERFORM LOCAL BLAST TESTING ON PRIMERS (PRIMER-DIMER FILTER) (Phuong Ma, sourced from Mark A. Jensen)
$user_result = 0;

my $local_database_name = 'primers';
my $local_database_data = 'data/Primer3Test.fa';
my $query_file = "data/primers.fa";
my $local_blast_output_file = 'data/primer_local_blast.txt';
my $local_method = [ '-num_alignments' => 10 ];

while ($user_result == 0) {
    print "Would you like to perform a local BLAST search on the primers? (y/n) ";
    my $user_prompt = <STDIN>;
    chomp $user_prompt;
    
    if ($user_prompt eq "y") {
        print "Please wait while your local BLAST search is being performed. This may take up to 10 minutes.\n";
        
        my @local_blast_params = (-prog_dir => $blast_directory, #flag to specify directory where blast programs are
                                  -db_name => $local_database_name,
                                  -db_data => $local_database_data,
                                  -create => 1);

        my @local_result_params = (-query => $query_file,
                                   -outfile => $local_blast_output_file,
                                   -method_args => $local_method
                                   -task => "blastn-short" );
    
        my $local_blast = Bio::Tools::Run::StandAloneBlastPlus->new(@local_blast_params);
        my $local_result = $local_blast->blastn(@local_result_params);
        $local_blast->cleanup;
        
        print "Your results are in $local_blast_output_file. \n\n";
        $user_result = 1;
    } elsif ($user_prompt eq "n") {
        $user_result = 1;
    } else {
        print "Please enter a valid answer!\n";
        $user_result = 0;
    }
}

print "\n";
print "Thank you for using the DyNAmic Viewer! Your results are in the data directory. \n";

#*******************************************************************************************************#
#CREATE DYNAMIC VIEWER
#GBrowse?



#*******************************************************************************************************#
#DISPLAY OUTPUT OF SUITABLE PRIMERS -> SEQUENCE, POSITION (DYNAMIC VIEWER?), TM
#DISPLAY OUTPUT OF THEIR AMPLIMERS -> SEQUENCE, BLAST RESULTS


