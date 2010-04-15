#!/usr/bin/perl
# combl.pl     Emanuel Heitlinger     2009/07/14 10:31:47

use warnings;      # avoid d'oh! bugs
use strict;        # avoid d'oh! bugs
$|=1;
use Data::Dumper;  # easy printing, sometimes only used during coding
use Getopt::Long;  # support options+arguments
use Pod::Usage;    # avoid redundant &Usage()
use File::Basename;
use Bio::SeqIO;
use Bio::SearchIO;
use Cwd;

# definition of variables needed for the help, man and debug output
my ($opt_help, $opt_man);
my $source_dir=cwd;

#own variables for command line args (optional)

my $infile;
my $outfile;
my @db;
my @program;
my $rmfalse;
my $evalue;

#save the call for documentation before we get the options
my $call = join (" ", @ARGV);

GetOptions(
	   'infile=s'      => \$infile,
	   'database=s'      => \@db,
	   'program=s' => \@program, 
           'evalue=s'  => \$evalue,
	   'outfile=s' => \$outfile,
	   'rmfalse=s' => \$rmfalse,
	   'help!'     => \$opt_help,                                # for all the standard 
	   'man!'      => \$opt_man,                                 # command-line  args
	  ) or pod2usage(-verbose => 1) && exit;                    # or exit printing usage
pod2usage(-verbose => 1) && exit if defined $opt_help;        # obligatory arguments wrongly
pod2usage(-verbose => 2) && exit if defined $opt_man;

# check existence of obligatory commmand-line arguments (optional)
pod2usage(-verbose => 1) && exit unless $infile;
pod2usage(-verbose => 1) && exit unless @db;
pod2usage(-verbose => 1) && exit if scalar(@db)!=scalar(@program);

#default for evalue
$evalue = 1e-3 unless $evalue;

my @db_base;
my $inbase = basename $infile;

for (@db) {
  push (@db_base, basename $_);
}

$outfile ="$inbase"."_vs_".join("and",@db_base).".csv" unless $outfile;

open (LOG, ">>logfile.cobl") or die "can't open logfile: $!";
select LOG;

if (-e "$outfile") {
  system ("rm $outfile");
  print "\n DELTED old outputfile $outfile\n";
}

# now we have all input, code starts here

my @all = &read_infile();
my %numreads_of = %{shift(@all)};
my %gc_of= %{shift(@all)};
my %seq_of = %{shift(@all)};
my %length_of = %{shift(@all)};


my %db_of;
my $pid = undef;
for (my $i=0; $i<$#db+1; $i++){
  if (defined($pid = fork)){
    if ($pid) {
      $db_of{$pid} =$db_base[$i];
    }
    else {
      &blast("$program[$i]", "$db[$i]", "$db_base[$i]");
      exit;
    }
  }
}
for (my $i=0 ;$i<$#db+1; $i++){
  wait;
}

my %evalue_of;
my %bscore_of;
my %hit_desc_of;
for my $db (values  %db_of) {
  my $br = new Bio::SearchIO (-file => ".temp_"."$db" , -format=>'blast'); 
  while (my $result = $br->next_result){
    my $query_acc = $result->query_accession;
    while (my $hit = $result->next_hit){ 			
      $hit_desc_of{$query_acc}{$db}=  $hit->description;
      $evalue_of{$query_acc}{$db}=  $hit->expect; 
      $bscore_of{$query_acc}{$db}=  $hit->bits; 
    }
  }
}

my $out_inf = join("", @db_base);
open (OUTT, ">$outfile") or die "cant open outfile:$!";
select OUTT;
my $line = "seq_name\tdbhit\tevalue\tbitscore\tannot\tnumreads\tgc\tlength\tsequence\n";
for my $seq_name(@all) {
  print "$line";
  my $max_score = 0;                                         
  if (grep {defined($_)} values %{$evalue_of{$seq_name}}){  # run this block if any database has a evalue
    for my $db_name(@db_base) {
      my $act_score = $bscore_of{$seq_name}{$db_name};          #get the actual e value of the different dbs
      if (defined $evalue_of{$seq_name}{$db_name} and $max_score <= $act_score) {	
	$line = "$seq_name\t$db_name\t$evalue_of{$seq_name}{$db_name}\t$bscore_of{$seq_name}{$db_name}\t$hit_desc_of{$seq_name}{$db_name}\t$numreads_of{$seq_name}\t$gc_of{$seq_name}\t$length_of{$seq_name}\t$seq_of{$seq_name}\n"; 
	$max_score = $act_score;                                     
      }
    }
  } 
  else {
    $line= "$seq_name\tNo_hit\tnone\tnone\tnone\t$numreads_of{$seq_name}\t$gc_of{$seq_name}\t$length_of{$seq_name}\t$seq_of{$seq_name}\n"; 
  }
}
close(OUTT);
select LOG;

print "\n output writte to $outfile\n\n";

unless ($rmfalse) { 
  for (values %db_of) {
    system ("rm .temp_"."$_");
  } 

}

###############################subroutines################################################################
sub read_infile{
  my @all;
  my %numreads_of;
  my %gc_of;
  my %length_of;
  my %sequence_of;
  my $ins_count = 0;
  my $in  = Bio::SeqIO->new(-file => "$infile" , '-format' => 'Fasta') or die "can't open input file";
  my $ins_s_count =0;
  while ( my $seq = $in->next_seq() ) {
    push(@all, $seq->id);
    my $info=$seq->desc;
    if ($info =~ m/length=(\d+).*numreads=(?<numreads>\d+).*/){
      $numreads_of{$seq->id}= $+{numreads};
    }
    else {
      $numreads_of{$seq->id}= 1;
    }      
    $gc_of{$seq->id}= &gc_percent($seq->seq); 
    $length_of{$seq->id}= length($seq->seq); 
    $sequence_of{$seq->id}= $seq->seq;  
    ++$ins_count;
  }
  print " read in $ins_count sequences from $infile\n"; 
  return(\%numreads_of, \%gc_of, \%sequence_of, \%length_of,  @all);
}
######
sub blast{
  my @what=@_;
  print "\n\tblasting against $what[2]\n";
  my $cmd = "blastall -p $what[0] -d $what[1] -i $infile -e $evalue -a 2 -b 1 -v 1 -o .temp_"."$what[2]";
  system($cmd); 
  print "\n\tblasting against $what[2] finished";
}
#####
sub gc_percent{
  my $sequence=$_[0];
  my $gc_count;
  my @DNA=split('',$sequence);
  for my $base(@DNA){
    if ($base =~ /^(G|C)/i){
      ++$gc_count;
    }
  }
  my $gc_percent = ($gc_count/length($sequence))*100;
  $gc_percent = sprintf("%.2f", $gc_percent);
  return($gc_percent);
}



__END__

=head1 NAME

blacomp.pl

=head1 SYNOPSIS

cobl.pl -i [input-file] -db [full_path/db1] -db [full_path/db2] -p [blastx] -p [blastn] -o [outfile]

=head1 DESCRIPTION

Stub documentation for blacomp.pl, 
created by template.el.

It looks like the author of this script was negligent 
enough to leave the stub unedited.

=head1 ARGUMENTS

 --i         input fasta-file
 --db        database(s) to blast against
 --outfile   output file (will be overwritten if exists) DFAULT: concat infile and databases
 --program   programm e.g. blastx or blastn to use for blasting
             db and programm must have same length!!!
 --help      print Options and Arguments 
 --man       print complete man page

=head1 AUTHOR

Emanuel Heitlinger, emanuelheitlinger@gmail.com

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2009 by Emanuel Heitlinger

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
