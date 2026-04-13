#!/usr/bin/perl

use strict;
use File::Basename;

my $depthSAM = $ARGV[0];
my $tsvNF = $ARGV[1];
my $finalFile = $ARGV[2];
my $consensus = $ARGV[3];
my $bamStats = $ARGV[4];


sub updateReport {
  my @inputs = @_;
  my $id = $inputs[0];
  $id=~s/\.txt|//g;
  $id=~s/.+\.//g;
  $id=~s/(_[^_]+)_/$1./g;

  my $line = `grep -w $id $inputs[1]`;
  chomp($line);
  my @parameters = split(/\t/,$line);

  my $depth = `tail -n 1 $inputs[0]`;
  chomp($depth);
  my @samtoolDepth = split(/\t/,$depth);

  my $sequence = `grep -v \"\>\" $inputs[2]`;
  chomp($sequence);
  $sequence=~s/\n//g;
  my @bases = split(//,$sequence);
  my ($A,$T,$C,$G,$N)=0;
  foreach my $base (@bases){
    if(uc($base) eq "A"){
        $A++;;
    }elsif(uc($base) eq "T"){
        $T++;
    }elsif(uc($base) eq "C"){
        $C++;
    }elsif(uc($base) eq "G"){
        $G++;
    }elsif(uc($base) eq "N"){
        $N++;
    }
  }
  my $knowBases=$A+$T+$C+$G;

  # Skip this virus if no BAM-mapped reads or consensus has no callable bases
  # (no position passed the quality/depth/agreement thresholds)
  if ($samtoolDepth[3] == 0 || $knowBases == 0) {
    return "";
  }

  # Read BAM-derived read identity and mean length
  my $bamstats_line = `cat $inputs[3]`;
  chomp($bamstats_line);
  my @bamstats = split(/\t/, $bamstats_line);
  my $mean_identity = $bamstats[0];
  my $mean_length = $bamstats[1];

  # All metrics derived from BAM alignment:
  #   samtoolDepth[3] = numreads    (mapped reads in BAM)
  #   samtoolDepth[5] = coverage    (horizontal coverage %)
  #   samtoolDepth[6] = meandepth   (vertical coverage)
  #   samtoolDepth[7] = meanbaseq   (mean base quality)
  #   mean_identity   = mean read identity from NM tags
  #   mean_length     = mean mapped read length
  my $results = "$parameters[1]\t$parameters[2]\t$parameters[3]\t$parameters[4]\t$samtoolDepth[3]\t$parameters[6]\t$parameters[7]\t$samtoolDepth[5]\t$samtoolDepth[6]\t$samtoolDepth[5]\t0\t$mean_identity\t$mean_length\t$samtoolDepth[7]\n";

  return $results;
}

if(-e $finalFile){
  open(FIN,">>$finalFile");
  my $update = updateReport($depthSAM,$tsvNF,$consensus,$bamStats);
  print FIN $update;
  close FIN;

}else {
  open(FIN,">$finalFile");
  open(TSV,"$tsvNF");
  my $line;
  while($line=<TSV>){
    chomp($line);
    if($line=~m/VirusName/){
      my @temp=split(/\t/,$line);
      print FIN "$temp[0]\t$temp[1]\t$temp[2]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t\"Coverage\"\t\"DepthAverage\"\t\"HorizontalCoverage\"\t\"N_content\"\t$temp[9]\t$temp[10]\tMeanBaseQuality\n";
    }
  }
  close TSV;
  my $update = updateReport($depthSAM,$tsvNF,$consensus,$bamStats);
  print FIN $update;
  close FIN;
}
