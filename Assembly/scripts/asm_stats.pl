#!/usr/bin/env perl

use File::Spec;
use strict;
use warnings;
use Getopt::Long;

my %stats;
my $debug = 0;

GetOptions('v|verbose' => \$debug );

my $BUSCO_dir = 'BUSCO';
my $BUSCO_pep = 'BUSCO_pep';
my $read_map_stat = 'mapping_report';
my $SAMPLES = 'samples.csv';
my $dir = shift || 'genomes';
my @header;
my %header_seen;
my (@mappingrun,@busco_rerun);
opendir(DIR,$dir) || die $!;
my $first = 1;
foreach my $file ( readdir(DIR) ) {
    next unless ( $file =~ /(\S+)(\.fasta)?\.stats.txt$/);
    my $stem = $1;
    $stem =~ s/\.sorted//;
	my $stemclean = $stem;
	$stemclean =~ s/\.AAFTF$//;
    #warn("$file ($dir)\n");
    open(my $fh => "$dir/$file") || die "cannot open $dir/$file: $!";
    while(<$fh>) {
	next if /^\s+$/;
	s/^\s+//;
	chomp;
	if ( /\s*(.+)\s+=\s+(\d+(\.\d+)?)/ ) {
	    my ($name,$val) = ($1,$2);	    
	    $name =~ s/\s*$//;
	    $name =~ s/\s+/_/g;
	    $stats{$stem}->{$name} = $val;

	    if( ! $header_seen{$name} ) {
		push @header, $name;
		$header_seen{$name} = 1;
	    }
	}
    }
    if ( -d $BUSCO_dir ) {
	if ( $first ) { 
	    push @header, qw(BUSCO_Complete BUSCO_Single BUSCO_Duplicate
			     BUSCO_Fragmented BUSCO_Missing BUSCO_NumGenes
		);
	}
	my $buscosub = File::Spec->catdir($BUSCO_dir,$stem);
	if ( -d $buscosub ) {
	    opendir(BUSCOD, $buscosub);
	    my @busco_files;
	    foreach my $file ( readdir(BUSCOD) ) {
			if ($file =~ /short_summary.specific.([^\.]+)\.\S+\.txt/) {
				push @busco_files, File::Spec->catfile($buscosub,$file);
			}
	    }	
	    if ( @busco_files ) {
		my $busco_file = shift @busco_files; # not sure what to do if there are multiple in same dir, don't think this will happen

		open(my $fh => $busco_file) || die $!;		
		while(<$fh>) {	 
		    if (/^\s+C:(\d+\.\d+)\%\[S:(\d+\.\d+)%,D:(\d+\.\d+)%\],F:(\d+\.\d+)%,M:(\d+\.\d+)%,n:(\d+)/ ) {
			$stats{$stem}->{"BUSCO_Complete"} = $1;
			$stats{$stem}->{"BUSCO_Single"} = $2;
			$stats{$stem}->{"BUSCO_Duplicate"} = $3;
			$stats{$stem}->{"BUSCO_Fragmented"} = $4;
			$stats{$stem}->{"BUSCO_Missing"} = $5;
			$stats{$stem}->{"BUSCO_NumGenes"} = $6;
		    } 
		}
	    } else {
			warn("Cannot find BUSCO result in $buscosub\n");
			
			my $r =	`tail -n +2 $SAMPLES | grep -n $stemclean`;
			my ($id) = ($r =~ /^(\d+):.+/);
			push @busco_rerun, $id;
			#warn($id);
	    }
	} else {
	    warn("BUSCO not run yet on $buscosub\n");
	    my $r = `tail -n +2 $SAMPLES | grep -n $stemclean`;
	    my ($id) = ($r =~ /^(\d+):.+/);
	    push @busco_rerun, $id;
	    #warn($id);
	}
    }

    if ( -d $BUSCO_pep ) {
	if ( $first ) { 
	    push @header, qw(BUSCOP_Complete BUSCOP_Single BUSCOP_Duplicate
			     BUSCOP_Fragmented BUSCOP_Missing BUSCOP_NumGenes
		);
	}
	my $buscosub = File::Spec->catdir($BUSCO_pep,$stem);

	if ( -d $buscosub ) {
	    opendir(BUSCOD, $buscosub);
	    my @busco_files;
	    foreach my $file ( readdir(BUSCOD) ) {
			if ($file =~ /short_summary.specific.([^\.]+)\.\S+\.txt/) {
				push @busco_files, File::Spec->catfile($buscosub,$file);
			}
	    }	
	    if ( @busco_files ) {
		my $busco_file = shift @busco_files; # not sure what to do if there are multiple in same dir, don't think this will happen

		open(my $fh => $busco_file) || die $!;				
		while(<$fh>) {
		    if (/^\s+C:(\d+\.\d+)\%\[S:(\d+\.\d+)%,D:(\d+\.\d+)%\],F:(\d+\.\d+)%,M:(\d+\.\d+)%,n:(\d+)/ ) {
			$stats{$stem}->{"BUSCOP_Complete"} = $1;
			$stats{$stem}->{"BUSCOP_Single"} = $2;
			$stats{$stem}->{"BUSCOP_Duplicate"} = $3;
			$stats{$stem}->{"BUSCOP_Fragmented"} = $4;
			$stats{$stem}->{"BUSCOP_Missing"} = $5;
			$stats{$stem}->{"BUSCOP_NumGenes"} = $6;
		    } 
		}

	    } else {
		warn("Cannot find BUSCO result files in $buscosub\n");
	    }
	}
    }

    if ( -d $read_map_stat ) {
		my $sumstatfile = File::Spec->catfile($read_map_stat,
					      sprintf("%s.bbmap_summary.txt",$stem));
		warn("sumstat is $sumstatfile\n") if $debug;
		if ( -f $sumstatfile ) {
			open(my $fh => $sumstatfile) || die "Cannot open $sumstatfile: $!";
			my $read_dir = 0;
			my $base_count = 0;
			$stats{$stem}->{'Mapped_reads'} = 0;
	    	while(<$fh>) {
				if( /Read\s+(\d+)\s+data:/) {
					$read_dir = $1;
				} elsif( $read_dir && /^mapped:\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+)/) {
					$base_count += $4;
					$stats{$stem}->{'Mapped_reads'} += $2;
				}  elsif( /^Reads Used:\s+(\S+)/) {
					$stats{$stem}->{'Total_Reads'} = $1;
				}
	   		}
			$stats{$stem}->{'Average_Coverage'} = 0;
			if (exists $stats{$stem}->{'TOTAL_LENGTH'} ) {
				$stats{$stem}->{'Average_Coverage'}  = sprintf("%.1f",$base_count / $stats{$stem}->{'TOTAL_LENGTH'});
			}
			if( $first )  {
				push @header, ('Total_Reads',
						'Mapped_reads',			   
						'Average_Coverage');
			}
		} else {
			warn("cannot find $sumstatfile\n");
			my $r = `tail -n +2 $SAMPLES | grep -n $stemclean`;
	    	my ($id) = ($r =~ /^(\d+):.+/);
			push @mappingrun, $id;
		}
    }
    $first = 0;
}
print join("\t", qw(SampleID), @header), "\n";
foreach my $sp ( sort keys %stats ) {
    print join("\t", $sp, map { $stats{$sp}->{$_} || 'NA' } @header), "\n";
}
if ( scalar(@busco_rerun) > 0 ) {
	warn("BUSCO rerun: -a",join(",",sort { $a <=> $b } @busco_rerun),"\n");
}
if ( scalar(@mappingrun) > 0 ) {
	warn("mapping rerun: -a",join(",",sort { $a <=> $b } @mappingrun),"\n");
}