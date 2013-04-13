#!/usr/bin/perl -w
use strict;
use Getopt::Long;

#----------------globals------------------------------------------------------------------
our ( $SV, $FASTA, $COVERAGE );
GetOptions(
	'sv=s'       => \$SV,
	'sequence=s' => \$FASTA,
	'coverage=s' => \$COVERAGE
);

#@chr is 1-indexing,$chr[index] ="0NNNNNNNN..."  N=A/T/C/G
my @chr;
my $chr_num = 0;

my @chr_link;

# $chr_link[current]->{'next'} =
#  -1:  printed
#		null: index + 1 or -1 according to the sign from SV
#		other: jumping to other potion

my @sv = ();

# $sv[index]->{'chromosome'} = 1 or 2
# $sv[index]->{'first'} = the first breakpoint
# $sv[index]->{'second'} =  the second breakpoint
my @inverse_sv;

#$inverse_sv[chr_index]{sv}= index of sv
my @sorted;

#array of SVs, sorted in ascending order and chromasome index order

my @coverage;

# $coverage[chr#]->{pos}=coverage

#--------------------main------------------------------------------------------------------------------------------------------------------------------------------
$| = 1;
read_fasta();

read_sv();

foreach ( 0 .. $#sv ) {
	$inverse_sv[ $sv[$_]->{'chromosome'} ]{ abs( $sv[$_]->{'first'} ) }  = $_;
	$inverse_sv[ $sv[$_]->{'chromosome'} ]{ abs( $sv[$_]->{'second'} ) } = $_;
}

foreach ( 1 .. $chr_num ) {
	push( @sorted, -1 );    #delimiter to separate the chromosome
	my @tmp = sort { $a <=> $b } ( keys %{ $inverse_sv[$_] } );
	push( @sorted, @tmp );
}

#print join "\n", @sorted;
#print "\n\n\n";

print "Reading coverage file";
foreach ( 1 .. $chr_num ) {
	read_coverage($_);
}
print "\nDone\n";

locate_breakpoint();

foreach ( 0 .. $#sv ) {
	print $sv[$_]->{'chromosome'} . "\t"
	  . $sv[$_]->{'first'} . "\t"
	  . $sv[$_]->{'second'} . "\t"
	  . $sv[$_]->{'first_bp'} . "\t"
	  . $sv[$_]->{'second_bp'} . "\n";
}

foreach ( 1 .. $chr_num ) {
	@chr_link = ();

	print "Linking breakpoints for chromosome $_...";
	make_link($_);
	print "Done\n";

	print "Generating output for chomosome $_...";
	print_chr($_);
	print "Done\n";
}

#----------------Read the files-----------------------------------------------------------
sub read_fasta {
	print "Reading FASTA file";
	open( CHR, "$FASTA" )
	  or die("cannot found the fasta file!");
	my $count = 0;
	while (<CHR>) {
		$count++;
		my $line = trim_line($_);
		if ( $line =~ />/ ) {

			#swith chromatin
			$chr_num++;
			$chr[$chr_num] = 0;
			next;
		}
		$chr[$chr_num] = $chr[$chr_num] . $line;
		print "." if ( $count % 5000 == 0 );
	}
	close(CHR);
	print "\nDone\n";
}

sub read_sv {
	print "Reading SV file...";
	open( SV, "$SV" ) or die("cannot find the sv file!");
	while (<SV>) {
		chomp;
		s/\s//;
		next if ( /^#/ || /^\s*$/ );    #annotation, skip
		if (m/([\+\-])chr(\d+):(\d+)\s*([\+\-])chr\d+:(\d+)/) {
			my ( $sign1, $chr, $pos1, $sign2, $pos2 ) = ( $1, $2, $3, $4, $5 );
			$sign1 = ( ( $sign1 eq "+" ) ? 1 : -1 );
			$sign2 = ( ( $sign2 eq "+" ) ? 1 : -1 );
			push(
				@sv,
				{
					'first'      => $sign1 * $pos1,
					'second'     => $sign2 * $pos2,
					'chromosome' => $chr
				}
			);
		}
	}
	close(SV);
	print "\nDone\n";
}

sub read_coverage {
	open( COV, "$COVERAGE" )
	  or die("Cannot fount the coverage file!");
	my $chr_index = shift;
	my $line;
	my $temp;
	my $chr;
	my $start;
	my $end;
	my $intensity;
	my $count=0;
	while (<COV>) {
		$line = $_;
		( $chr, $temp, $temp, $start, $end, $intensity ) =
		  split( "\t", $line );
		foreach ( $start .. $end ) {
			if ( $chr eq 'chr' . $chr_index ) {
				$coverage[$chr_index]->{$_} = $intensity;
			}
		}
		$count++;
		print "." if($count%100000==0);
	}
	close(COV);
}

#-------------------functions------------------------------------------------------
sub print_chr {
	my $chr_index = shift @_;
	open( OUT, ">>out.fasta" );
	if ( $chr_index != 1 ) {
		print OUT "\n>chr$chr_index\n";
	}
	else {
		print OUT ">chr$chr_index\n";
	}
	my $direction = 1;
	my $cur       = 1;
	while ( $cur != length( $chr[$chr_index] ) ) {
		if ( $direction == 1 ) {
			print OUT substr( $chr[$chr_index], $cur, 1 );
		}
		else {
			##print the reverse complement
			my $base = substr( $chr[$chr_index], $cur, 1 );
			print OUT "T" if ( $base eq "A" or $base eq "a" );
			print OUT "A" if ( $base eq "T" or $base eq "t" );
			print OUT "G" if ( $base eq "C" or $base eq "c" );
			print OUT "C" if ( $base eq "G" or $base eq "g" );
		}
		if ( $chr_link[$cur]->{'next'} && $chr_link[$cur]->{'next'} != -1 ) {
			$direction = -$direction
			  if ( $chr_link[$cur]->{'switch'}
				&& $chr_link[$cur]->{'switch'} == 1 );

			#jump
			my $next = $chr_link[$cur]->{'next'};
			$chr_link[$cur]->{'next'}  = -1;      #mark as printed
			$chr_link[$next]->{'next'} = -1;      #mark as printed
			$cur                       = $next;
		}
		else {

			#one by one
			$chr_link[$cur]->{'next'} = -1;       #mark as printed
			while ($chr_link[$cur]->{'next'}
				&& $chr_link[$cur]->{'next'} == -1 )
			{
				$cur += $direction;
			}
		}
	}

	close OUT;
}

sub make_link {
	my $chr_index = shift @_;
	foreach my $sv_tmp (@sv) {
		if ( $sv_tmp->{'chromosome'} == $chr_index ) {
			if ( $sv_tmp->{'first'} * $sv_tmp->{'second'} < 0 ) {

				#switch sign
				$chr_link[ abs( $sv_tmp->{'first_bp'} ) ]->{'switch'}  = 1;
				$chr_link[ abs( $sv_tmp->{'second_bp'} ) ]->{'switch'} = 1;
			}
			$chr_link[ abs( $sv_tmp->{'first_bp'} ) ]->{'next'} =
			  abs( $sv_tmp->{'second_bp'} );
			$chr_link[ abs( $sv_tmp->{'second_bp'} ) ]->{'next'} =
			  abs( $sv_tmp->{'first_bp'} );
		}
	}
}

sub trim_line {
	my $line = shift;
	chomp $line;
	$line =~ s/\s+$//;
	$line =~ s/^\s+//;
	$line;
}

sub locate_breakpoint {

  #take in all the svs from 1 chromosome
  #sort them in ascending order
  #the real break point is within a window of (b-500, b+500)
  #if |b1-b2|<1000, the two region overlap, they have almost the same breakpoint
  #locate the brak point b in this 1000 bases, and b1=b, b2=b+1
	print "Locating exact breakpoints...";
	my $chr_index = 0;
	foreach my $bp (@sorted) {

		#switch chromosome
		if ( $bp == -1 ) {
			$chr_index++;
			next;
		}

		#find the minimum as real breakpoint
		my $real_bp = $bp - 500;
		foreach my $i ( $bp - 499 .. $bp + 500 ) {
			$real_bp = $i
			  if ( $coverage[$chr_index]->{$i} <
				$coverage[$chr_index]->{$real_bp} );
		}
		my $tmp_index = $inverse_sv[$chr_index]{$bp};
		if ( abs( $sv[$tmp_index]->{'first'} ) eq $bp ) {
			$sv[$tmp_index]->{'first_bp'} = $real_bp;
		}
		else {
			$sv[$tmp_index]->{'second_bp'} = $real_bp;
		}

	}

	#make sure no overlap
	$chr_index = 0;
	foreach my $i ( 0 .. $#sorted ) {
		if ( $sorted[$i] == -1 ) {
			$chr_index++;
			next;
		}
		next if ( $sorted[ $i - 1 ] == -1 );
		my $tmp_index     = $inverse_sv[$chr_index]{ $sorted[$i] };
		my $tmp_index_pre = $inverse_sv[$chr_index]{ $sorted[ $i - 1 ] };

		if ( $sv[$tmp_index]->{'first_bp'} == $sv[$tmp_index_pre]->{'first_bp'}
			|| $sv[$tmp_index]->{'second_bp'} ==
			$sv[$tmp_index_pre]->{'first_bp'} )
		{
			$sv[$tmp_index_pre]->{'first_bp'}--;
		}
		elsif (
			$sv[$tmp_index]->{'first_bp'} == $sv[$tmp_index_pre]->{'second_bp'}
			|| $sv[$tmp_index]->{'second_bp'} ==
			$sv[$tmp_index_pre]->{'second_bp'} )
		{
			$sv[$tmp_index_pre]->{'second_bp'}--;
		}
	}
	print "done\n";
}

