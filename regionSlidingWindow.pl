#!/usr/bin/perl -w 

##takes in normalised wig file performs sliding window analysis over regions listed in ROI file.
##treats each ROI individually so that no windows are merged into one region
##avoids file read per region by seperating roi into distinct non-overlapping lists

use File::Sort qw(sort_file);
use Getopt::Long;
use strict;

#get user options
my %opt = (
	out => "windows.tab", 
	window => 50,
	step => 20,
	chrcol => 0,
	startcol => 1,
	endcol => 2,
	hitcol => 3,
	adjust => 1,
	);
            
GetOptions(\%opt, "window=i", "step=i", "norm=s", "regions=s", "out=s" );

print STDOUT ("Window = $opt{window}\nStep = $opt{step}\n"); #print options to info

##validate options
if ($opt{window}<1||$opt{step}<1){
	die "use of invalid option values";
}


open(ROI, "$opt{regions}" ) or die "cannot open $opt{regions}";

my %chromosomes; #represents the chromosomes covered in input files

#
# Read peaks file and get chromosome numbers
#

my $sortc;
my $sorts;

while(<ROI>){
	chomp;
	next if /#/;
	my($c,$s) = (split /\t/)[$opt{chrcol},$opt{startcol}]; #get chr numbers
	if (!exists $chromosomes{$c}){
		$chromosomes{$c} = 0;
		$sortc=$c;
		$sorts=$s;
	}
	elsif($c ne $sortc || $s<$sorts){
		print STDERR "Interval file is not sorted by chromosome and start position\n";
		exit;
	}
}
close ROI;

open(SORT, "$opt{norm}" ) or die "cannot open $opt{norm}";

my %chrsort; #represents the chromosomes covered in input files


#
# Read wig file and check sorting
#

$sortc="";
$sorts=0;

while(<SORT>){
	chomp;
	next if /track/;
	my($c,$s) = (split /\s/)[$opt{chrcol},$opt{startcol}]; #get chr numbers
	if (!exists $chrsort{$c}){
		$chrsort{$c} = 0;
		$sortc=$c;
		$sorts=$s;
	}
	elsif($c ne $sortc || $s<$sorts){
		print STDERR "Wiggle file is not sorted by chromosome and start position\n";
		exit;
	}
}
close SORT;
%chrsort=();

#variables for each chromo
my $score = 0; #output score for window, total number of reads or avg?
my @starts = (); #arrays containing start and end positions of each region
my @ends = ();
my @windows = (); #number of windows per region
my %norm = (); #represents normalised values at each base pair covered in a specific region e.g. 1 chromosome


#open output files
open(OUT, ">$opt{out}" ) or die "cannot write to $opt{out} file";

#for each chromosome read wig file and store locations in hash
foreach my $chromo (sort {$a cmp $b} keys %chromosomes){ #for each chromosome

	open(ROI, "$opt{regions}" ) or die "cannot open $opt{regions}";

	#reset variables for each chromo
	@starts = (); #arrays containing start and end positions of each region
	@ends = ();
	@windows = ();

	#for all regions in peak file
	while(my $line = <ROI>){
		chomp ($line);
		next if $line =~ m/#/;
		my ($cReg, $sReg, $eReg) = (split /\t/, $line)[
					      $opt{chrcol},
					      $opt{startcol}, 
					      $opt{endcol}
					      ]; #Chr Start Stop Number-of-hits
			
		if ($cReg eq $chromo) #if region is in current working chromosome
		{		
				push (@starts,$sReg); #add start and end positions
				push (@ends,$eReg); #add start and end positions
				my $final = ((($eReg - $sReg)-$opt{window})/$opt{step})+1;
				push (@windows,$final);	
		}

	}
	close ROI;

	
	
################################separate overlapping regions

	my @startList = ();
	my @endList = ();
	my @windowList = ();
	push @startList, [ @starts ]; #put arrays into AoA
	push @endList, [ @ends ];
	push @windowList, [ @windows ];

	my $count = 0;
	my $finished = 0;	

	while ($finished == 0){ #loop until no new temp arrays created
		my @sTempArray = ();
		my @eTempArray = ();
		my @wTempArray = ();

		my $totalRegions = @{$startList[$count]}; #total in current list
		my $prevEnd = $startList[$count][0]; #set prev end to start of list
		$finished = 1;
		for (my $i=0;$i<$totalRegions;$i++){
			if ($startList[$count][$i] < $prevEnd){ #if region overlaps previous
				push (@sTempArray, $startList[$count][$i]); #put overlapping elements in new array
				push (@eTempArray, $endList[$count][$i]);
				push (@wTempArray, $windowList[$count][$i]);
				delete $startList[$count][$i]; #remove overlapping elements from original array
				delete $endList[$count][$i];
				delete $windowList[$count][$i];
				$finished = 0; #not finished loop as new list created
			}
			else{
				$prevEnd = $endList[$count][$i]; #set previous end if region remains in current list
			}
		}
		
		push @startList, [ @sTempArray ]; #store temp array in AoA
		push @endList, [ @eTempArray ];
		push @windowList, [ @wTempArray ];
		$count ++; #move to the next list
	}



#################

	my $all = @startList-1;
	for (my $a=0;$a<$all;$a++){ #loop through lists
		#variables for each list
		my $totalRegions = @{$startList[$a]};
		my $count = 0;	
		my $regionStart = $startList[$a][$count];
		my $regionEnd = $endList[$a][$count];
		my $chrFound = 0;
		my $firstRead = 1;
		%norm = ();

		open(NORM, "$opt{norm}" ) or die "cannot write to $opt{norm} file";
		while(<NORM>){
			next if /track/; #ignore wiggle track
			chomp;
	
			my($c, $s, $e, $v) = (split /\s/)[
						      $opt{chrcol},	
						      $opt{startcol}, 
						      $opt{endcol}, 
						      $opt{hitcol}
						      ]; # Start Stop Number-of-hits

			if ($c eq $chromo){ #if region is in current working chromosome
			
				$chrFound = 1;	#chr has been found in input file and currently being read

				#find the first wig read that intersects with a peak region and set all previous regions to score 0
				if($s>=$regionEnd){			
					while ($s>=$regionEnd && $count < $totalRegions){ #if count = total regions then no more regions
						if (!exists($norm{$regionStart}))
						{					
							$norm{$regionStart} ||= 0; #if new region then set start to 0
						}
						#perform sliding window analysis over region
						&window ($chromo,$startList[$a][$count],$endList[$a][$count],$windowList[$a][$count]);
						$count ++;	
						while (!exists $startList[$a][$count]&&$count < $totalRegions){
							$count ++; ##skip past deleted regions in list
						}
						if ($count < $totalRegions){
							$regionStart = $startList[$a][$count]; #next region
							$regionEnd = $endList[$a][$count];
							%norm = ();
							$firstRead = 1;
						}
					}
					
				}
				#if read start is within region
				if ($s>=$regionStart && $count < $totalRegions){
					if ($firstRead == 1){ #if first read in a new region
							$norm{$regionStart} = 0; #set 0 from start of region to first read
							$firstRead = 0;
					}
					$norm{$s} ||=0;
					$norm{$s} += $v; #set read start

					#set read end
					if ($e<$regionEnd){
						$norm{$e} = 0; #within region
					}
					else{
						#perform sliding window analysis over region
						&window ($chromo,$startList[$a][$count],$endList[$a][$count],$windowList[$a][$count]);
						$count ++;
						while (!exists $startList[$a][$count]&&$count < $totalRegions){
							$count ++; ##skip past deleted regions in list
						}
						if ($count < $totalRegions){
							$regionStart = $startList[$a][$count]; # get next region if one exists
							$regionEnd = $endList[$a][$count];
							%norm = ();
							$firstRead = 1;

							while ($e > $regionEnd && $count < $totalRegions){ #check if current read overlaps new regions
								$norm{$regionStart} ||= $v; #set region to read score
								#perform sliding window analysis over region
								&window ($chromo,$startList[$a][$count],$endList[$a][$count],$windowList[$a][$count]);
								$count ++;
								while (!exists $startList[$a][$count]&&$count < $totalRegions){
									$count ++; ##skip past deleted regions in list
								}
								if ($count < $totalRegions){
									$regionStart = $startList[$a][$count]; #next region
									$regionEnd = $endList[$a][$count];
									%norm = ();
									$firstRead = 1;
								}
							}

							if ($e >= $regionStart && $count < $totalRegions){ #if read ends mid region set region start to read score
								$norm{$regionStart} = $v;
								$norm{$e} = 0; #set end of read
								$firstRead = 0;
							}
						}
					}
				}
				elsif ($e>$regionStart && $count < $totalRegions){ #if read start before region start but still intersects
					$norm{$regionStart} ||=0;
					$norm{$regionStart} += $v; #set read start as start of region
				
					if ($firstRead == 1){ #if first read in a new region
						$firstRead = 0;
					}
					#set read end
					if ($e<$regionEnd){
						$norm{$e} = 0; #within region								
					}
					else{
						#perform sliding window analysis over region
						&window ($chromo,$startList[$a][$count],$endList[$a][$count],$windowList[$a][$count]);
						$count ++;
						while (!exists $startList[$a][$count]&&$count < $totalRegions){
							$count ++; ##skip past deleted regions in list
						}
						if ($count < $totalRegions){
							$regionStart = $startList[$a][$count]; # get next region if one exists
							$regionEnd = $endList[$a][$count];
							%norm = ();
							$firstRead = 1;

							while ($e > $regionEnd && $count < $totalRegions){ #check if current read overlaps new regions
								$norm{$regionStart} ||= $v; #set region to read score
								#perform sliding window analysis over region
								&window ($chromo,$startList[$a][$count],$endList[$a][$count],$windowList[$a][$count]);
								$count ++;
								while (!exists $startList[$a][$count]&&$count < $totalRegions){
									$count ++; ##skip past deleted regions in list
								}
								if ($count < $totalRegions){
									$regionStart = $startList[$a][$count]; #next region
									$regionEnd = $endList[$a][$count];
									%norm = ();
									$firstRead = 1;
								}
							}

							if ($e >= $regionStart && $count < $totalRegions){ #if read ends mid region set region start to read score
								$norm{$regionStart} = $v;
								$norm{$e} = 0; #set end of read
								$firstRead = 0;
							}
						}
					}
				}
				if ($count == $totalRegions){
					last;  #end loop if all peak regions already covered
				}
			}
			elsif ($chrFound == 1){
				last; #end loop if all reads in chr have been read (assumes input file is sorted by chr)
			}
		}	

		#if there are still regions left
		if ($count < $totalRegions){
			#perform sliding window analysis over region
			&window ($chromo,$startList[$a][$count],$endList[$a][$count],$windowList[$a][$count]);

			#set any remaining regions to start=0, end=-1	
			for (my $i=$count+1; $i<$totalRegions;$i++){
				%norm = ();
				if (exists $startList[$a][$i]){ #skip past deleted regions in sublist
					$norm{$startList[$a][$i]} ||= 0;
					#perform sliding window analysis over region
					&window ($chromo,$startList[$a][$i],$endList[$a][$i],$windowList[$a][$i]);
				}
			}			
		}	
	}
	close NORM;
}
close OUT;



###############sub routines !!
### would it be easier to set all ends to -1 in one step?



#sliding window analysis
sub window {	

	my $chromo = $_[0]; #chromosome number
	my $rStart = $_[1]; #start of region
	my $rEnd = $_[2]; #end of region
	my $rWindows = $_[3]; #no of windows


	my $baseScore = 0; #hits at each base
	my $windowStart = $rStart; #start of window
	my $windowEnd = $windowStart + $opt{window}-1; #end of window
			 		
	print OUT "$chromo\t$rStart\t$rEnd";	

	for (my $w=1;$w<=$rWindows;$w++){ #while not final window
			
		$score = 0;
		
		#sliding window analysis
		for (my $b=$windowStart; $b<=$windowEnd; $b ++)
		{
							
			#set the baseScore from first base or a change in score
			if (exists $norm{$b}){
					$baseScore = $norm{$b};  
			}


			if ($b == $windowStart + $opt{step})
			{
				$norm{$b} = $baseScore; #set base score for start of next window (unless final window)
			}
			$score += $baseScore; #add to the score for the window
			#print "$rStart\t$rEnd\t$chromo\t$b\t$baseScore\n";
		}
	
		my $avgScore = $score / $opt{window}; #calculate reads per base
		print OUT ("\t$avgScore"); #print outputs

		#set new window positions
		$windowStart += $opt{step}; 
		$windowEnd = $windowStart + $opt{window} -1;
	}

	print OUT "\n";
}
