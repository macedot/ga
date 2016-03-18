#!/usr/bin/perl
################################################################################
################################################################################
use strict;
use Time::HiRes qw( time );
use POSIX qw/ceil floor/;
use threads;
use threads::shared;
use Thread qw(async);
use List::Util qw[min max sum];
#use Math::BigInt;
use Data::Dumper;
################################################################################
# Parameters
################################################################################
my $iterationnumber 	= 100; 							# number of iterations
my $crossoverprob 		= 0.3; 								# crossover probability
my $mutationprob 		= 0.012; 							# mutation probability
my $fitnesspow 			= 4.0; 								# pow value for fitness formula
my $indnumber 			= 90; 								# number of individuals (solutions)
my $holdoutval 			= 0.8; 								# percentage of inputs that will be used for cross validation
my $infile 				= "./usd_brl_big_pocket0.1.txt"; 	# indicators results input file
my $initconst 			= 0.4; 								# initialize: chance of bit=1
my $avgenvelope 		= 0.25; 							# initialize: chance of bit=1
################################################################################
# Variables
################################################################################
my $linecount = -1;
my $maxval = 0;
my $end;
my $indsize; # Size of individuals - number of genes
my $crossresultssize = 0;
my $trainingsize = 0;
my $loadedresultssize = 0;
my $currentiteration = 1; # current iteration
# merit function variables
my $param_c = 1;
my $param_l = 1;
my $idxFitness;
my $idxFitnessAdj;
my $idxUsedind;
my @training;
my @crossresults;
my @loadedresults;
my @indINT;
################################################################################
################################################################################
my $BRUTE_FORCE = 0;
my $MAX_THREADS = 4;
my $DEBUG_MODE = 0;
my $MAX_SUBJECT_VALUE = 0; ## RESET AT INPUT FILE READING!!!
################################################################################
################################################################################
#system(clear);
my $start;
$start = time();
print "GA - START [$start]\n" if $DEBUG_MODE;
################################################################################
# Reading parameters
################################################################################
# Reading input file
if ( $ARGV[0] ){
	$infile = $ARGV[0];
}
# Reading population size. If empty, default value, $indnumber, will be used
if ( $ARGV[1] ){
	$indnumber = $ARGV[1];
}
# Reading number of iterations. If empty, default value, $iterationnumber, will be used
if ( $ARGV[2] ){
	$iterationnumber = $ARGV[2];
}
# Reading avg envelope. If empty, default value, $avgenvelope, will be used
if ( $ARGV[3] ){
	$avgenvelope = $ARGV[3];
}
################################################################################
# Printing GA conditions
################################################################################
sub printgaconditions {
	#return unless $DEBUG_MODE;
	print ">>> GA conditions\n";
	print "Number of iterations ..........: $iterationnumber\n";
	print "Crossover probability .........: $crossoverprob\n";
	print "Mutation probability ..........: $mutationprob\n";
	print "Holdout probability ...........: $holdoutval\n";
	print "Fitness pow ...................: $fitnesspow\n";
	print "Population size ...............: $indnumber\n";
	print "Number of genes ...............: $indsize\n";
	print "Input file ....................: $infile\n";
	print "Total data set size ...........: $loadedresultssize\n";
	print "Training data set size ........: $trainingsize\n";
	print "Cross validation data set size : $crossresultssize\n";
	print "Avg envelope ..................: $avgenvelope\n";
	print "Init const ....................: $initconst\n";
	print "\n";
}
################################################################################
################################################################################
sub printholdoutconditions {
	return unless $DEBUG_MODE;
	print "Training data set size:\t\t$trainingsize\n";
	print "Cross validation data set size:\t$crossresultssize\n";
}
################################################################################
# populates a matrix with the indicator results calculated in the 
# scalar(@fields) = 1 + $#fields
################################################################################
sub popresults {
#	my $line;
	my @fields = ();
	my $realize = 0;

	@loadedresults = ();
	open(INFILE,$infile) or die("Cant open file:$!");
	while (<INFILE>){
		my $line = $_;
		chomp($line);
		@fields = split ',', $line;
		my @origFields = @fields;
		# get realized value;
		$realize = $fields[-1];

		# print "\n\n";
		# print join('', @fields),"|$realize|",scalar(@fields),"\n";

		# overwrite with TTCI (The Trusted Coin Estimator) at ACTUAL LAST COLUMN;
		$fields[$#fields] = rand() > 0.5 ? 1 : 0;

		# print join('', @fields),"|$realize|",scalar(@fields),"\n";

		# Creates "Do-Contra" estimator as a NEW column;
		my $val = 0;
		foreach (0..30) { $val = $val | $fields[$_]; };
		$fields[scalar(@fields)] = $val ? 0 : 1;

		# print join('', @fields),"|$realize|",scalar(@fields),"\n";

		# Restores realized value as a NEW COLUMN at the end;
		$fields[scalar(@fields)] = $realize;

		# print join('', @fields),"|$realize|",scalar(@fields),"\n";

		unless ($val) {
			print scalar(scalar(@fields)),"|",$#fields,"\n";
			print "\n";
			print join('', @fields), "|",scalar(@fields),"|$val\n";
			# exit(0);
		}

		# @fields = @origFields;

		$linecount++;
		for my $colcount (0..$#fields){
			$loadedresults[$linecount][$colcount] = $fields[$colcount];
		}
	}
	close(INFILE);
	# Defining individual size based on input file
	# $indsize = $#fields;
	$indsize = $#fields;
	print "indsize -> {{$indsize}}\n";
	# RELATIVE TO IND ARRAY!!!
	# RESULT MUST SUM $indsize!!
	# pos == 0 -> subject as integer;
	$idxFitness    = 1;
	$idxFitnessAdj = 2;
	$idxUsedind    = 3;
	my $thePow = $#fields;
	$MAX_SUBJECT_VALUE = (1 << $thePow);
	$MAX_SUBJECT_VALUE--;
	print "{{$MAX_SUBJECT_VALUE}}",($MAX_SUBJECT_VALUE - 33554431),"\n\n";
	# Defining loaded results size
	$loadedresultssize = scalar(@loadedresults);
}
################################################################################
# Holdout method to fullfill traning and crossvalidation with loadedresults
################################################################################
sub holdout_OLD {
	@training = ();
	@crossresults = ();
	my $trainingrow = 0;
	my $crossvalrow = 0;
	my $myrand = 0;
	for my $row (0..$#loadedresults){
#		next unless $BRUTE_FORCE || $loadedresults[$row][$indsize];
		$myrand=rand();
#		next unless $BRUTE_FORCE || $loadedresults[$row][$indsize];
		if ($BRUTE_FORCE || $myrand > $holdoutval){
#		if ($BRUTE_FORCE || rand() > $holdoutval){
			@{$training[$trainingrow++]} = @{$loadedresults[$row]};
		}
#		next unless $loadedresults[$row][$indsize];
		if ( $myrand <= $holdoutval ){
#		if ( rand() <= $holdoutval ){
			@{$crossresults[$crossvalrow++]} = @{$loadedresults[$row]};
		}
	}
	# ASSERT
	@{$training[0]}     = @{$loadedresults[floor(rand() * $#loadedresults)]} unless $trainingrow;
	@{$crossresults[0]} = @{$loadedresults[floor(rand() * $#loadedresults)]} unless $crossvalrow;
	# Defining training data set size
	$trainingsize = scalar(@training);
	# Defining cross validation data set size
	$crossresultssize = scalar(@crossresults);
}

sub holdout {
	@training = ();
	@crossresults = ();

	my $trainingrow = 0;
	my $crossvalrow = 0;

	#my $splitPos = floor(rand() * $#loadedresults);
	my $splitPos = floor($holdoutval * $#loadedresults);

	for my $row (0..$splitPos){
		@{$training[$trainingrow++]} = @{$loadedresults[$row]};
	}

	for my $row ($splitPos..$#loadedresults){
		@{$crossresults[$crossvalrow++]} = @{$loadedresults[$row]};
	}

	# ASSERT
	@{$training[0]}     = @{$loadedresults[floor(rand() * $#loadedresults)]} unless $trainingrow;
	@{$crossresults[0]} = @{$loadedresults[floor(rand() * $#loadedresults)]} unless $crossvalrow;

	# Defining training data set size
	$trainingsize = scalar(@training);
	# Defining cross validation data set size
	$crossresultssize = scalar(@crossresults);
}
################################################################################
# Calculates performance of each indicator if used separately
################################################################################
sub indicatorstest {
	#return if $BRUTE_FORCE;
	
	my $totalrights = 0;
	my $totalratio = 0;
	print "Correct predictions per indicator\n" if $DEBUG_MODE;
	print "Indicator\tTotal ($loadedresultssize)\n" if $DEBUG_MODE;

	for my $indicator (0..($indsize-1)){
		$totalrights = 0;
		for my $row (0..$#crossresults){
			# Compare if it predicted correctly
			if ( $crossresults[$row][$indsize] == $crossresults[$row][$indicator] ){
				$totalrights++;
			}
		}
		$totalratio = $totalrights / $loadedresultssize;
		$totalratio = $totalratio * 100;
		printf "\t$indicator\t$totalrights\t%.2f%\n",$totalratio if $DEBUG_MODE;
	}
	print "\n";
}
################################################################################
# Print arrays
################################################################################
sub printarray {
	return unless $DEBUG_MODE;
	my @array = @_;
	for my $indn (0..$#array){
		print "\tLine\t$indn:\t";
		printf "%s", dec2bin($array[$indn][0]);
		printf "|%5s", $array[$indn][$idxFitness];
		printf "|%5s", $array[$indn][$idxFitnessAdj];
		printf "|%2s", $array[$indn][$idxUsedind];
		print "\n";
	}
	print "\n";
}
################################################################################
# Initialize subjects
################################################################################
sub BRUTE_FORCE {
	print ">>> BRUTE FORCE <<<\n";
	# print dec2bin(1234567890), "\n";
	# print length(dec2bin(1234567890)), "\n";
	# print NumberOfSetBits(1234567890), "\n";
	# exit 0;
	my $min = 1;
	#my $max = 100000; ## 2^25 - 1
	my $max = $MAX_SUBJECT_VALUE; ## 2^25 - 1
	my $first = $min;
	my $limit = $max;
	my $chunk = min($limit, 1 + int(($limit - $first) / $MAX_THREADS));
	my @threadList = ();
	#print "[$first, $chunk, $limit]\n";
	while ($first < $limit) {
		my $last = min($first + $chunk, $limit);
		my $thd = async {
			print "[$first, $chunk, $last]\n";
			my $bestInd = 0;
			my $bestIndNum = $indsize+1;
			my $bestFit = $trainingsize * $trainingsize;
			for (my $subject = $first; $subject <= $last; $subject++) {
				my $usedind = NumberOfSetBits($subject);
#				my $fitnessValue = subjectFitness_avg($subject, $usedind);
				my $fitnessValue = subjectFitness_merit($subject, $usedind);
				if ($bestFit > $fitnessValue) {
					$bestFit = $fitnessValue;
					$bestInd = $subject;
					$bestIndNum = $usedind;
				} elsif ($bestFit == $fitnessValue) {
					if ($bestIndNum > $usedind) {
						$bestInd = $subject;
						$bestIndNum = $usedind;
					}
				}
			}
			print "bestInd = $bestInd [" . dec2bin($bestInd) . "] bestFit = $bestFit; bestIndNum = $bestIndNum\n";
		};
		push(@threadList, $thd);
		$first = $last
	}
	map { $_->join(); } @threadList;
	printf("Elapsed time: %.2f\n", time() - $start);
	exit(0);
}
sub initialize {
	&BRUTE_FORCE() if $BRUTE_FORCE;
	print ">>> Initializing individuals\n" if $DEBUG_MODE;
	my $const = $initconst;
	for my $indn (0..$indnumber-1){
		$indINT[$indn][0] = int(rand() * $MAX_SUBJECT_VALUE);
		#$indINT[$indn][0] = 0;
		#$indINT[$indn][0] = $MAX_SUBJECT_VALUE;
		$indINT[$indn][$idxFitness] = 0;
		$indINT[$indn][$idxFitnessAdj] = 0;
		$indINT[$indn][$idxUsedind] = 0;
		$const = 1.0 - $const;
	}
}
################################################################################
# Normalize fitness column for next generation selection through roulette wheel method
################################################################################
sub fitnessadj_merit {
	my $diffsum=0.00;
	my $fitnessval=0.00;
	my $raw_fitadj=0.00;
	# Getting the min of fitness column
	my $minFit = 999999999;
	foreach(0..$#indINT) {
		if ($indINT[$_][$idxFitness] < $minFit) {
			$minFit = $indINT[$_][$idxFitness];
		}
	}
	for my $row (0..$indnumber-1){
		$raw_fitadj = ( $indINT[$row][$idxFitness] - $minFit );
		$indINT[$row][$idxFitnessAdj] = $raw_fitadj ** $fitnesspow;
		$indINT[$row][$idxUsedind] = usedindnum($row);
		# Getting the sum of adjusted Fitness with pow;
		$diffsum += $indINT[$row][$idxFitnessAdj];
	}
	return $diffsum;
}
################################################################################
# Normalize fitness column for next generation selection through roulette wheel method
################################################################################
sub fitnessadj {
	my $diffsum=0.00;
	my $fitnessval=0.00;
	my $diffabs=0.00;
	# Getting the max of fitness column
	my $maxFit = 0;
	foreach(0..$#indINT) {
		if ($indINT[$_][$idxFitness] > $maxFit) {
			$maxFit = $indINT[$_][$idxFitness];
		}
	}
	# abs of diff between fitness and its max
	for my $row (0..$indnumber-1){
		$diffabs = ($maxFit - $indINT[$row][$idxFitness]);
		$indINT[$row][$idxFitnessAdj] = $diffabs ** $fitnesspow;
		$indINT[$row][$idxUsedind] = usedindnum($row);
		# Getting the sum of adjusted Fitness with pow;
		$diffsum += $indINT[$row][$idxFitnessAdj];
	}
	return $diffsum;
}
################################################################################
# Returns number of used indicators by individual which row is passed as first parameter
################################################################################
sub usedindnum {
	my $indrow = $_[0];
	return NumberOfSetBits($indINT[$indrow][0]);
}
################################################################################
# Prediction function for a given indicator and result row. AVG method.
# Returns 2 if number of indicators used equals 0
# Returns 17 in case of errors in the logic
################################################################################
sub prediction_avg {
	my $arTest  = $_[0];
	my $arInd   = $_[1];
	my $usedind = $_[2];
	# print join('', @{$arTest}),"\n";
	# print join('', @{$arInd}),"\n";
	# Sum of used indicators results
	my $indSum = 0;
	for my $col (0..($indsize-1)){
		# next unless @{$arInd}[$col];
		# next unless @{$arTest}[$col];
		$indSum += ( @{$arInd}[$col] * @{$arTest}[$col] );
	}
	return $indSum unless $indSum;
	my $avg = $indSum / $usedind;
	#print "$avg = $indSum / $usedind;\n";
	# Mutable envelope
#	my $avgenvelope = 1 / $usedind;
	# The 2 lines below are equivalent
#	my $avgenvelope = ceil($usedind / 2.0) / $usedind;
#	my $avgenvelope = ceil($avgconst * $usedind) / $usedind;
##	my $avgenvelope = 0.5;
	if ( $avg >= $avgenvelope ) {
		return 1;
	} 
	return 0;
}
################################################################################
################################################################################
sub subjectFitness_merit {
	my $subject = shift;
	my $usedind = shift;
	my $fitnessValue = 0.00;
	my %confusionTable;
	$confusionTable{'11'} = 0;
	$confusionTable{'10'} = 0;
	$confusionTable{'01'} = 0;
	$confusionTable{'00'} = 0;
	my $confstrg = 0;
	my $myprediction;
	my @arSubject = split('', dec2bin($subject));
	# For each training row...
	for my $trow (0..$#training){
#		next unless $training[$trow][$indsize];
		# Adds 1 if observed tendency is different than prediction
		$myprediction = prediction_avg($training[$trow], \@arSubject, $usedind);
#		chomp $training[$trow][$indsize];
#		printf "%x\n",$training[$trow][$indsize];
		$confstrg=$training[$trow][$indsize].$myprediction;
		$confusionTable{$confstrg}++;
#		print "[$confstrg]\n";
#		print "prediction=$myprediction\n";
#		print "actual=$training[$trow][$indsize]\n";
	}
	# print "param_c=$param_c\n";
	# print "param_l=$param_l\n";
	# print "confTable[11]=$confusionTable{'11'}\n";
	# print "confTable[01]=$confusionTable{'01'}\n";
	$fitnessValue = ($param_c * $confusionTable{'11'}) - ($param_l * $confusionTable{'01'});
	
#	if ($confusionTable{'11}'} ==0 && $confusionTable{'01'} == 0) {
#		$fitnessValue = 0;
#	}
#	print Dumper(\%confusionTable),"\n";
#	print "fitnessValue=$fitnessValue\n";
#	exit 1;
	return $fitnessValue;
}
################################################################################
################################################################################
sub subjectFitness_avg {
	my $subject = shift;
	my $usedind = shift;
	my $fitnessValue = 0.00;
	my $myprediction;
	my @arSubject = split('', dec2bin($subject));
	# For each training row...
	for my $trow (0..$#training){
#		next unless $training[$trow][$indsize];
		# Adds 1 if observed tendency is different than prediction
		$myprediction = prediction_avg($training[$trow], \@arSubject, $usedind);
		$fitnessValue += ($training[$trow][$indsize] != $myprediction);
	}
	return $fitnessValue;
}
################################################################################
# Fitness function AVG method
################################################################################
sub calcfitness_avg {
	my $first = 0;
	my $limit = $#indINT;
	my $chunk = min($limit, 1 + int(($limit - $first) / $MAX_THREADS));
	my @threadList = ();
	my %theFit :shared;
	%theFit = {};
	my $subject;
	my $usedind;
	while ($first < $limit) {
		my $last = min($first + $chunk, $limit);
		my $thd = async {
			my $fitnessValue;
			for my $idxSubject ($first..$last) {
				# If prediction == 2 (that happens when no indicators are selected)
				# assign worst score possible. "* 2" because diff between 1 and -1 is 2.
				$fitnessValue = (-1) * $trainingsize;
				$subject = $indINT[$idxSubject][0];
				if ($subject > 0) {
					$usedind = NumberOfSetBits($subject);
					$fitnessValue = subjectFitness_merit($subject, $usedind);
				}
				{ lock(%theFit); $theFit{$idxSubject} = $fitnessValue; }
			}
		};
		push(@threadList, $thd);
		$first = $last
	}
	map { $_->join(); } @threadList;
	for my $idxSubject (0..$limit){
		unless (exists($theFit{$idxSubject})) {
			print "\t-> INVALID: $idxSubject\n";
			next;
		}
		$indINT[$idxSubject][$idxFitness] = $theFit{$idxSubject};
	}
	my $max = &fitnessadj_merit;
	&sortarray();
	return $max;
}
################################################################################
#Sort ind array by its last field (adjusted fitness)
################################################################################
sub sortarray {
	my @array = sort cmpfunc @indINT;
	@indINT = @array;
}
################################################################################
################################################################################
sub cmpfunc {
	$a->[$idxFitnessAdj] <=> $b->[$idxFitnessAdj]
	|| $a->[$idxUsedind] <=> $b->[$idxUsedind];
}
################################################################################
# Select new population
################################################################################
sub subjectSelection {
	my $max = $_[0];
	return 0 unless $max;
	my $sum = 0;
	my $curfitadj = 0;
	my $myrand = rand() * $max;
	for my $subject (0..$#indINT){
		$curfitadj = $indINT[$subject][$idxFitnessAdj];
		if ( ($sum <= $myrand) && ($myrand <= ($curfitadj + $sum) ) ) {
			return $subject;
		}
		$sum += $curfitadj;
	}
	return $#indINT;
}
sub selection {
	my $max = $_[0];
	my $curfitadj = 0;
	my @indauxINT = ();
	# loop ($indnumber-2) times to get (same number - 1) subjects for next generation
	for my $count (0..$indnumber-2){
		my $subject = subjectSelection($max);
		#print "\t->INVALID: $subject\n" unless exists $ind[$subject];
		@{$indauxINT[$count]} = @{$indINT[$subject]};
	}
	# Best individual goes automatically to new population
	@{$indauxINT[$#indINT]} = @{$indINT[$#indINT]};
	@indINT = @indauxINT;
	&sortarray();
}
################################################################################
# Crossover
################################################################################
sub crossover {
	my $first = 0;
	my $second = 0;
	my $addone = 0;
	my $pairs = int($indnumber / 2.0);
	# In case of odd number of individuals, add 50% chances of crosssing over pairs shifted by 1
	# ie., in case of 5 individuals, add 50% chances of crossing over paiasr (2-3) and (4-5)
	# without the code below, only pairs (1-2) and (3-4) would have a chance of crossing over
	my $remainder = $indnumber % 2;
	$addone = 1 if $remainder && rand() > 0.5;
	# For possible pairs loop...
	for my $count (1..$pairs) {
		if ( rand() <= $crossoverprob ) {
			$first = (2*$count)-1+$addone;
			$second = (2*$count)-2+$addone;
			my $myrandpos = floor(rand() * ($indsize));
			print "\tCrossover between $first and $second at position $myrandpos\n" if $DEBUG_MODE;
			# From position pos till the last data position (scalar-3) swap values
			my $input1 = $indINT[$first][0];
			my $input2 = $indINT[$second][0];
			my $mask1 = ((0xffffffffffffffff >> 64 * $myrandpos) << 64 * $myrandpos);
			my $mask2 = 0xffffffffffffffff ^ $mask1;
			$indINT[$first][0]  = ($input1 & $mask1) ^ ($input2 & $mask2);
			$indINT[$second][0] = ($input1 & $mask2) ^ ($input2 & $mask1);
		}
	}
	print "\n" if $DEBUG_MODE;
}
################################################################################
# Mutation
################################################################################
sub mutation {
	for my $subject (0..$#indINT){
		for my $col (0..($indsize-1)){
			if ( rand() <= $mutationprob ) {
				print "\tMutation in individual $subject at position $col\n" if $DEBUG_MODE;
				#$ind[$subject][$col] ^= 1;
				my $mask = (1 << $col);
				$indINT[$subject][0] ^= $mask;
			}
		}
	}
	print "\n"  if $DEBUG_MODE;
}
################################################################################
################################################################################
sub bin2dec {
	# my $str = substr("0" x 64 . $_[0], -64);
    # return unpack("Q>", pack("B64", $str));
	return oct("0b" . shift)
}
sub dec2bin {
	my $dec = $_[0];
	my $str = unpack("B64", pack("Q>", $dec));
	my $bin = substr $str, -$indsize;
	return $bin;
}
sub NumberOfSetBits {
	return (dec2bin($_[0]) =~ tr/1//);
}
################################################################################
# Cross Validation
################################################################################
sub crossvalidation_merit {
	my $iteration = $_[0];
	my $chosenindnum = 0;
	my $percentage = 0;
	my $prediction = 0;
	my $crosspercentage = 0;
	my $intSubject = $indINT[-1][0];
	my $fitSubject = $indINT[-1][$idxFitness];
	my $strSubject = dec2bin($intSubject);
	my @arSubject = split('', $strSubject);
	my %confusionTable;
	my $confstrg = 0;
	
	# The Chosen One
	printf "I%d # %s [%d] {%d}", $iteration, $strSubject, $intSubject, $fitSubject;
	# Calculating the sum of the columns of the best individual, ie., number of technical indicators chosen
	$chosenindnum = usedindnum(-1);
	printf " #i: %2d/%2d", $chosenindnum, $indsize;
	# Calculating the number of distinct subjects in the population
	my %cardinal = ();
	for my $subject (0..$#indINT){
		$cardinal{$indINT[$subject][0]}++;
	}
	printf " #s: %2d/%2d", scalar keys %cardinal, $indnumber;
	# Counts how many results are reproduced by the best individual
	for my $crow (0..$#crossresults){
	$prediction = prediction_avg(\@{$crossresults[$crow]}, \@arSubject, $chosenindnum);
			$confstrg=$crossresults[$crow][$indsize].$prediction;
			$confusionTable{$confstrg}++;
	}
	
	my $totalACC = 0;
	foreach my $elem ( sort {$b cmp $a} keys %confusionTable) {
		printf " | %s:%03d", $elem, $confusionTable{$elem};
		$totalACC += $confusionTable{$elem};
	}
	print " [$totalACC]";

	my $rights = $confusionTable{'11'} + $confusionTable{'00'};

	my $totalTPR = $confusionTable{'11'} + $confusionTable{'10'}; # TPR
	my $totalPPV = $confusionTable{'11'} + $confusionTable{'01'}; # PPV

	my $percACC = 0; $percACC = 100 * $rights / $totalACC if $totalACC > 0;
	my $percTRP = 0; $percTRP = 100 * $confusionTable{'11'} / $totalTPR if $totalTPR > 0;
	my $percPPV = 0; $percPPV = 100 * $confusionTable{'11'} / $totalPPV if $totalPPV > 0;
	printf " # ACC: %.2f% TPR: %.2f% PPV: %.2f%", $percACC, $percTRP, $percPPV;

	print "\n";
	print "\n" if $DEBUG_MODE;
}
################################################################################
# Cross Validation
################################################################################
sub crossvalidation_avg {
	my $iteration = $_[0];
	my $chosenindnum = 0;
	my $percentage = 0;
	my $prediction = 0;
	my $crosspercentage = 0;
	my $rights = 0;
	my $intSubject = $indINT[-1][0];
	my $fitSubject = $indINT[-1][$idxFitness];
	my $strSubject = dec2bin($intSubject);
	my @arSubject = split('', $strSubject);
	
	# The Chosen One
	printf "I%d #B: %s [%s] {%d} ", $iteration, $strSubject, $intSubject, $fitSubject;
	# Calculating the sum of the columns of the best individual, ie., number of technical indicators chosen
	$chosenindnum = usedindnum(-1);
	printf " #i: %2d/%2d", $chosenindnum, $indsize;
	# Calculating the number of distinct subjects in the population
	my %cardinal = ();
	for my $subject (0..$#indINT){
		$cardinal{$indINT[$subject][0]}++;
	}
	printf " #s: %2d/%2d", scalar keys %cardinal, $indnumber;
	# Counts how many results are reproduced by the best individual
	# For each cross results record
	if ($intSubject) {
		for my $crow (0..$#crossresults) {
			$prediction = prediction_avg(\@{$crossresults[$crow]}, \@arSubject, $chosenindnum);
			if ( $prediction == $crossresults[$crow][$indsize] ){
				$rights++;
			}
		}
	}
	# for my $crow (0..($trainingsize-1)) {
	# 	$prediction = prediction_avg(\@training, $crow, -1);
	# 	if ( $prediction == $training[$crow][$idxFitness] ){
	# 		$rights++;
	# 	}
	# }
	# for my $crow (0..($loadedresultssize-1)) {
	# 	$prediction = prediction_avg(\@loadedresults, $crow, -1);
	# 	if ( $prediction == $loadedresults[$crow][$idxFitness] ){
	# 		$rights++;
	# 	}
	# }
	$crosspercentage = $rights / $crossresultssize;
	$crosspercentage = $crosspercentage * 100;
	printf " # Correct $rights/$crossresultssize (or %.2f%)\n", $crosspercentage;
	print "\n\n" if $DEBUG_MODE;
}
################################################
##### MAIN ########## MAIN ########## MAIN #####
################################################
&popresults;
print ">>> Input file loaded\n";
#&printarray(@loadedresults);
&holdout;
print ">>> Holdout executed\n\n" if $DEBUG_MODE;
&printgaconditions;
print ">>> Training data\n";
#&printarray(@training);
print ">>> Validation data\n";
#&printarray(@crossresults);
print ">>> Performance of each indicator if used separately\n" if $DEBUG_MODE;
&indicatorstest;
&initialize();
print ">>> First Population\n" if $DEBUG_MODE;
&printarray(@indINT);
print "##########################################\n" if $DEBUG_MODE;
print "### Iteration #1\n\n" if $DEBUG_MODE;
print "I$currentiteration\tCalc Fitness\n" if $DEBUG_MODE;
$maxval = &calcfitness_avg();
print "I1\tPopulation after fitness\n" if $DEBUG_MODE;
&printarray(@indINT);
print "I$currentiteration\tCross validation results\n" if $DEBUG_MODE;
&crossvalidation_merit($currentiteration);
for my $currentiteration (2..$iterationnumber){
	print "##########################################\n" if $DEBUG_MODE;
	print "### Iteration #" . $currentiteration . "\n\n" if $DEBUG_MODE;
##
#	print ">>> Performance of each indicator if used separately\n" if $DEBUG_MODE;
#	&indicatorstest;
##
	print "I$currentiteration\tSelection\n" if $DEBUG_MODE;
	&selection($maxval);
	# print "New population -  after selection\n";
	# &printarray(@ind);
	print "I$currentiteration\tCrossover\n" if $DEBUG_MODE;
	&crossover;
	# print "New population - after crossover\n";
	# &printarray(@ind);
	print "I$currentiteration\tMutation\n" if $DEBUG_MODE;
	&mutation;
	# print "New population -  after mutation\n";
	# &printarray(@ind);
	# &holdout;
	# print ">>> Holdout executed\n\n" if $DEBUG_MODE;
	# &printholdoutconditions;
	print "I$currentiteration\tCalc Fitness\n" if $DEBUG_MODE;
	$maxval = &calcfitness_avg();
	print "I$currentiteration\tNew population after fitness\n" if $DEBUG_MODE;
	&printarray(@indINT);
	print "I$currentiteration\tCross validation with best individual - results\n" if $DEBUG_MODE;
	&crossvalidation_merit($currentiteration);
}
$end = time();
print "##########################################\n" if $DEBUG_MODE;
printf("Elapsed time: %.2f\n", $end - $start)  if $DEBUG_MODE;
print "Exiting...\n" if $DEBUG_MODE;
print "GA - END\n" if $DEBUG_MODE;
warn sprintf "%s\n", $indINT[-1][0];
# The answer is 42
exit 42-42;
