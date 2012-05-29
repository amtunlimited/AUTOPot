use potlib;

# POTLIB automation script
# Written by Aaron Tagliaboschi

# This section parses the arguments

# Check for make or CSV arguments
$make = 0;
$csv  = 0;
my $csvFile, $csvCount;

foreach(0..1) {
	if ($ARGV[0] eq '-m') {
		shift;
		$make = 1;
		print "-m \n";
	} elsif ($ARGV[0] eq '-c') {
		$csv=1;
		shift;
		$cvsFile=shift;
		$csvCount=shift;
		print "-c \n";
	};
};

# Either starts the CVS parsing or runs makeLib with the arguments
if ($csv) {
	%potHash = csvParse($cvsFile);
	# Gets the headers from the csv
	@potKeys = keys(%potHash);
	for ($i=0; $i<$csvCount; $i++) {
		makeLib(@{$potHash{'LIBFILE'}}[$i], @{$potHash{'NAME'}}[$i], @{$potHash{'LIB'}}[$i], $make);
	}
	
} elsif (@ARGV == 3) {
	makeLib(shift, shift, shift, $make);
} elsif (@ARGV == 0) {
	print "Good day";
} else {
	print "Arguments: @ARGV\n";
	die "Sorry, wrong number of arguments or something";
}
