# The main lib-creating sub
# arguments:
#	file: the filename
#	name: the eventual name of the folder, function, and excecutable
#	lib:  the name of the lib-folder

sub makeLib {
	$file = shift;
	$name = shift;
	$lib  = shift;
	$make = shift;
	print("Now making $name\n");
	
	# Copy the particular library and remane it the function name
	system("cp", $lib, $name, "-r");
	# And copy the lib source into that
	system("cp", $file, $name."/src/potlib.f");

	# Now for the fun part...
	chdir("./".$name);
	# Text replacement!
	
	# First we go through the make file
	fileFindReplace("Makefile", "NAME", $name);
	
	# Then the ML template
	fileFindReplace("src/pot.tm", "FUNC", $name);
	
	# run make if the "-m" was up
	if ($make) {
		print "ta da!\n";
		system("make");
	}
	
	# ...and that's it! Easy, wasn't it? A little too easy...
	
	chdir("../");
}

# Opens a file (fisrt argument) and basically find-and-replaces
# all instances of the second argument with the third
sub fileFindReplace {
	
	$file = $_[0];
	$newFile = $_[0]."~";
	
	open(FILEIN,  "<", $file);
	open(FILEOUT, ">", $newFile);
	
	while (<FILEIN>) {
		$replace = $_;
		$replace =~ s/$_[1]/$_[2]/ge;
		print FILEOUT $replace;
	}
	
	close FILEIN;
	close FILEOUT;
	
	system("rm", $file);
	system("mv", $newFile, $file);
}

# This subroutine parses a CSV file and returns a hash

sub csvParse{
	# Puts the file line by line in the array @csv
	open CSVIN, "<", shift;
	@csv = <CSVIN>;
	close CSVIN;
	
	# Gets the hash keys from the first line
	@keys = split ',', shift(@csv);
	
	# Removes trailing whitespace including \n (gave me a whole heap of trouble)
	foreach (@keys) {
		$_ =~ s/^\s+//;
		$_ =~ s/\s+$//;
	};
	
	# Makes the hash using the keys
	%csvHash = ();
	foreach (@keys) {
		$csvHash{$_} = [];
	}
	
	# Puts the rest of the file in the hash
	$i=0;
	foreach (@csv) {
		@row = split ',', $_;
		$j=0;
		foreach (@row) {	
			$_ =~ s/^\s+//;
			$_ =~ s/\s+$//;
			$csvHash{@keys[$j]}[$i]=$_;
			$j++;
		}
		$i++
	}
	
	return %csvHash;
}

# Now for te windows version. It's mostly copypasta, except 
# for changing the commands and switching the slashes

sub makeLibWin {
	$file = shift;
	$name = shift;
	print("Now making $name\n");
	
	print "1\n";
	# Copy the particular library and remane it the function name
	system("mkdir", $name);
	system("xcopy", "POTLIBWIN\\\*", $name, "/s", "/i");
	
	print "2\n";
	# And copy the lib source into that
	system("copy", $file, $name."\\src\\potlib.f");

	print "3\n";
	# Now for the fun part...
	chdir(".\\".$name);
	# Text replacement!
	
	print "4\n";
	# First we go through the make file
	fileFindReplace("Makefile", "NAME", $name);
	
	print "5\n";
	# Then the ML template
	fileFindReplace("src\\pot.tm", "FUNC", $name);
	
	print "6\n";
	print "ta da!\n";
	system("gmake");
	
	
	# ...and that's it! Easy, wasn't it? A little too easy...
	
	chdir("..\\");
}

# Windows version. Again, just correcting slashes and commands
sub fileFindReplace {
	
	$file = $_[0];
	$newFile = $_[0]."~";
	
	open(FILEIN,  "<", $file);
	open(FILEOUT, ">", $newFile);
	
	while (<FILEIN>) {
		$replace = $_;
		$replace =~ s/$_[1]/$_[2]/ge;
		print FILEOUT $replace;
	}
	
	close FILEIN;
	close FILEOUT;
	
	system("del", $file);
	system("move", $newFile, $file);
}

# Return true for the lib-user's sake
1;
