#!/usr/bin/awk -f

#chr7    127471196  127472363  Pos1  0  +

BEGIN {
	# Define field seperators
	FS"\t"; 
	OFS="\t"
}

# Skip the header
/NR==1/ {next;}

# Skip LRG entries
/LRG.*/ {next;}

# For all non-comment lines
{
	chromosome = "chr"$5;
	start      = $3-1;
	end        = $4;
	name       = $2;
	## Write out the line
	print chromosome, start, end, name;
}
