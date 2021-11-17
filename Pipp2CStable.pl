#! /usr/bin/perl



unless(@ARGV==2) {
  print "input file name and output file name are needed as input\n";
  exit;
}

my $filename = $ARGV[0];
my $output = $ARGV[1];
 

unless(open(FILE_DATA, $filename)) {
    exit;
}

unless (open(OUT, ">$output")) {
    print "can not open file";
    exit;
}

	printf OUT ("%-5s %-5s %-7s %-7s %-7s %-7s %-7s %-7s %-7s %-7s %-7s %-7s %-7s %-7s\n","Id", "Res","HN","HA1","HA2","HB1","HB2","HG1","HG2","HD1","HD2","HE","HH", "N"); 
$HA1 =$HA2 = $HB1 = $HB2= $HN=$HG1 =$HG2=$HD1=$HD2=$HE = $HH =$N ='_';



while(<FILE_DATA>) {
    chomp;
    if ($_ =~ /RES_ID/) {
        @column = split(' ', $_);
        ($res_id) = $column[1];
	#print $res_id, "\n";
    }

    if ($_ =~ /RES_TYPE/) {
        @column = split(' ', $_);
        ($res_type) = $column[1];
	#print $res_type,"\n";
    }

    
    

    if ($_ =~/^   [N|H]/) {
	@column = split('\(', $_);
	@cc = split(' ', $column[0]);
	#print $cc[0], "\t", $cc[1], "\n";

	if($cc[0]=~/\|/) {
	    @tmp = split('\|', $cc[0]);
	    $sign = $tmp[0];
	}
	else {
	    $sign =$cc[0];
	}


	if ($sign =~ /HA/) {
	    if ($sign =~/HA1/) {
		$HA1 = $cc[1];
	    }
	
	    elsif ($sign =~/HA2/) {
		$HA2 = $cc[1];
	    }
	    else {
		$HA1 = $cc[1];
	    }
	}

	if ($sign =~ /HB/) {
	    if ($sign =~/HB1/) {
		$HB1 = $cc[1];
	    }
	
	    elsif ($sign =~/HB2/) {
		$HB2 = $cc[1];
	    }
	    else{
		$HB1 = $cc[1];
	    }
	}

	if ($sign =~ /HG/) {
	    if ($sign =~/HG1/) {
		$HG1 = $cc[1];
	    }
	
	    elsif ($sign =~/HG2/) {
		$HG2 = $cc[1];
	    }
	    else{
		$HG1 = $cc[1];
	    }
	}

	if ($sign =~ /HD/) {
	    if ($sign =~/HD1/) {
		$HD1 = $cc[1];
	    }
	
	    elsif ($sign =~/HD2/) {
		$HD2 = $cc[1];
	    }
	    else {
		$HD1 = $cc[1];
	    }
	}

	if ($sign =~ /HE/) {
	    $HE = $cc[1];
	}
	
	if ($sign =~ /HH/) {
	    $HH = $cc[1];
	}

	if ($sign =~ /^N$/) {
	    $N = $cc[1];
	}

    	if ($sign =~ /HN/) {
	    $HN = $cc[1];
	}
    }

    if ($_ =~ /END_RES_DEF/) {
	printf OUT ("%-5d %-5s %-7.3f %-7.3f %-7.3f %-7.3f %-7.3f %-7.3f %-7.3f %-7.3f %-7.3f %-7.3f %-7.3f %-7.3f\n", $res_id, $res_type, $HN, $HA1, $HA2, $HB1, $HB2, $HG1, $HG2, $HD1, $HD2, $HE, $HH, $N); 
$HA1 =$HA2 = $HB1 = $HB2= $HN=$HG1 =$HG2=$HD1=$HD2=$HE = $HH =$N ='_';
    }
 }   



