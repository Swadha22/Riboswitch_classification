use warnings;
open(FH,"riboswitches.csv") || die "Can not open the file\n";
open(FHR, "Rfam.txt") || die "Can not open the file\n";
my $data;
my @comb=qw(AA AC AG AU CA CC CG CU GA GC GG GU UA UC UG UU A C G U);
my $matrix_first_row = join("\t",@comb);
while(<FHR>)
{$data .= $_}
my @rfam = split(/>/,$data);
@rfam = sort @rfam;
while(<FH>)
{
	if(/^common/)
	{next;}
	chomp;
	($a,$b)=split(/\t+/);
	$hash{$b}=$a;
}
foreach my $key(sort keys %hash)
{
	#print "$key\t$hash{$key}\n";
}
#Task1=>To print seqs of same ids to 16 difffrent different files.
foreach my $ID(sort keys %hash)
{
	my $fh=ID;
	open(ID,">$ID.txt");
	push(@files, "$ID.txt");
	foreach my $entry(@rfam)
	{
		%hash=();
		if($entry =~ /^$ID/)
		{
		print ID ">$entry";
		}
	}
	close(ID);
}
#End of the task1.
@files = sort @files;
#print "@files", "\n";
foreach my $file(@files)
{
	my %hash= ();
	my $fdata = '';
	my $fh = FH;
	my $op = OP;
	my $line ='';
	my $flag =0;
	open($fh,$file);
	$file =~ s/\.txt//g;
	open($op,">$file.freq.txt");
	print OP "# File Name: $file.freq.txt\n\n";
#If you want to give any other word in place of 'BACTERIA', edit here!
	print OP "BACTERIA\t$matrix_first_row\n";
	while(<$fh>)
	{
	$fdata .= $_;
	}
	@fdata = split('>',$fdata);
	for(my $i=1; $i<scalar(@fdata); $i++)
	{
		$fdata[$i] =~ s/^(\w+\d+.+$)//m;
		$header = $1;
		$fdata[$i] =~ s/\s*//g;
#Following hash table is having fasta header as key and  the sequence as value for all the sequences #in the file.
		$hash{$header}= $fdata[$i]; 
#print "$header$hash{$header}\n"    #To print header and seq for each file.
	}
	foreach my $key(sort keys %hash)
	{				#Variables initialized for each sequence.
		my %hash_dinuc = ();
		my %hash_mono_nuc = ();
		my $dna=$hash{$key};
		my $sum1 = 0;
		my $sum2 = 0;
		my $freq =0;
		my @tmp1=();
		my @tmp2 = ();
#Task2=>Dinucleotides Frequency Calculations.
		for(my $i=0; $i<length($dna)-1; $i++)
		{
			my $dinucleotide = substr($dna,$i,2);
#unless condition to avoid repetition of dinucletides.
			unless(exists $hash_dinuc{$dinucleotide})
			{
				$count1 = () = $dna =~ /$dinucleotide/g;
				if($count1>0)
				{		
				$hash_dinuc{$dinucleotide}= $count1; #another hash table for each sequence.
				}
				else
				{
				$hash_dinuc{$dinucleotide}=0;
				}
				$sum1 = $sum1 +$count1;			
			}
		}
	foreach my $key(sort keys %hash_dinuc)
	{
		my $count = $hash_dinuc{$key};
		my $percent = ($hash_dinuc{$key}/$sum1); 
		   $percent = sprintf("%.2f",$percent);  #Calculating the freq.
			push(@tmp1,$key);
			push(@tmp2,$count);     #output.
	}
#End of Task2.
#Task3=>Mononucleotides Frequency Calculations.
	for(my $i=0; $i<length($dna); $i++)
	{
		my $single_nucleotide = substr($dna,$i,1);
#unless condition to avoid repetition of mononucletides.
		unless(exists $hash_mono_nuc{$single_nucleotide})
			{
			$count2 = () = $dna =~ /$single_nucleotide/g;
			
			$hash_mono_nuc{$single_nucleotide}= $count2; #another hash table for each sequence.
			$sum2 = $sum2 +$count2;			
			}
	}
foreach my $key(sort keys %hash_mono_nuc)
	{
		my $count = $hash_mono_nuc{$key};
		my $percent = ($hash_mono_nuc{$key}/$sum2);
		$percent = sprintf("%.2f",$percent);   #Calculating the freq.
			push(@tmp1,$key);
			push(@tmp2,$count);         #output.
	}
#End of Task3.
my $output2 = join("\t",@tmp2);

#print OP "$key\t$output1\n";
print OP "$key\t$output2\n";
}
print OP "########End of file $file###########\n";
}
exit;
