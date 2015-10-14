# Author : Shujia Huang
# Date   : 2014-05-29 21:01:35
#!/use/bin/perl
use strict;
use warnings;

my ($vcf) = @ARGV;

open I, ( ( $vcf =~ /\.gz$/ ) ? "gzip -dc $vcf |" : $vcf ) or die "Cannot open $vcf : $!\n";
while ( <I> ) {

	chomp;
	if (/^#/) { print "$_\n"; next; }
	my @tmp = split;
	my @tag = split /:/, $tmp[8];
	my %fmat;
	for ( my $i = 0; $i < @tag; ++$i ) { $fmat{$tag[$i]} = $i; }
	next if $tmp[6] =~ /FALSE/;

	my @ales = split /,/, $tmp[4];
	my $isKnown = 0;
    $isKnown = 1 if ( $tmp[2] =~ m/^rs/ );

	my %hit = ();
    my $hitN= 0;
    my $tot = 0;
	my $fail= 0;

	for ( my $i = 9; $i < @tmp; ++$i ) {

		my @samp = split /:/, $tmp[$i];
        ++$tot;

		next if @samp < $fmat{'VT'} + 1 or @samp < $fmat{'QR'} + 1;
		next if (@samp == 1) or ($samp[$fmat{'VT'}] !~ /DEL/ and $samp[$fmat{'VT'}] !~ /INS/);

        $samp[$fmat{'VT'}] = 'INS' if $samp[$fmat{'VT'}] eq 'SINS';
        $samp[$fmat{'VT'}] = 'DEL' if $samp[$fmat{'VT'}] eq 'SDEL';
       
		my @reg = split /,/, $samp[$fmat{'QR'}]; 
		$fail   = 1 if @reg > 1; # Variant Type changed!

		my $k = $samp[$fmat{'TR'}] . ":" . $samp[$fmat{'VT'}] . ":" . $samp[$fmat{'VS'}];
		if (($samp[$fmat{'QR'}] ne ".") and ($samp[$fmat{'GT'}] ne './.')) {
			++$hit{$k};
			++$hitN;
		}
	}
	my @val = values %hit;
	print "$_\n" if ($isKnown and (@val == 1 && $val[0] >= 2 && $hitN < $tot && !$fail)); 
}
close I;

