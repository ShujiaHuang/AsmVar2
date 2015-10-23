# Author : Shujia Huang
# Date   : 2015-10-23
#
# This program is used for extracting varant information.
#!/usr/bin/perl
use strict;
use warnings;
use File::Basename qw/dirname/; 
use lib dirname($0)."/lib"; 
use AsmvarVCFtools;

die qq/perl $0 [vcf] > output\n/ if @ARGV == 0;

my ($vcf_input_file) = @ARGV;

open I, ($vcf_input_file =~ /\.gz$/) ? "gzip -dc $vcf_input_file|" : 
                                        $vcf_input_file or die "[ERROR] $!\n";
while (<I>) {

    chomp;
    next if /^##/;
    
    if (/^#CHROM/) {
        my @col = split;
        print join("\t", @col[0,1,3,4], "SVTYPE\tSVSIZE", @col[9..$#col]), "\n";
    } else {
        my @col = split;
        next if $col[6] ne 'PASS';

        my $svtype = AsmvarVCFtools::GetDataInSpInfo('SVTYPE', \$col[7]); 
        my $svsize = AsmvarVCFtools::GetDataInSpInfo('SVSIZE', \$col[7]); 
        for (my $i = 9; $i < @col; ++$i) {
            my $gt = (split /:/, $col[$i])[0];
            if ($gt !~ /\./) {
                $col[$i] = (split /\//, $gt)[0] + (split /\//, $gt)[1];
            } else {
                $col[$i] = '-';
            }
        }
        print join("\t", @col[0,1], 
                   length($col[3]), 
                   length((split /,/, $col[4])[0]), 
                   $svtype, 
                   $svsize, 
                   @col[9..$#col]), "\n";
    }
}
close I;
