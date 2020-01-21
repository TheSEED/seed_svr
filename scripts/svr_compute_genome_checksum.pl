

# Usage: svr_compute_genome_checksum contigs-file

use strict;
use Digest::MD5 'md5_hex';

my $ctx = Digest::MD5->new();
my $cur;
my $len;
my @sums;
while (<>)
{
    if (/^>(\S+)/)
    {
	if ($cur)
	{
	    push(@sums, [$cur, $len, $ctx->hexdigest]);
	}
	$cur = $1;
	$len = 0;
    }
    else
    {
	s/\s+//g;
	$ctx->add(uc($_));
	$len += length($_);
    }
}
if ($cur)
{
    push(@sums, [$cur, $len, $ctx->hexdigest]);
}

my $genome_sig = md5_hex(join("", sort map { $_->[2] } @sums));
print "genome_checksum\t$genome_sig\n";

foreach my $ent (sort { $a->[0] cmp $b->[0] } @sums)
{
    print join("\t", "contig_checksum", $ent->[0], $ent->[2]), "\n";
    print join("\t", "contig_length", $ent->[0], $ent->[1]), "\n";
}

