use Data::Dumper;
use strict;
use JSON::XS;
use Bio::KBase::AppService::SchedulerDB;
use Bio::P3::Workspace::WorkspaceClientExt;

my $ws = Bio::P3::Workspace::WorkspaceClientExt->new;

#
# Given results from analyze-assembly-data.pl, pull the workspace info again and
# determine input sizes.
#

my $db = Bio::KBase::AppService::SchedulerDB->new;

my @work;
my %data;
while (<>)
{
    chomp;
    my($id, $owner, $size, $elap, $mem) = split(/\t/);
    my($h, $m, $s) = split(/:/, $elap);
    my $time = $h * 3600 + $m * 60 + $s;

    my $val = [$id, $owner, $size, $elap, $mem, $time];
    push(@work, $val);
    $data{$id} = $val;
}

while (@work)
{
    my @batch = splice(@work, 0, 100);

    my $qry = join(", ", );
    my $qs = join(", ", map { "?" } @batch);
    
    my $res = $db->dbh->selectall_arrayref(qq(SELECT id, output_path, output_file, params
					      FROM TaskWithActiveJob
					      WHERE id in ($qs)), undef, map { $_->[0] } @batch);
    for my $e (@$res)
    {
	my($id, $path, $file, $params) = @$e;

	my $dat = eval { decode_json($params); };
	next unless $dat;

	my $comp = 0;
	my $uncomp = 0;
	
	my $pl = $dat->{paired_end_libs};
	my $sl = $dat->{single_end_libs};
	my $sra = $dat->{srr_ids};

	for my $pe (@$pl)
	{
	    accum($pe->{read1}, \$comp, \$uncomp);
	    accum($pe->{read2}, \$comp, \$uncomp);
	}
	for my $se (@$sl)
	{
	    accum($se->{read}, \$comp, \$uncomp);
	}
	for my $srr (@$sra)
	{
	    accum_sra($srr, \$comp, \$uncomp);
	}

	my(undef, $owner, $size, $elap, $mem, $time) = @{$data{$id}};
	
	print join("\t", $id, $owner, $size, $elap, $mem, $time, $comp, $uncomp), "\n";
    }
}

sub accum
{
    my($path, $comp, $uncomp) = @_;
    my @a = $ws->stat($path);
    my $sz = $a[7];
    my $what = $ws->file_is_gzipped($path) ? $comp : $uncomp;
    $$what += $sz;
}

sub accum_sra
{
    my($id, $comp, $uncomp) = @_;
    open(P, "-|", "p3-sra", "--id", $id, "--metaonly") or die "open p3 sra failed $id $!";
    while (<P>)
    {
	if (/total_bases\D+(\d+)/)
	{
	    $$uncomp += $1;
	    close(P);
	    return;
	}
    }
    close(P);
    die "couldn't find bases for $id\n";
}
