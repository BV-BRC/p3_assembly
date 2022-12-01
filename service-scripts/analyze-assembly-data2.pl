use Data::Dumper;
use strict;
use JSON::XS;
use Bio::KBase::AppService::SchedulerDB;
use Bio::P3::Workspace::WorkspaceClientExt;

my $ws = Bio::P3::Workspace::WorkspaceClientExt->new;

#
# Given results from analyze-assembly-data.pl, pull the workspace info again and
# determine assembler used.
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
	
	my $path1 = "$path/.$file/details/${file}_run_details.json";
	my $path2 = "$path/.$file/${file}_run_details.json";
	my $path3 = "$path/.$file/${file}_run_details.txt";
	my $txt;
	open(my $fh, ">", \$txt) or die $!;
	eval { $ws->copy_files_to_handles(1, undef, [[$path1, $fh]], { admin => 1 }); };
	if ($@)
	{
	    eval { $ws->copy_files_to_handles(1, undef, [[$path2, $fh]], { admin => 1 }); };
	}
	if ($@)
	{
	    eval { $ws->copy_files_to_handles(1, undef, [[$path3, $fh]], { admin => 1 }); };
	}
	
	warn "$path1: $@ " if $@;
	next unless $txt;
	
	my $dat = eval { decode_json($txt); };
	next unless $dat;
	
	next if @{$dat->{problem}};

	my $asm = $dat->{assembly}->{assembler};

	my(undef, $owner, $size, $elap, $mem, $time) = @{$data{$id}};
	
	print join("\t", $id, $owner, $size, $elap, $mem, $asm), "\n";
    }
}


