#
# The Genome Assembly application, version 2, using p3_assembly.
#

=head1 NAME

App-GenomeAssembly2 - assemble a set of reads

=head1 SYNOPSIS

    App-GenomeAssembly [--preflight] service-url app-definition parameters

=head1 DESCRIPTION

Assemble a set of reads.

=head2 PREFLIGHT INFORMATION

On a preflight request, we will generate a JSON object with the following
key/value pairs:

=over 4

=item ram

Requested memory. For standard run we request 128GB.

=item cpu

Requested CPUs.
    
=cut

use strict;
use Cwd;
use Carp;
use Data::Dumper;
use File::Temp;
use File::Basename;
use IPC::Run 'run';
use POSIX;
use File::Slurp;
use JSON::XS;

use Bio::KBase::AppService::AppScript;
use Bio::KBase::AppService::ReadSet;

my $script = Bio::KBase::AppService::AppScript->new(\&assemble, \&preflight);

my $download_path;

my $MAX_BASES = 1e10;

my $rc = $script->run(\@ARGV);

exit $rc;


sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;

    print STDERR "preflight GenomeAssembly2 ", Dumper($params, $app);

    my $token = $app->token();
    my $ws = $app->workspace();

    my $readset;
    eval {
	$readset = Bio::KBase::AppService::ReadSet->create_from_asssembly_params($params);
    };
    if ($@)
    {
	die "Error parsing assembly parameters: $@";
    }

    my($ok, $errs, $comp_size, $uncomp_size) = $readset->validate($ws);

    if (!$ok)
    {
	die "Reads as defined in parameters failed to validate. Errors:\n\t" . join("\n\t", @$errs) . "\n";
    }
    print STDERR "comp=$comp_size uncomp=$uncomp_size\n";

    my $est_comp = $comp_size + 0.75 * $uncomp_size;
    $est_comp /= 1e6;
    #
    # Estimated conservative rate is 10sec/MB for compressed data under 1.5G, 4sec/GM for data over that.
    my $est_time = int($est_comp < 1500 ? (10 * $est_comp) : (4 * $est_comp));

    my %plats = map { $_->{platform} => 1 } $readset->libraries;

    if ($plats{nanopore} || $plats{pacbio})
    {
	$est_time *= 5;
    }

    #
    # Unicycler etc is a brave new world.
    #
    $est_time *= 10;
    $est_time = 3600 if $est_time < 3600;

    # Estimated compressed storage based on input compressed size, converted at 75% compression estimate.
    my $est_storage = int(1.3e6 * $est_comp / 0.75);

    my $est_cpu = 12;

    my $est_ram = "128G";
    if ($est_time < 6*60*60)
    {
	$est_cpu = 12;
	$est_ram = "48000M";
    }

    return {
	cpu => $est_cpu,
	memory => $est_ram,
	runtime => $est_time,
	storage => $est_storage,
    };
}

sub assemble
{
    my($app, $app_def, $raw_params, $params) = @_;

    print "Begin assembly ", Dumper($app_def, $raw_params, $params);

    #
    # We need to arrange our path. Hardcode for now. Think about parameterizing.
    #
    my @path_additions = ("$ENV{KB_RUNTIME}/spades-4.0.0/bin",
			  "$ENV{KB_RUNTIME}/samtools-1.20/bin"); #,  "$ENV{KB_RUNTIME}/bowtie2-v2.2.9/bin");

    my $token = $app->token();
    my $ws = $app->workspace();

    my $cleanup = $params->{debug} > 0 ? 0 : 1;

    my $tmpdir = File::Temp->newdir( CLEANUP => $cleanup );
    print STDERR "Debug=$params->{debug_level} cleanup=$cleanup tmpdir=$tmpdir\n";
    $download_path = $tmpdir;

    my $asm_out = "$tmpdir/assembly";
    mkdir($asm_out) or die "cannot mkdir $asm_out: $!";
    my $stage_dir = "$tmpdir/staging";
    mkdir($stage_dir) or die "cannot mkdir $tmpdir/staging: $!";
    
    my $readset = Bio::KBase::AppService::ReadSet->create_from_asssembly_params($params, 1);

    my($ok, $errs) = $readset->validate($ws);

    if (!$ok)
    {
	die "Reads as defined in parameters failed to validate. Errors:\n\t" . join("\n\t", @$errs);
    }

    $readset->localize_libraries($stage_dir);

    $readset->stage_in($ws);

    my @params = $readset->build_p3_assembly_arguments();

    #
    # If we are running under Slurm, pick up our memory and CPU limits.
    #
    my $mem = $ENV{P3_ALLOCATED_MEMORY};
    my $cpu = $ENV{P3_ALLOCATED_CPU};

    if ($mem)
    {
	my $bytes;
	my %fac = (k => 1024, m => 1024*1024, g => 1024*1024*1024, t => 1024*1024*1024*1024 );
	my($val, $suffix) = $mem =~ /^(.*)([mkgt])$/i;
	if ($suffix)
	{
	    $bytes = $val * $fac{lc($suffix)};
	}
	else
	{
	    $bytes = $mem;
	}
	$mem = int($bytes / (1024*1024*1024));
	push(@params, "-m", $mem);
    }

    #
    # TESTING
    #
    # push(@params, "--kmers", 103);
    
    push(@params, "-t", $cpu) if $cpu;

    my $prefix = $params->{output_file} . "_";
    $prefix =~ tr/ /_/;

    push(@params, "--recipe", $params->{recipe});

    #
    # Polishing options
    #
    push(@params, "--racon_iterations", $params->{racon_iter});
    push(@params, "--pilon_iterations", $params->{pilon_iter});

    # 
    # Limit maximum allowed bases across all read files
    #
    if ($params->{max_bases})
    {
	push(@params, '--max_bases', $params->{max_bases});
    }

    #
    # Trimming options.
    #
    # We use TrimGalore here.
    #
    if ($params->{trim})
    {
	push(@params, '--trim');
    }

    #
    # Normalize options.
    #
    # We use BBNorm here.
    #
    if ($params->{normalize})
    {
	push(@params, '--normalize');
    }

    # We use filtlong here.
    #
    if ($params->{filtlong})
    {
	push(@params, '--filtlong');
    }

    # Target depth for down-sampling by bbnorm and/or filtlong
    #
    if ($params->{target_depth})
    {
	push(@params, '--target_depth', $params->{target_depth});
    }

    #
    # Contig filtering options
    #
    push(@params, "--min_contig_length", $params->{min_contig_len});

     push(@params, "--min_contig_coverage", $params->{min_contig_cov});

    #
    # Genome size; used as parameter for canu and flye
    # also used for filtlong downsampling to target_depth multiple of estimated genome size
    #
    push(@params, "--genome_size", $params->{genome_size});

    #
    # MAX_BASES is the upper limit on how much data gets analyzed
    #
    if (exists $params->{max_bases}) {
        push(@params, "--max_bases", $params->{max_bases});
        }
    elsif ($MAX_BASES) {
        push(@params, "--max_bases", $MAX_BASES);
    }

    #
    # Request bandage plots
    #
    push(@params, "--bandage");
    push(@params, "--maxContigsForBandage", 100); # bandage plots of assemblies with too many contigs are not helpful

    #
    # We change directory to $asm_out since p3x-assembly
    # writes its output to the current working directory.
    #
    
    my @cmd = ("p3x-assembly",
	       "--prefix", "$prefix",
	       "--path-prefix", @path_additions,
	       @params);

    print "Start assembler: @cmd\n";
    print Dumper(\@cmd, \@path_additions);

    open(ASM_OUT, ">", "$asm_out/p3x-assembly.stdout") or die "Cannot write $asm_out/p3x-assembly.stdout: $!";
    open(ASM_ERR, ">", "$asm_out/p3x-assembly.stderr") or die "Cannot write $asm_out/p3x-assembly.stderr: $!";
    my $asm_ok = run(\@cmd,
		     init => sub {
			 chdir $asm_out or die "Cannot chdir $asm_out: $!";
			 # we use --path-prefix now
#			 $ENV{PATH} = join(":", @path_additions, $ENV{PATH});
		     },
		     ">", sub { my($dat) = @_; print $dat; print ASM_OUT $dat; },
		     "2>", sub { my($dat) = @_; print STDERR $dat; print ASM_ERR $dat; },
		     # ">", "$asm_out/p3x-assembly.stdout",
		     # "2>", "$asm_out/p3x-assembly.stderr",
		    );
    my $asm_rc = $?;

    my $output_folder = $app->result_folder();

    #
    # The p3x-assembly script will generate a folder called "p3_assembly_work". This will in turn
    # contain a folder named "save" containing all output which should be copied back to the workspace.
    #

    my $type_map = {
	html => "html",
	fasta => "contigs",
	run_details => 'txt',
	txt => 'txt',
	log => 'txt',
	stdout => "txt",
	stderr => "txt",
	jpg => 'jpg',
    };

    my $save_path = "$asm_out/p3_assembly_work/save";

    if (-s "$asm_out/p3x-assembly.stdout")
    {
	$ws->save_file_to_file("$asm_out/p3x-assembly.stdout", {}, "$output_folder/p3x-assembly.stdout", 'txt', 1, 1);
    }
    if (-s "$asm_out/p3x-assembly.stderr")
    {
	$ws->save_file_to_file("$asm_out/p3x-assembly.stderr", {}, "$output_folder/p3x-assembly.stderr", 'txt', 1, 1);
    }

    if (opendir(DIR, $save_path))
    {
	while (my $f = readdir(DIR))
	{
	    next if $f =~ /^\./;
	    my $path = "$save_path/$f";
	    my($suffix) = $f =~ /\.([^.]+)$/;
	    my $type = $type_map->{$suffix} // 'txt';
	    
	    if (-f $path)
	    {
		$ws->save_file_to_file($path, {}, "$output_folder/$f", $type, 1, 1);
	    }
	    else
	    {
		$ws->upload_folder($path, "$output_folder", { type_map => $type_map });
	    }
	}
	closedir(DIR);
    }
    else
    {
	warn "Cannot opendir $save_path: $!";
    }

    if (!$asm_ok)
    {
	die "Assembler failed with rc=$asm_rc";
    }
}


