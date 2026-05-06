#!/usr/bin/perl

# based on original work by Egon A. Ozer 
# source: (https://github.com/egonozer/script_micro_wgs/blob/main/scripts/assembly_qc_filter.pl)


# MIT License

# Copyright (c) 2020 Egon A. Ozer
# Copyright (c) 2026 Hassan Ghayas

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

my $version = 0.2;

use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;
use File::Spec::Functions qw ( catfile path );
use File::Temp qw( tempdir );
use Time::HiRes qw( time );
use Getopt::Long;

my $start_time = time();

# set defaults
my $cfile;
my $read1;
my $read2;
my $min_cov     = 5;
my $min_len     = 200;
my $bwa_loc;
my $blast_loc;
my $makeblastdb_loc;
my $cont_file   = "/usr/local/share/phiX.fasta";
my $cont_id     = 98;
my $threads     = 8;
my $pref        = 'output';

my $usage = "
assembly_qc_filter.pl v$version

code based on original assembly_qc_filter.pl by:
    Egon A. Ozer (https://github.com/egonozer/script_micro_wgs/blob/main/scripts/assembly_qc_filter.pl)

1) Align reads back to assembly contigs to determine coverage stats.
2) Remove contigs that are too short, have low coverage, or align to contaminants (e.g. phiX).

Required:
  -c    assembly file of contigs in fasta format
  -1    forward read file
  -2    reverse read file
  
Optional:
  -m    Minimum average fold coverage of a contig to keep (default: $min_cov)
  -l    Minimum length of a contig to keep, in bases (default: $min_len)
  -a    path to bwa aligner executable (default: search in PATH)
  -b    path to blastn executable (default: search in PATH)
  -p    path to contaminant sequence file (default: $cont_file)
  -x    minimum percent identity for contaminant alignment (default: $cont_id)
  -t    threads (default: $threads)
  -o    output prefix (default: '$pref')
  
  Example:
  assembly_qc_filter.pl -c <assembly.fasta> -1 <read1> -2 <read2> -o <output prefix> -p <phiX.fasta/contamination.fasta>
";

Getopt::Long::Configure qw(gnu_getopt);
GetOptions(
    'c=s' => \$cfile,
    '1=s' => \$read1,
    '2=s' => \$read2,
    'm=f' => \$min_cov,
    'l=i' => \$min_len,
    'a=s' => \$bwa_loc,
    'b=s' => \$blast_loc,
    'p=s' => \$cont_file,
    'x=f' => \$cont_id,
    't=i' => \$threads,
    'o=s' => \$pref
) or die "$usage";
die "$usage" unless $cfile and $read1 and $read2;

# Setup temp directory for automatic cleanup
my $tmpdir = tempdir( CLEANUP => 1 );

# Check executables
$bwa_loc = find_exe("bwa", $bwa_loc);
$blast_loc = find_exe("blastn", $blast_loc);
my $blast_bin_dir = dirname($blast_loc);
$makeblastdb_loc = catfile($blast_bin_dir, "makeblastdb");
unless (-x $makeblastdb_loc) {
    $makeblastdb_loc = find_exe("makeblastdb");
}

## 1: Parse Contigs
print STDERR "Counting contigs...\n";
my @order;
my @pre_contig_lengs;
my ($tot_gc, $tot_non_n) = (0, 0);
my %starts;
my %stops;

open (my $cin, "<$cfile") or die "ERROR: Can't open $cfile: $!\n";
my $id;
my $current_seq = "";
while (my $line = <$cin>){
    chomp $line;
    if ($line =~ m/^>/){
        if ($id){
            process_contig_data($id, $current_seq);
        }
        $id = substr($line, 1);
        $id =~ s/\s.*$//;
        $current_seq = "";
        next;
    }
    $line =~ s/\s//g;
    $current_seq .= $line;
}
process_contig_data($id, $current_seq) if $id;
close ($cin);

sub process_contig_data {
    my ($cid, $seq) = @_;
    my $leng = length($seq);
    my $non_n = ($seq =~ tr/ACGTacgt/ACGTacgt/);
    my $gc = ($seq =~ tr/GCgc/GCgc/);
    push @order, [$cid, $leng, $non_n, $gc, $seq];
    push @pre_contig_lengs, $leng;
    @{$starts{$cid}} = (0) x ($leng + 2); # +2 to handle 1-based indexing and end logic
    @{$stops{$cid}} = (0) x ($leng + 2);
    $tot_non_n += $non_n;
    $tot_gc += $gc;
}

die "ERROR: contig file $cfile has no fasta records\n" unless @order;

# Pre-filter stats
my @sorted_pre_lengs = sort{$a <=> $b}@pre_contig_lengs;
my $pre_num = scalar @order;
my $pre_size = 0; $pre_size += $_ for @pre_contig_lengs;
my ($pre_min, $pre_max) = ($sorted_pre_lengs[0], $sorted_pre_lengs[-1]);
my ($pre_avg, $pre_stdev) = average(\@pre_contig_lengs);
my $pre_med = median(\@sorted_pre_lengs);
my $pre_n50 = n50(\@sorted_pre_lengs);
my $pre_gcpct = $tot_non_n > 0 ? sprintf("%.2f", 100 * ($tot_gc / $tot_non_n)) : "0.00";

## 2: Indexing and Alignment
print STDERR "Indexing contigs...\n";
my $bwa_idx = catfile($tmpdir, "bwa_idx");
system("$bwa_loc index -p $bwa_idx $cfile > /dev/null 2>&1") == 0 or die "ERROR: bwa index failed\n";

print STDERR "Aligning reads...\n";
open (my $bin, "$bwa_loc mem -t $threads $bwa_idx $read1 $read2 2>/dev/null | ") or die "ERROR: Can't run bwa mem: $!\n";
my ($unaligned, $readcount, $sec) = (0, 0, 0);
my @tlens;

while (my $line = <$bin>){
    next if $line =~ m/^@/;
    my @tmp = split("\t", $line);
    my ($flag, $rid, $pos, $mapq, $cigar, $tlen) = ($tmp[1], $tmp[2], $tmp[3], $tmp[4], $tmp[5], $tmp[8]);
    
    if ($flag & 256) { $sec++; next; }
    $readcount++;
    print STDERR "\rAligning ($readcount)" if $readcount % 10000 == 0;
    
    if ($flag & 4){ $unaligned++; next; }
    if ($tlen > 0){ push @tlens, $tlen; }

    # Simple CIGAR parsing for end position
    my $end = $pos;
    while ($cigar =~ /(\d+)([MDX=])/g) {
        $end += $1;
    }
    $end--; # position is inclusive
    
    if (exists $starts{$rid}) {
        $starts{$rid}[$pos]++;
        $stops{$rid}[$end]++;
    }
}
close ($bin);
print STDERR "\rAligning ($readcount) - Done\n";

my $pct_unaligned = $readcount > 0 ? sprintf("%.2f", 100 * ($unaligned / $readcount)) : "0.00";
my ($avg_ins, $stdev_ins) = average(\@tlens);
print "#total_reads\tunaligned_reads\tpct_unaligned\tavg_ins_size\tstdev_ins_size\n";
print "$readcount\t$unaligned\t$pct_unaligned\t$avg_ins\t$stdev_ins\n\n";

## 3: Calculate Coverages
print STDERR "Calculating coverages...\n";
my @results;
my @all_pre_base_covs;

foreach my $slice (@order){
    my ($id, $leng, $non_n, $gc) = @{$slice};
    my $gcpct = $non_n > 0 ? sprintf("%.2f", 100 * ($gc / $non_n)) : "0.00";
    
    my $current_cov = 0;
    my @contig_base_covs;
    for my $i (1 .. $leng){
        $current_cov += $starts{$id}[$i];
        push @contig_base_covs, $current_cov;
        push @all_pre_base_covs, $current_cov;
        $current_cov -= $stops{$id}[$i];
    }
    
    my @sorted_covs = sort{$a <=> $b}@contig_base_covs;
    my ($min, $max) = ($sorted_covs[0] // 0, $sorted_covs[-1] // 0);
    my ($avg, $stdev) = average(\@contig_base_covs);
    my $med = median(\@sorted_covs);
    
    push @results, [$id, $leng, $gcpct, $avg, $stdev, $med, $min, $max, \@sorted_covs];
}

my @sorted_all_pre_covs = sort{$a <=> $b}@all_pre_base_covs;
my ($pre_total_cov_min, $pre_total_cov_max) = ($sorted_all_pre_covs[0] // 0, $sorted_all_pre_covs[-1] // 0);
my ($pre_total_cov_avg, $pre_total_cov_stdev) = average(\@all_pre_base_covs);
my $pre_total_cov_med = median(\@sorted_all_pre_covs);

## 4: BLAST for Contaminants (Optimized)
my %hit_info;
if ($cont_file && -s $cont_file) {
    # print STDERR "Running batch BLAST against contaminants...\n";
    my $blast_db = catfile($tmpdir, "contam_db");
    system("$makeblastdb_loc -in $cont_file -dbtype nucl -out $blast_db > /dev/null 2>&1");
    
    open (my $blout, "$blast_loc -query $cfile -db $blast_db -outfmt 6 -num_threads $threads | ") or die "ERROR: Blast failed\n";
    while (my $line = <$blout>){
        chomp $line;
        my @f = split("\t", $line);
        my ($qid, $sid, $pident, $qstart, $qend) = ($f[0], $f[1], $f[2], $f[6], $f[7]);
        if ($pident >= $cont_id) {
            $hit_info{$qid} .= "$sid($qstart-$qend,$pident%) ";
        }
    }
    close $blout;
}

## 5: Filtering and Output
print STDERR "Filtering and outputting...\n";
open (my $outf, ">$pref.filtered.fasta");
open (my $outc, ">$pref.contig_stats.txt");
print $outc "id\tlength\tgc_pct\tavg_cov\tstdev_cov\tmed_cov\tmin_cov\tmax_cov\tfiltered\n";

my (@post_lengs, @post_base_covs, @filt_lengs, @filt_base_covs);
my ($post_gc, $post_non_n, $filt_gc, $filt_non_n) = (0, 0, 0, 0);
my @contaminated_ids;

for my $i (0 .. $#results) {
    my ($id, $leng, $gcpct, $avg, $stdev, $med, $min, $max, $base_covs_ref) = @{$results[$i]};
    my ($x1, $x2, $non_n, $gc, $seq) = @{$order[$i]};
    
    if ($leng >= $min_len && $avg >= $min_cov) {
        if ($hit_info{$id}) {
            my $info = $hit_info{$id}; $info =~ s/\s$//;
            print $outc "$id\t$leng\t$gcpct\t$avg\t$stdev\t$med\t$min\t$max\tY ($info)\n";
            push @contaminated_ids, $id;
        } else {
            print $outf ">$id\n$seq\n";
            print $outc "$id\t$leng\t$gcpct\t$avg\t$stdev\t$med\t$min\t$max\tN\n";
            push @post_lengs, $leng;
            push @post_base_covs, @$base_covs_ref;
            $post_gc += $gc; $post_non_n += $non_n;
        }
    } else {
        print $outc "$id\t$leng\t$gcpct\t$avg\t$stdev\t$med\t$min\t$max\tY\n";
        push @filt_lengs, $leng;
        push @filt_base_covs, @$base_covs_ref;
        $filt_gc += $gc; $filt_non_n += $non_n;
    }
}
close $outf; close $outc;

# Summary Calcs
my @sorted_post_lengs = sort{$a <=> $b}@post_lengs;
my $post_num = scalar @post_lengs;
my $post_size = 0; $post_size += $_ for @post_lengs;
my ($post_min, $post_max) = ($sorted_post_lengs[0] // 0, $sorted_post_lengs[-1] // 0);
my ($post_avg, $post_stdev) = average(\@post_lengs);
my $post_med = median(\@sorted_post_lengs);
my $post_n50 = n50(\@sorted_post_lengs);
my $post_gcpct = $post_non_n > 0 ? sprintf("%.2f", 100 * ($post_gc / $post_non_n)) : "0.00";

my @sorted_post_covs = sort{$a <=> $b}@post_base_covs;
my ($post_total_cov_min, $post_total_cov_max) = ($sorted_post_covs[0] // 0, $sorted_post_covs[-1] // 0);
my ($post_total_cov_avg, $post_total_cov_stdev) = average(\@post_base_covs);
my $post_total_cov_med = median(\@sorted_post_covs);

# Filt Summary
my @sorted_filt_lengs = sort{$a <=> $b}@filt_lengs;
my $filt_num = scalar @filt_lengs;
my $filt_size = 0; $filt_size += $_ for @filt_lengs;
my ($filt_min, $filt_max) = ($sorted_filt_lengs[0] // 0, $sorted_filt_lengs[-1] // 0);
my ($filt_avg, $filt_stdev) = average(\@filt_lengs);
my $filt_med = median(\@sorted_filt_lengs);
my $filt_n50 = n50(\@sorted_filt_lengs);
my $filt_gcpct = $filt_non_n > 0 ? sprintf("%.2f", 100 * ($filt_gc / $filt_non_n)) : "0.00";

my @sorted_filt_covs = sort{$a <=> $b}@filt_base_covs;
my ($filt_total_cov_min, $filt_total_cov_max) = ($sorted_filt_covs[0] // 0, $sorted_filt_covs[-1] // 0);
my ($filt_total_cov_avg, $filt_total_cov_stdev) = average(\@filt_base_covs);
my $filt_total_cov_med = median(\@sorted_filt_covs);

# Final Prints
print "#pre/post\tnum\tsize\tavg\tstdev\tmedian\tmin\tmax\tn50\tgc\tavg_cov\tstdev_cov\tmedian_cov\tmin_cov\tmax_cov\n";
printf "pre\t%d\t%d\t%.4f\t%.4f\t%s\t%d\t%d\t%d\t%s\t%.4f\t%.4f\t%s\t%d\t%d\n", 
    $pre_num, $pre_size, $pre_avg, $pre_stdev, $pre_med, $pre_min, $pre_max, $pre_n50, $pre_gcpct, $pre_total_cov_avg, $pre_total_cov_stdev, $pre_total_cov_med, $pre_total_cov_min, $pre_total_cov_max;
printf "post\t%d\t%d\t%.4f\t%.4f\t%s\t%d\t%d\t%d\t%s\t%.4f\t%.4f\t%s\t%d\t%d\n", 
    $post_num, $post_size, $post_avg, $post_stdev, $post_med, $post_min, $post_max, $post_n50, $post_gcpct, $post_total_cov_avg, $post_total_cov_stdev, $post_total_cov_med, $post_total_cov_min, $post_total_cov_max;
printf "filt\t%d\t%d\t%.4f\t%.4f\t%s\t%d\t%d\t%d\t%s\t%.4f\t%.4f\t%s\t%d\t%d\n", 
    $filt_num, $filt_size, $filt_avg, $filt_stdev, $filt_med, $filt_min, $filt_max, $filt_n50, $filt_gcpct, $filt_total_cov_avg, $filt_total_cov_stdev, $filt_total_cov_med, $filt_total_cov_min, $filt_total_cov_max;

print "Contaminated contigs removed: ", scalar @contaminated_ids, "\n";

open(my $outs, ">$pref.assembly_stats.tsv");
print $outs "seqID\tcontigs\tsize\tlargest_contig\tn50\tGC\tavg_cov\n";
print $outs "$pref\t$post_num\t$post_size\t$post_max\t$post_n50\t$post_gcpct\t$post_total_cov_avg\n";
close($outs);

my $end_time = time();
printf STDERR "\nTotal Execution Time: %.2f seconds\n", ($end_time - $start_time);

# Helper Subs
sub find_exe {
    my ($exe, $provided) = @_;
    if ($provided) {
        my $p = abs_path($provided);
        die "ERROR: $exe not found at $p\n" unless -x $p;
        return $p;
    }
    for my $dir (path()) {
        my $f = catfile($dir, $exe);
        return $f if -x $f;
        if ($^O eq 'MSWin32') {
            for my $ext (qw(.exe .bat .cmd)) {
                my $fe = "$f$ext";
                return $fe if -x $fe;
            }
        }
    }
    die "ERROR: $exe not found in PATH. Please specify with flags.\n";
}

sub average {
    my $pts = shift;
    return ("0.0000", "0.0000") unless @$pts;
    my $sum = 0; $sum += $_ for @$pts;
    my $avg = $sum / @$pts;
    my $sq_sum = 0; $sq_sum += ($_ - $avg)**2 for @$pts;
    my $stdev = sqrt($sq_sum / @$pts);
    return (sprintf("%.4f", $avg), sprintf("%.4f", $stdev));
}

sub median {
    my $pts = shift;
    my $count = scalar @$pts;
    return "NA" unless $count;
    # Array must be sorted
    if ($count % 2) { return $pts->[int($count/2)]; }
    return ($pts->[($count/2)-1] + $pts->[$count/2]) / 2;
}

sub n50 {
    my $pts = shift;
    my $total = 0; $total += $_ for @$pts;
    my $running = 0;
    # Array must be sorted ASC for this logic with reverse traversal or DESC
    for (reverse @$pts) {
        $running += $_;
        return $_ if $running >= $total / 2;
    }
    return 0;
}
