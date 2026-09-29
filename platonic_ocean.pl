#!/usr/bin/perl

use strict;
use warnings;
use feature qw(state);

use Cwd qw(abs_path);
use File::Find qw(find);
use File::Path qw(make_path remove_tree);
use Getopt::Long qw(GetOptions);
use List::Util qw(min max);
use POSIX qw(acos asin ceil floor);
use Time::HiRes qw(time);

use constant PI       => 4 * atan2(1, 1);
use constant EARTH_KM => 6371.0088;

# Output raster.  Unlike an ordinary 2:1 equirectangular image this has a real
# row for each pole: y=0 is +90 degrees, y=1800 is -90 degrees.
use constant XW => 3600;
use constant XH => 1801;

# Ocean-vs-land classification raster.  We build this ourselves directly from
# OSM's coastline-derived land polygons.  0.1 degree cells are deliberately
# coarse: the optimizer's winning points should be hundreds of km offshore, and
# the separate coastline-distance calculation catches small islands anyway.
use constant MASK_W => 3600;
use constant MASK_H => 1800;

# Coastline line segments are sampled along great circles at this interval for
# nearest-coast distance.  Maximum sampling-only error is about half a step.
use constant COAST_STEP_DEG => 0.05;

# 3-D spatial hash cell width in unit-sphere chord coordinates.
use constant COAST_CELL => 0.04;
use constant COAST_NCELL => 51;   # enough for [-1,+1] at width 0.04
use constant CHECKPOINT_SECONDS => 3600;

my $cache_dir = 'osm-cache';
my $out_dir   = '.';
my $trials    = 100_000;
my $keep      = 12;
my $seed      = 20260921;
my $refresh   = 0;
my $help      = 0;

GetOptions(
    'cache=s'   => \$cache_dir,
    'out=s'     => \$out_dir,
    'trials=i'  => \$trials,
    'keep=i'    => \$keep,
    'seed=i'    => \$seed,
    'refresh!'  => \$refresh,
    'help|h'    => \$help,
) or usage(2);

usage(0) if $help;
die "--trials must be >= 1\n" if $trials < 1;
die "--keep must be >= 1\n"   if $keep < 1;

make_path($cache_dir) unless -d $cache_dir;
make_path($out_dir)   unless -d $out_dir;

check_tools();


my $LAND_URL  = 'https://osmdata.openstreetmap.de/download/land-polygons-split-4326.zip';
my $COAST_URL = 'https://osmdata.openstreetmap.de/download/coastlines-split-4326.zip';

my $land_zip  = "$cache_dir/land-polygons-split-4326.zip";
my $coast_zip = "$cache_dir/coastlines-split-4326.zip";
my $land_dir  = "$cache_dir/land-polygons-split-4326";
my $coast_dir = "$cache_dir/coastlines-split-4326";

my $mask_raw  = "$cache_dir/ocean-mask-" . MASK_W . "x" . MASK_H . ".raw";
my $coast_raw = "$cache_dir/coast-samples-f32.raw";
my $base_raw  = "$cache_dir/base-xpm-" . XW . "x" . XH . ".raw";

if ($refresh) {
    unlink $land_zip  if -e $land_zip;
    unlink $coast_zip if -e $coast_zip;
    remove_tree($land_dir)  if -d $land_dir;
    remove_tree($coast_dir) if -d $coast_dir;
    unlink $mask_raw  if -e $mask_raw;
    unlink $coast_raw if -e $coast_raw;
    unlink $base_raw  if -e $base_raw;
}

fetch_if_missing($LAND_URL,  $land_zip);
fetch_if_missing($COAST_URL, $coast_zip);
unpack_if_missing($land_zip,  $land_dir);
unpack_if_missing($coast_zip, $coast_dir);

my $land_shp  = find_first_shp($land_dir);
my $coast_shp = find_first_shp($coast_dir);

print "land polygons: $land_shp\n";
print "coastlines:    $coast_shp\n";

build_ocean_mask_from_shp($land_shp, $mask_raw)
    if $refresh || !-s $mask_raw;
my ($mask_w, $mask_h, $ocean_mask) = read_raw_mask($mask_raw);
die "unexpected ocean mask size ${mask_w}x${mask_h}\n"
    unless $mask_w == MASK_W && $mask_h == MASK_H;

build_coast_cache_from_shp($coast_shp, $coast_raw, $base_raw)
    if $refresh || !-s $coast_raw || !-s $base_raw;

my $base_xpm = read_exact_file($base_raw, XW * XH);
my %coast_grid;
my $coast_samples = load_coast_sample_cache($coast_raw, \%coast_grid);
print "coast samples: $coast_samples\n";

my %solids = platonic_solids();
my @order = qw(tetrahedron cube octahedron dodecahedron icosahedron);

for my $solid_index (0 .. $#order) {
    my $name = $order[$solid_index];
    my $verts = $solids{$name};
    normalize_vertices($verts);
    my $edges = find_edges($verts);

    my %expect = (
        tetrahedron  => 6,
        cube         => 12,
        octahedron   => 12,
        dodecahedron => 30,
        icosahedron  => 30,
    );
    die "$name edge detection got " . scalar(@$edges) . ", expected $expect{$name}\n"
        unless @$edges == $expect{$name};

    my $csv = "$out_dir/$name.csv";
    my $xpm = "$out_dir/$name.xpm";

    if (-s $csv && -s $xpm) {
        print "\n$name: already complete; keeping $csv and $xpm\n";
        next;
    }

    print "\n$name: " . scalar(@$verts) . " vertices, " . scalar(@$edges) . " edges\n";
    my ($q, $score) = optimize_orientation($name, $verts, $solid_index);
    my $rverts = rotate_vertices($verts, $q);

    write_csv($csv, $name, $q, $score, $rverts);
    write_xpm($xpm, $name, $base_xpm, $rverts, $edges);
    clear_checkpoint($name);

    printf "%s maximin: %.3f km\n", $name, $score;
    print "  $csv\n  $xpm\n";
}

print "\nDone: 10 result files written to $out_dir\n";
exit 0;

sub usage {
    my ($status) = @_;
    print <<'USAGE';
usage: platonic_ocean.pl [options]

Find maximin-ocean orientations for all five Platonic solids using current
OpenStreetMap coastline data.  Produces exactly these result files:

  tetrahedron.csv      tetrahedron.xpm
  cube.csv             cube.xpm
  octahedron.csv       octahedron.xpm
  dodecahedron.csv     dodecahedron.xpm
  icosahedron.csv      icosahedron.xpm

Options:
  --trials N       random orientations per solid (default 100000)
  --keep N         top random candidates locally refined (default 12)
  --seed N         deterministic PRNG seed (default 20260921)
  --cache DIR      download/derived-data cache (default osm-cache)
  --out DIR        output directory (default .)
  --refresh        re-download/rebuild cached OSM data
  --help

Search checkpoints are written once per hour as .<solid>.checkpoint in --out.
On restart, an unfinished solid resumes from its saved random-search trial.
Completed solids with both CSV and XPM outputs are skipped.

External commands required:
  curl or wget, unzip

The OSM downloads are:
  land-polygons-split-4326.zip    coastline-derived continents/islands
  coastlines-split-4326.zip       coastline lines

The Shapefiles are parsed directly in Perl; GDAL is not required.
USAGE
    exit $status;
}

sub command_exists {
    my ($cmd) = @_;
    for my $d (split /:/, ($ENV{PATH} // '')) {
        return 1 if -x "$d/$cmd";
    }
    return 0;
}

sub check_tools {
    my @missing;
    push @missing, 'unzip' unless command_exists('unzip');
    push @missing, 'curl-or-wget'
        unless command_exists('curl') || command_exists('wget');
    if (@missing) {
        die "missing required command(s): @missing\n";
    }
}

sub run_cmd {
    my (@cmd) = @_;
    print '+ ', join(' ', map { shell_display($_) } @cmd), "\n";
    system(@cmd) == 0 or die "command failed ($?): @cmd\n";
}

sub shell_display {
    my ($s) = @_;
    return $s if $s =~ m{^[A-Za-z0-9_./:=+,-]+$};
    $s =~ s/'/'\\''/g;
    return "'$s'";
}

sub fetch_if_missing {
    my ($url, $path) = @_;
    return if -s $path;

    my $tmp = "$path.part";
    unlink $tmp if -e $tmp;

    if (command_exists('curl')) {
        run_cmd('curl', '-L', '--fail', '--retry', '3', '-o', $tmp, $url);
    }
    else {
        run_cmd('wget', '-O', $tmp, $url);
    }

    rename($tmp, $path) or die "rename $tmp -> $path: $!\n";
}

sub unpack_if_missing {
    my ($zip, $dir) = @_;
    return if -d $dir && find_first_shp_or_undef($dir);
    remove_tree($dir) if -d $dir;
    make_path($dir);
    run_cmd('unzip', '-q', $zip, '-d', $dir);
}

sub find_first_shp_or_undef {
    my ($dir) = @_;
    my $found;
    find({
        wanted => sub {
            return if $found;
            $found = $File::Find::name if -f $_ && /\.shp\z/i;
        },
        no_chdir => 1,
    }, $dir);
    return $found;
}

sub find_first_shp {
    my ($dir) = @_;
    my $p = find_first_shp_or_undef($dir);
    die "no .shp found under $dir\n" unless $p;
    return $p;
}

sub read_exact {
    my ($fh, $n, $what) = @_;
    my $buf = '';
    while (length($buf) < $n) {
        my $got = read($fh, my $chunk, $n - length($buf));
        die "read $what: $!\n" unless defined $got;
        die "unexpected EOF reading $what\n" if $got == 0;
        $buf .= $chunk;
    }
    return $buf;
}

sub read_exact_file {
    my ($path, $bytes) = @_;
    open my $fh, '<:raw', $path or die "open $path: $!\n";
    my $buf = read_exact($fh, $bytes, $path);
    my $extra = read($fh, my $junk, 1);
    die "$path: file is larger than expected\n" if defined($extra) && $extra;
    close $fh;
    return $buf;
}

sub walk_shapefile {
    my ($path, $callback) = @_;
    open my $fh, '<:raw', $path or die "open $path: $!\n";

    my $hdr = read_exact($fh, 100, "$path header");
    my $code = unpack('N', substr($hdr, 0, 4));
    die "$path: bad Shapefile magic $code\n" unless $code == 9994;

    my $file_type = unpack('V', substr($hdr, 32, 4));
    print "  shapefile type $file_type: $path\n";

    my $records = 0;
    while (1) {
        my $got = read($fh, my $rh, 8);
        die "read $path record header: $!\n" unless defined $got;
        last if $got == 0;
        die "$path: short record header\n" unless $got == 8;

        my ($record_no, $words) = unpack('N2', $rh);
        my $bytes = $words * 2;
        my $rec = read_exact($fh, $bytes, "$path record $record_no");
        $records++;

        my $type = unpack('V', substr($rec, 0, 4));
        next if $type == 0; # Null shape

        # PolyLine / Polygon, including Z and M variants.  Their XY prefix is
        # identical; any Z/M arrays occur after the points and can be ignored.
        unless ($type == 3 || $type == 5 ||
                $type == 13 || $type == 15 ||
                $type == 23 || $type == 25) {
            die "$path record $record_no: unsupported shape type $type\n";
        }

        die "$path record $record_no: short shape record\n" if length($rec) < 44;
        my @bbox = unpack('d<4', substr($rec, 4, 32));
        my ($nparts, $npoints) = unpack('V2', substr($rec, 36, 8));

        my $parts_bytes = 4 * $nparts;
        my $points_off = 44 + $parts_bytes;
        my $points_bytes = 16 * $npoints;
        die "$path record $record_no: malformed part/point counts\n"
            if $points_off + $points_bytes > length($rec);

        my @parts = $nparts
            ? unpack("V$nparts", substr($rec, 44, $parts_bytes))
            : ();
        my @xy = $npoints
            ? unpack("d<" . (2 * $npoints), substr($rec, $points_off, $points_bytes))
            : ();

        $callback->({
            record_no => $record_no,
            type      => $type,
            bbox      => \@bbox,       # xmin,ymin,xmax,ymax
            parts     => \@parts,
            xy        => \@xy,         # x0,y0,x1,y1,...
            npoints   => $npoints,
        });
    }

    close $fh;
    return $records;
}

sub shape_part_range {
    my ($shape, $part_index) = @_;
    my $start = $shape->{parts}[$part_index];
    my $end = ($part_index + 1 < @{$shape->{parts}})
        ? $shape->{parts}[$part_index + 1]
        : $shape->{npoints};
    return ($start, $end);
}

sub mask_scan_lat {
    my ($y) = @_;
    return 90.0 - ($y + 0.5) * 180.0 / MASK_H;
}

sub lon_to_mask_edge_x {
    my ($lon) = @_;
    return ($lon + 180.0) * MASK_W / 360.0 - 0.5;
}

sub rasterize_land_shape {
    my ($shape, $maskref) = @_;
    my %hits;

    for my $pi (0 .. $#{$shape->{parts}}) {
        my ($start, $end) = shape_part_range($shape, $pi);
        next if $end - $start < 3;

        for my $i ($start .. $end - 2) {
            my $j = $i + 1;
            my ($x1,$y1) = @{$shape->{xy}}[2*$i, 2*$i+1];
            my ($x2,$y2) = @{$shape->{xy}}[2*$j, 2*$j+1];
            next if $y1 == $y2;

            # The split WGS84 land polygons normally do not cross the date
            # line inside a record.  Handle an accidental short crossing anyway.
            if (abs($x2 - $x1) > 180.0) {
                $x2 += 360.0 if $x2 < $x1;
                $x1 += 360.0 if $x1 < $x2 - 180.0;
            }

            my $maxlat = max($y1,$y2);
            my $minlat = min($y1,$y2);
            my $ya = floor((90.0 - $maxlat) * MASK_H / 180.0) - 1;
            my $yb = ceil ((90.0 - $minlat) * MASK_H / 180.0) + 1;
            $ya = 0 if $ya < 0;
            $yb = MASK_H - 1 if $yb >= MASK_H;

            for my $row ($ya .. $yb) {
                my $lat = mask_scan_lat($row);

                # Standard half-open scanline crossing rule avoids counting
                # shared polygon vertices twice.
                next unless (($y1 <= $lat && $lat < $y2) ||
                             ($y2 <= $lat && $lat < $y1));

                my $t = ($lat - $y1) / ($y2 - $y1);
                my $lon = $x1 + $t * ($x2 - $x1);
                $lon -= 360.0 while $lon >= 180.0;
                $lon += 360.0 while $lon < -180.0;
                push @{$hits{$row}}, $lon;
            }
        }
    }

    for my $row (keys %hits) {
        my @xs = sort { $a <=> $b } @{$hits{$row}};
        next unless @xs >= 2;

        # Even/odd fill across all rings in this Polygon record.
        for (my $i = 0; $i + 1 < @xs; $i += 2) {
            my ($a,$b) = ($xs[$i], $xs[$i+1]);
            ($a,$b) = ($b,$a) if $a > $b;

            my $xa = ceil(lon_to_mask_edge_x($a));
            my $xb = floor(lon_to_mask_edge_x($b));
            $xa = 0 if $xa < 0;
            $xb = MASK_W - 1 if $xb >= MASK_W;
            next if $xb < $xa;

            substr($$maskref, $row * MASK_W + $xa, $xb - $xa + 1,
                   "\x00" x ($xb - $xa + 1));
        }
    }
}

sub build_ocean_mask_from_shp {
    my ($shp, $raw) = @_;
    print "building ocean/land mask directly from Shapefile...\n";

    # 1 = ocean, 0 = land.
    my $mask = "\x01" x (MASK_W * MASK_H);
    my $seen = 0;

    walk_shapefile($shp, sub {
        my ($shape) = @_;
        my $base_type = $shape->{type} % 10;
        die "$shp: expected Polygon shapes, got type $shape->{type}\n"
            unless $base_type == 5;
        rasterize_land_shape($shape, \$mask);
        $seen++;
        print "  land polygon records: $seen\r" if ($seen % 1000) == 0;
    });
    print "  land polygon records: $seen\n";

    open my $fh, '>:raw', $raw or die "open $raw: $!\n";
    print {$fh} $mask or die "write $raw: $!\n";
    close $fh or die "close $raw: $!\n";
}

sub read_raw_mask {
    my ($path) = @_;
    my $data = read_exact_file($path, MASK_W * MASK_H);
    return (MASK_W, MASK_H, $data);
}

sub append_coast_sample {
    my ($fh, $x,$y,$z) = @_;
    print {$fh} pack('f<3', $x,$y,$z) or die "write coast sample cache: $!\n";
}

sub process_coast_part {
    my ($shape, $pi, $sample_fh, $countref, $pixref) = @_;
    my ($start, $end) = shape_part_range($shape, $pi);
    return if $end <= $start;

    for my $i ($start .. $end - 2) {
        my $j = $i + 1;
        my $a = [$shape->{xy}[2*$i], $shape->{xy}[2*$i+1]];
        my $b = [$shape->{xy}[2*$j], $shape->{xy}[2*$j+1]];

        slerp_points($a, $b, COAST_STEP_DEG, sub {
            my ($x,$y,$z) = @_;
            append_coast_sample($sample_fh, $x,$y,$z);
            $$countref++;
            my ($lon,$lat) = vec_to_lonlat($x,$y,$z);
            set_geo_pixel($pixref, $lon,$lat, 'c');
        });
    }

    my $last = $end - 1;
    my $lon = $shape->{xy}[2*$last];
    my $lat = $shape->{xy}[2*$last+1];
    my ($x,$y,$z) = lonlat_to_vec($lon,$lat);
    append_coast_sample($sample_fh, $x,$y,$z);
    $$countref++;
    set_geo_pixel($pixref, $lon,$lat,'c');
}

sub build_coast_cache_from_shp {
    my ($shp, $samples_path, $base_path) = @_;
    print "building coastline cache directly from Shapefile...\n";

    my $pix = ' ' x (XW * XH);
    draw_grid(\$pix);

    my $tmp_samples = "$samples_path.part";
    unlink $tmp_samples if -e $tmp_samples;
    open my $sfh, '>:raw', $tmp_samples or die "open $tmp_samples: $!\n";

    my $samples = 0;
    my $records = 0;
    walk_shapefile($shp, sub {
        my ($shape) = @_;
        my $base_type = $shape->{type} % 10;
        die "$shp: expected PolyLine shapes, got type $shape->{type}\n"
            unless $base_type == 3;

        for my $pi (0 .. $#{$shape->{parts}}) {
            process_coast_part($shape, $pi, $sfh, \$samples, \$pix);
        }

        $records++;
        print "  coastline records: $records; samples: $samples\r"
            if ($records % 5000) == 0;
    });

    close $sfh or die "close $tmp_samples: $!\n";
    rename($tmp_samples, $samples_path)
        or die "rename $tmp_samples -> $samples_path: $!\n";

    open my $bfh, '>:raw', $base_path or die "open $base_path: $!\n";
    print {$bfh} $pix or die "write $base_path: $!\n";
    close $bfh or die "close $base_path: $!\n";

    print "  coastline records: $records; samples: $samples\n";
}

sub add_coast_sample_packed {
    my ($grid,$x,$y,$z) = @_;
    my $ix = coast_cell_coord($x);
    my $iy = coast_cell_coord($y);
    my $iz = coast_cell_coord($z);
    my $key = coast_key($ix,$iy,$iz);
    $grid->{$key} .= pack('f<3', $x,$y,$z);
}

sub load_coast_sample_cache {
    my ($path, $grid) = @_;
    open my $fh, '<:raw', $path or die "open $path: $!\n";
    my $count = 0;

    while (1) {
        my $got = read($fh, my $buf, 12 * 65536);
        die "read $path: $!\n" unless defined $got;
        last if $got == 0;
        die "$path: truncated float record\n" if $got % 12;

        my @f = unpack('f<*', $buf);
        for (my $i=0; $i<@f; $i+=3) {
            add_coast_sample_packed($grid, @f[$i,$i+1,$i+2]);
            $count++;
        }
    }

    close $fh;
    return $count;
}

sub clamp {
    my ($v, $lo, $hi) = @_;
    return $lo if $v < $lo;
    return $hi if $v > $hi;
    return $v;
}

sub lonlat_to_vec {
    my ($lon_deg, $lat_deg) = @_;
    my $lon = $lon_deg * PI / 180.0;
    my $lat = $lat_deg * PI / 180.0;
    my $cl = cos($lat);
    return (
        $cl * sin($lon),
        sin($lat),
        $cl * cos($lon),
    );
}

sub vec_to_lonlat {
    my ($x, $y, $z) = @_;
    my $lat = asin(clamp($y, -1, 1)) * 180.0 / PI;
    my $lon = atan2($x, $z) * 180.0 / PI;
    $lon = -180.0 if $lon >= 180.0;
    return ($lon, $lat);
}

sub normalize3 {
    my ($x, $y, $z) = @_;
    my $n = sqrt($x*$x + $y*$y + $z*$z);
    die "zero vector\n" if $n == 0;
    return ($x/$n, $y/$n, $z/$n);
}

sub slerp_points {
    my ($a, $b, $max_step_deg, $callback) = @_;
    my ($ax,$ay,$az) = lonlat_to_vec($a->[0], $a->[1]);
    my ($bx,$by,$bz) = lonlat_to_vec($b->[0], $b->[1]);
    my $dot = clamp($ax*$bx + $ay*$by + $az*$bz, -1, 1);
    my $ang = acos($dot);

    if ($ang < 1e-12) {
        $callback->($ax,$ay,$az);
        return;
    }

    my $steps = ceil(($ang * 180.0 / PI) / $max_step_deg);
    $steps = 1 if $steps < 1;
    my $sa = sin($ang);

    for my $i (0 .. $steps - 1) {
        my $t = $i / $steps;
        my $wa = sin((1-$t)*$ang) / $sa;
        my $wb = sin($t*$ang) / $sa;
        my ($x,$y,$z) = normalize3(
            $wa*$ax + $wb*$bx,
            $wa*$ay + $wb*$by,
            $wa*$az + $wb*$bz,
        );
        $callback->($x,$y,$z);
    }
}

sub coast_cell_coord {
    my ($v) = @_;
    my $i = floor(($v + 1.0) / COAST_CELL);
    $i = 0 if $i < 0;
    $i = COAST_NCELL - 1 if $i >= COAST_NCELL;
    return $i;
}

sub coast_key {
    my ($ix,$iy,$iz) = @_;
    return $ix + COAST_NCELL * ($iy + COAST_NCELL * $iz);
}

sub nearest_coast_km {
    my ($x,$y,$z) = @_;
    my $cx = coast_cell_coord($x);
    my $cy = coast_cell_coord($y);
    my $cz = coast_cell_coord($z);
    my $best2 = 4.0;

    for my $r (0 .. COAST_NCELL) {
        my $xmin = max(0, $cx-$r); my $xmax = min(COAST_NCELL-1, $cx+$r);
        my $ymin = max(0, $cy-$r); my $ymax = min(COAST_NCELL-1, $cy+$r);
        my $zmin = max(0, $cz-$r); my $zmax = min(COAST_NCELL-1, $cz+$r);

        for my $ix ($xmin .. $xmax) {
            for my $iy ($ymin .. $ymax) {
                for my $iz ($zmin .. $zmax) {
                    next if $r > 0 &&
                        abs($ix-$cx) < $r && abs($iy-$cy) < $r && abs($iz-$cz) < $r;
                    my $bucket = $coast_grid{coast_key($ix,$iy,$iz)} or next;
                    for (my $off=0; $off<length($bucket); $off+=12) {
                        my ($bx,$by,$bz) = unpack('f<3', substr($bucket,$off,12));
                        my $dx = $x - $bx;
                        my $dy = $y - $by;
                        my $dz = $z - $bz;
                        my $d2 = $dx*$dx + $dy*$dy + $dz*$dz;
                        $best2 = $d2 if $d2 < $best2;
                    }
                }
            }
        }

        if ($best2 < 4.0) {
            my $need = ceil(sqrt($best2) / COAST_CELL) + 2;
            if ($r >= $need) {
                my $chord = sqrt($best2);
                my $ang = 2.0 * asin(clamp($chord / 2.0, 0, 1));
                return EARTH_KM * $ang;
            }
        }
    }

    die "nearest-coast search failed\n";
}

sub is_ocean_vec {
    my ($x,$y,$z) = @_;
    my ($lon,$lat) = vec_to_lonlat($x,$y,$z);
    my $px = int(($lon + 180.0) / 360.0 * $mask_w);
    my $py = int((90.0 - $lat) / 180.0 * $mask_h);
    $px = 0 if $px < 0; $px = $mask_w-1 if $px >= $mask_w;
    $py = 0 if $py < 0; $py = $mask_h-1 if $py >= $mask_h;
    return vec($ocean_mask, $py*$mask_w + $px, 8) != 0;
}

sub platonic_solids {
    my $phi = (1 + sqrt(5)) / 2;
    my $ip  = 1 / $phi;

    my @tetra = (
        [ 1, 1, 1], [ 1,-1,-1], [-1, 1,-1], [-1,-1, 1],
    );

    my @cube;
    for my $x (-1,1) { for my $y (-1,1) { for my $z (-1,1) {
        push @cube, [$x,$y,$z];
    }}}

    my @octa = (
        [ 1,0,0], [-1,0,0], [0, 1,0], [0,-1,0], [0,0, 1], [0,0,-1],
    );

    my @icosa = (
        [-1, $phi,0], [ 1, $phi,0], [-1,-$phi,0], [ 1,-$phi,0],
        [0,-1, $phi], [0, 1, $phi], [0,-1,-$phi], [0, 1,-$phi],
        [$phi,0,-1], [$phi,0, 1], [-$phi,0,-1], [-$phi,0, 1],
    );

    my @dodeca;
    for my $x (-1,1) { for my $y (-1,1) { for my $z (-1,1) {
        push @dodeca, [$x,$y,$z];
    }}}
    for my $a (-1,1) { for my $b (-1,1) {
        push @dodeca, [0, $a*$ip, $b*$phi];
        push @dodeca, [$a*$ip, $b*$phi, 0];
        push @dodeca, [$a*$phi, 0, $b*$ip];
    }}

    return (
        tetrahedron  => \@tetra,
        cube         => \@cube,
        octahedron   => \@octa,
        dodecahedron => \@dodeca,
        icosahedron  => \@icosa,
    );
}

sub normalize_vertices {
    my ($verts) = @_;
    for my $v (@$verts) {
        @$v = normalize3(@$v);
    }
}

sub find_edges {
    my ($verts) = @_;
    my $mind2 = 1e99;
    for my $i (0 .. $#$verts-1) {
        for my $j ($i+1 .. $#$verts) {
            my $d2 = 0;
            for my $k (0..2) {
                my $d = $verts->[$i][$k] - $verts->[$j][$k];
                $d2 += $d*$d;
            }
            $mind2 = $d2 if $d2 > 1e-12 && $d2 < $mind2;
        }
    }

    my @edges;
    my $tol = $mind2 * 1e-7 + 1e-10;
    for my $i (0 .. $#$verts-1) {
        for my $j ($i+1 .. $#$verts) {
            my $d2 = 0;
            for my $k (0..2) {
                my $d = $verts->[$i][$k] - $verts->[$j][$k];
                $d2 += $d*$d;
            }
            push @edges, [$i,$j] if abs($d2-$mind2) <= $tol;
        }
    }
    return \@edges;
}

sub q_normalize {
    my ($q) = @_;
    my $n = sqrt($q->[0]**2 + $q->[1]**2 + $q->[2]**2 + $q->[3]**2);
    return [$q->[0]/$n, $q->[1]/$n, $q->[2]/$n, $q->[3]/$n];
}

sub q_mul {
    my ($a,$b) = @_;
    my ($aw,$ax,$ay,$az) = @$a;
    my ($bw,$bx,$by,$bz) = @$b;
    return [
        $aw*$bw - $ax*$bx - $ay*$by - $az*$bz,
        $aw*$bx + $ax*$bw + $ay*$bz - $az*$by,
        $aw*$by - $ax*$bz + $ay*$bw + $az*$bx,
        $aw*$bz + $ax*$by - $ay*$bx + $az*$bw,
    ];
}

sub q_axis_angle {
    my ($ax,$ay,$az,$deg) = @_;
    ($ax,$ay,$az) = normalize3($ax,$ay,$az);
    my $h = $deg * PI / 360.0;
    my $s = sin($h);
    return [cos($h), $ax*$s, $ay*$s, $az*$s];
}

sub q_random {
    my $u1 = rand();
    my $u2 = rand();
    my $u3 = rand();
    my $a = 2*PI*$u2;
    my $b = 2*PI*$u3;

    # Shoemake uniform random unit quaternion.  Stored as [w,x,y,z].
    my $x = sqrt(1-$u1) * sin($a);
    my $y = sqrt(1-$u1) * cos($a);
    my $z = sqrt($u1)   * sin($b);
    my $w = sqrt($u1)   * cos($b);
    return [$w,$x,$y,$z];
}

sub rotate_one {
    my ($q,$v) = @_;
    my ($w,$qx,$qy,$qz) = @$q;
    my ($x,$y,$z) = @$v;

    # v' = v + w*(2 qv x v) + qv x (2 qv x v)
    my $tx = 2 * ($qy*$z - $qz*$y);
    my $ty = 2 * ($qz*$x - $qx*$z);
    my $tz = 2 * ($qx*$y - $qy*$x);

    return (
        $x + $w*$tx + ($qy*$tz - $qz*$ty),
        $y + $w*$ty + ($qz*$tx - $qx*$tz),
        $z + $w*$tz + ($qx*$ty - $qy*$tx),
    );
}

sub rotate_vertices {
    my ($verts,$q) = @_;
    return [ map { [rotate_one($q,$_)] } @$verts ];
}

sub score_orientation {
    my ($verts,$q,$cutoff) = @_;
    my $worst = 1e99;

    for my $v (@$verts) {
        my ($x,$y,$z) = rotate_one($q,$v);
        return -1e30 unless is_ocean_vec($x,$y,$z);
        my $d = nearest_coast_km($x,$y,$z);
        $worst = $d if $d < $worst;
        return $worst if $worst <= $cutoff;
    }
    return $worst;
}

sub insert_top {
    my ($top,$score,$q,$limit) = @_;
    push @$top, { score=>$score, q=>[@$q] };
    @$top = sort { $b->{score} <=> $a->{score} } @$top;
    splice(@$top, $limit) if @$top > $limit;
}

sub checkpoint_path {
    my ($name) = @_;
    return "$out_dir/.$name.checkpoint";
}

sub save_checkpoint {
    my ($name,$solid_seed,$done,$top,$phase) = @_;
    my $path = checkpoint_path($name);
    my $tmp = "$path.tmp.$$";

    open my $fh, '>', $tmp or die "open $tmp: $!\n";
    print $fh "version\t1\n";
    print $fh "name\t$name\n";
    print $fh "seed\t$seed\n";
    print $fh "solid_seed\t$solid_seed\n";
    print $fh "trials\t$trials\n";
    print $fh "keep\t$keep\n";
    print $fh "done\t$done\n";
    print $fh "phase\t$phase\n";

    for my $cand (@$top) {
        printf $fh "candidate\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
            $cand->{score}, @{$cand->{q}};
    }

    close $fh or die "close $tmp: $!\n";
    rename $tmp, $path or die "rename $tmp -> $path: $!\n";
}

sub load_checkpoint {
    my ($name,$solid_seed) = @_;
    my $path = checkpoint_path($name);
    return unless -s $path;

    open my $fh, '<', $path or die "open $path: $!\n";

    my %meta;
    my @top;

    while (my $line = <$fh>) {
        chomp $line;
        my @f = split /\t/, $line;

        if ($f[0] eq 'candidate') {
            die "$path: malformed candidate\n" unless @f == 6;
            push @top, {
                score => 0 + $f[1],
                q => [ map { 0 + $_ } @f[2..5] ],
            };
        }
        else {
            $meta{$f[0]} = $f[1];
        }
    }
    close $fh;

    die "$path: unsupported checkpoint version\n"
        unless defined($meta{version}) && $meta{version} == 1;
    die "$path: checkpoint is for $meta{name}, not $name\n"
        unless defined($meta{name}) && $meta{name} eq $name;
    die "$path: --seed changed ($meta{seed} -> $seed); remove checkpoint or use matching --seed\n"
        unless defined($meta{seed}) && $meta{seed} == $seed;
    die "$path: --trials changed ($meta{trials} -> $trials); remove checkpoint or use matching --trials\n"
        unless defined($meta{trials}) && $meta{trials} == $trials;
    die "$path: --keep changed ($meta{keep} -> $keep); remove checkpoint or use matching --keep\n"
        unless defined($meta{keep}) && $meta{keep} == $keep;
    die "$path: internal solid seed changed\n"
        unless defined($meta{solid_seed}) && $meta{solid_seed} == $solid_seed;

    @top = sort { $b->{score} <=> $a->{score} } @top;
    splice(@top, $keep) if @top > $keep;

    my $done = int($meta{done} // 0);
    my $phase = $meta{phase} // 'random';

    return ($done,\@top,$phase);
}

sub clear_checkpoint {
    my ($name) = @_;
    my $path = checkpoint_path($name);
    unlink $path if -e $path;
}

sub solid_seed {
    my ($solid_index) = @_;
    # Stable independent PRNG stream for each solid.
    return ($seed + 1_000_003 * ($solid_index + 1)) & 0x7fffffff;
}

sub optimize_orientation {
    my ($name,$verts,$solid_index) = @_;
    my @top;

    my $solid_seed = solid_seed($solid_index);
    my $start_i = 0;
    my $phase = 'random';

    my @loaded = load_checkpoint($name,$solid_seed);
    if (@loaded) {
        ($start_i,my $topref,$phase) = @loaded;
        @top = @$topref;
        printf "  resuming checkpoint: %d/%d random trials complete; best %.3f km\n",
            $start_i, $trials, (@top ? $top[0]{score} : -1e30);
    }

    srand($solid_seed);

    if (!@top) {
        # Include the canonical orientation as a deterministic candidate.
        my $identity = [1,0,0,0];
        my $s0 = score_orientation($verts,$identity,-1e31);
        insert_top(\@top,$s0,$identity,$keep);
    }

    # Perl's rand() state is not portable/serializable.  Re-create the exact
    # stream by replaying only q_random() calls up to the checkpoint.  This is
    # intentionally cheap: no coastline scoring is repeated.
    for (1 .. $start_i) {
        q_random();
    }

    if ($phase eq 'random' && $start_i < $trials) {
        my $search_start = time();
        my $last_report  = $search_start;
        my $last_checkpoint = $search_start;
        my $report_every = 5.0;

        printf "  random %8d/%d  %6.2f%%  best %.3f km\n",
            $start_i, $trials, 100.0 * $start_i / $trials, $top[0]{score};

        for my $i ($start_i + 1 .. $trials) {
            my $threshold = @top >= $keep ? $top[-1]{score} : -1e31;
            my $q = q_random();
            my $score = score_orientation($verts,$q,$threshold);
            insert_top(\@top,$score,$q,$keep)
                if $score > $threshold || @top < $keep;

            my $now = time();

            if (($now - $last_checkpoint) >= CHECKPOINT_SECONDS) {
                save_checkpoint($name,$solid_seed,$i,\@top,'random');
                printf "  checkpoint saved at %d/%d (%.2f%%)\n",
                    $i, $trials, 100.0 * $i / $trials;
                $last_checkpoint = $now;
            }

            if (($now - $last_report) >= $report_every || $i == $trials) {
                my $this_run_done = $i - $start_i;
                my $elapsed = $now - $search_start;
                my $rate = $elapsed > 0 ? $this_run_done / $elapsed : 0;
                my $remain = $rate > 0 ? ($trials - $i) / $rate : 0;
                my $pct = 100.0 * $i / $trials;

                my $eta;
                if ($rate > 0) {
                    my $hours = int($remain / 3600);
                    my $mins  = int(($remain - $hours * 3600) / 60);
                    my $secs  = int($remain - $hours * 3600 - $mins * 60 + 0.5);
                    if ($hours) {
                        $eta = sprintf("%dh%02dm%02ds", $hours, $mins, $secs);
                    }
                    elsif ($mins) {
                        $eta = sprintf("%dm%02ds", $mins, $secs);
                    }
                    else {
                        $eta = sprintf("%ds", $secs);
                    }
                }
                else {
                    $eta = "?";
                }

                printf "  random %8d/%d  %6.2f%%  best %.3f km  %.1f/s  ETA %s\n",
                    $i, $trials, $pct, $top[0]{score}, $rate, $eta;

                $last_report = $now;
            }
        }

        save_checkpoint($name,$solid_seed,$trials,\@top,'refine');
        $phase = 'refine';
    }
    elsif ($phase eq 'random' && $start_i >= $trials) {
        save_checkpoint($name,$solid_seed,$trials,\@top,'refine');
        $phase = 'refine';
    }

    my @steps = (5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01);
    my $best = $top[0];

    for my $cand (@top) {
        my $q = [@{$cand->{q}}];
        my $score = score_orientation($verts,$q,-1e31);

        for my $deg (@steps) {
            my $changed = 1;
            my $rounds = 0;
            while ($changed && $rounds++ < 50) {
                $changed = 0;
                my @axes = (
                    [ 1,0,0], [-1,0,0], [0, 1,0], [0,-1,0], [0,0, 1], [0,0,-1],
                    normalize_axis(1,1,0), normalize_axis(1,-1,0),
                    normalize_axis(1,0,1), normalize_axis(1,0,-1),
                    normalize_axis(0,1,1), normalize_axis(0,1,-1),
                );
                for my $a (@axes) {
                    my $dq = q_axis_angle(@$a,$deg);
                    my $try = q_normalize(q_mul($dq,$q));
                    my $s = score_orientation($verts,$try,$score);
                    if ($s > $score + 1e-7) {
                        $q = $try;
                        $score = $s;
                        $changed = 1;
                    }
                }
            }
        }

        if ($score > $best->{score}) {
            $best = { score=>$score, q=>$q };
        }
    }

    printf "  refined       best %.3f km\n", $best->{score};
    return ($best->{q}, $best->{score});
}

sub normalize_axis {
    my ($x,$y,$z) = @_;
    my @v = normalize3($x,$y,$z);
    return \@v;
}

sub geo_to_pixel {
    my ($lon,$lat) = @_;
    while ($lon < -180) { $lon += 360 }
    while ($lon >= 180) { $lon -= 360 }
    $lat = clamp($lat,-90,90);

    my $x = int(($lon + 180.0) * 10.0 + 0.5);
    $x %= XW;
    my $y = int((90.0 - $lat) * 10.0 + 0.5);
    $y = 0 if $y < 0;
    $y = XH-1 if $y >= XH;
    return ($x,$y);
}

sub set_pixel {
    my ($pixref,$x,$y,$ch) = @_;
    return if $x < 0 || $x >= XW || $y < 0 || $y >= XH;
    substr($$pixref, $y*XW+$x, 1, $ch);
}

sub set_geo_pixel {
    my ($pixref,$lon,$lat,$ch) = @_;
    my ($x,$y) = geo_to_pixel($lon,$lat);
    set_pixel($pixref,$x,$y,$ch);
}

sub draw_grid {
    my ($pixref) = @_;

    # Horizontal 10-degree latitude lines.  At each pole all longitudes are one
    # physical point, so emit only one pixel there.
    for (my $y=0; $y<XH; $y+=100) {
        if ($y == 0 || $y == XH-1) {
            set_pixel($pixref, XW/2, $y, 'b');
        }
        else {
            substr($$pixref, $y*XW, XW, 'b' x XW);
        }
    }

    # Vertical 10-degree meridians, excluding the two pole rows.
    for (my $x=0; $x<XW; $x+=100) {
        for my $y (1 .. XH-2) {
            set_pixel($pixref,$x,$y,'b');
        }
    }
}

sub draw_vec_arc {
    my ($pixref,$a,$b,$ch) = @_;
    my $dot = clamp($a->[0]*$b->[0] + $a->[1]*$b->[1] + $a->[2]*$b->[2], -1, 1);
    my $ang = acos($dot);
    return if $ang < 1e-12;
    my $steps = ceil(($ang * 180.0 / PI) / 0.05);
    my $sa = sin($ang);

    for my $i (0 .. $steps) {
        my $t = $i / $steps;
        my ($x,$y,$z);
        if ($i == 0) {
            ($x,$y,$z) = @$a;
        }
        elsif ($i == $steps) {
            ($x,$y,$z) = @$b;
        }
        else {
            my $wa = sin((1-$t)*$ang)/$sa;
            my $wb = sin($t*$ang)/$sa;
            ($x,$y,$z) = normalize3(
                $wa*$a->[0] + $wb*$b->[0],
                $wa*$a->[1] + $wb*$b->[1],
                $wa*$a->[2] + $wb*$b->[2],
            );
        }
        my ($lon,$lat) = vec_to_lonlat($x,$y,$z);
        set_geo_pixel($pixref,$lon,$lat,$ch);
    }
}

sub write_csv {
    my ($path,$name,$q,$score,$verts) = @_;
    open my $fh, '>', $path or die "open $path: $!\n";
    print $fh "solid,vertex,latitude,longitude,coast_distance_km,ocean,maximin_km,qw,qx,qy,qz\n";
    for my $i (0 .. $#$verts) {
        my ($x,$y,$z) = @{$verts->[$i]};
        my ($lon,$lat) = vec_to_lonlat($x,$y,$z);
        my $d = nearest_coast_km($x,$y,$z);
        my $ocean = is_ocean_vec($x,$y,$z) ? 1 : 0;
        printf $fh "%s,%d,%.8f,%.8f,%.3f,%d,%.3f,%.12g,%.12g,%.12g,%.12g\n",
            $name,$i,$lat,$lon,$d,$ocean,$score,@$q;
    }
    close $fh;
}

sub write_xpm {
    my ($path,$name,$base,$verts,$edges) = @_;
    my $pix = $base;

    # Polyhedron edges as great-circle arcs in light gray.
    for my $e (@$edges) {
        draw_vec_arc(\$pix, $verts->[$e->[0]], $verts->[$e->[1]], 'e');
    }

    # Control points last, so they remain visible where lines meet.
    for my $v (@$verts) {
        my ($lon,$lat) = vec_to_lonlat(@$v);
        set_geo_pixel(\$pix,$lon,$lat,'o');
    }

    open my $fh, '>', $path or die "open $path: $!\n";
    print $fh "/* XPM */\n";
    print $fh "/* $name; coastline data (c) OpenStreetMap contributors, ODbL */\n";
    print $fh "static char * XFACE[] = {\n";
    print $fh '"', XW, ' ', XH, " 5 1\",\n";
    print $fh '"  c black",', "\n";
    print $fh '"b c blue",', "\n";
    print $fh '"c c #A9A9A9",', "\n";  # DarkGray coastline
    print $fh '"e c #D3D3D3",', "\n";  # LightGray polyhedron edges
    print $fh '"o c white",', "\n";

    for my $y (0 .. XH-1) {
        my $row = substr($pix, $y*XW, XW);
        print $fh '"', $row, '"';
        print $fh ',' if $y != XH-1;
        print $fh "\n";
    }
    print $fh "};\n";
    close $fh;
}

# vim:set ai softtabstop=4 shiftwidth=4 tabstop=4 expandtab: ff=unix
