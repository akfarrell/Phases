#! /opt/local/bin/perl
#
# calculate direct travel time for 1-D velocity model at local/regional distances for flat earth
#
# requires input velocity file ; distance ; source depth ; phase [P|S] ; and verbosity
# tt_direct.pl velfile 100 10 P 0
#
# velocity file: layer thickness, vp, [vs] for arbitrary number of layers and last being half space
#
# version 1.0 01/08/13 jb
#
use Math::Trig;
$pi = 3.14159265;

$verb = 0;

# input: velocity model, distance and P or S phase to calculate travel time for.
# Should be parsed from a main routine ...
$filein = @ARGV[0];
$dist = @ARGV[1];
$depth = @ARGV[2];
$phas = @ARGV[3];
$verb = @ARGV[4];
#print "$filein $dist $depth $phas\n";
# read veolocity model ...
$i=1;
open (FILEIN, $filein) || die "Can't open in-file $!";
while (<FILEIN>)
   {
   chomp;
   ($D[$i],$vp[$i],$vs[$i]) = split(' ',$_);
   if ($vs[$i] == 0) { $vs[$i] = $vp[$i]/1.732 ; }
   #print "$i $D[$i] $vp[$i] $vs[$i]\n";
   $i++;
   }
$i--;

# Calculate S travel time?
if (substr($phas,0,1) eq "S" || substr($phas,0,1) eq "s")
   {
   if ($verb > 0) {print "phase is S\n";}
   $j=1;
   while ($j <= $i) { $vp[$j] = $vs[$j] ; $j++ ; }
   }

# depth = 0: ray runs simply along surface ...
if ($depth == 0)
   {
   $ttd = $dist / $vp[1];
   #$tof = (atan2 $dist,($depth+.0001))/$pi*180;
   $tof = 90;
   if ($verb > 0) { printf "%6.2f %4.0f %7.3f %6.2f %6.2f %6.3f\n",0,$depth,0,$tof,$dist,$ttd;
      printf "%6.2f %4.0f %7.3f\n",$dist,$depth,$ttd; }
   }
#
# event occurs in first layer, tt-curve is simple hyperbola ...
elsif ($depth > 0 && $depth <= $D[1])
   {
   $ttd = sqrt ($depth * $depth + $dist * $dist) / $vp[1];
   $tof = 90 + (atan2 $depth,$dist)/$pi*180;
   if ($verb > 0) { printf "%6.2f %4.0f %7.3f %6.2f %6.2f %6.3f\n",0,$depth,0,$tof,$dist,$ttd;
      printf "%6.2f %4.0f %7.3f\n",$dist,$depth,$ttd; }
   }
#
# event is in second layer or deeper ...
else
   {
   # which layer is source? -- $sl is source layer, $TD is vertical travel path in source layer
   $sl = 1;
   $TD = 0;
   while ($depth > $TD && $sl < $i) { $TD += $D[$sl] ; $sl++; }
   $sl--;
   if ($TD < $depth) { $sl = $i ; }  # source in half space
   # vertical path length in source layer
   $l = 0;
   $TD = 0;
   while ($l < $sl) { $TD += $D[$l] ; $l++; }
   $l-- ; # $l layers above source layer
   $TD = $depth - $TD;
   #print "$depth $TD $sl $l $i\n";
   # find correct take off angle for up-going ray $tup -- relative to a downgoing ray this angle is $pi-$tof
   $tup = 0;
   $step = 4;
   $tol = 0.001;
   while (abs ($xtest-$dist) > $tol) {
      $tupr = $tup / 180. * $pi;   #convert to radians
      $xtest = $TD * tan ($tupr);
      for ($nl = 1; $nl < $sl; $nl++) {
         $ang[$nl] = asin ( sin($tupr) / $vp[$sl] * $vp[$nl]) ;
         $xtest += $D[$nl] * tan ($ang[$nl]);
         }
      #print "$xtest $dist $tol $tup $step $cnt\n";
      if (($xtest > $dist) || ($xtest < 0)) { $tup -= $step ; $step = $step/2. ; }
      $tup += $step;
      if ($tup > 90) { $tup = 90.; }  # only upgoing rays are considered ...
      $cnt ++;
      if ($cnt >= 100) { if ($verb > 0) {print "$xtest $dist\n";print "Cannot find correct distance. Give up after $cnt trials\n";} exit ;}
      }
   # horizontal distance travelled and time spent in layer nl (nl < sl)
   for ($nl = ($sl-1); $nl > 0; $nl--) {
       $xl[$nl] = $D[$nl] * tan($ang[$nl]) ;
       $tl[$nl] = sqrt ($xl[$nl] * $xl[$nl] + $D[$nl] * $D[$nl])/$vp[$nl];
       }
   # horizontal distanced travelled and time spent on source layer sl
   $xl[$sl] = $TD * tan($tupr);
   $tl[$sl] = sqrt ($xl[$sl] * $xl[$sl] + $TD * $TD)/$vp[$sl];
   $ang[$sl] = $tupr;
   # cdist: cumulative distance on x direction from source ; cdepth: z-value corresponing to cdist along ray (together they
   # are the ray path ... ; $ttd: cumulative time from source along ray ...
   $cdist = 0;
   $cdepth = $depth;
   $ttd = 0;
   # express take off angle relative to vertical-down = 0 instead of vertical-up, 180 degree flip
   for ($nl = $sl ; $nl > 0 ; $nl--) {
       if ($verb > 0) {printf "%6.2f %4.0f %7.3f %6.2f %6.2f %6.3f\n",$cdist,$cdepth,$ttd,180-$ang[$nl]*180/$pi,$xl[$nl],$tl[$nl];}
       if ($nl == $sl) { $cdepth -= $TD ; }
       else { $cdepth -= $D[$nl]; }
       $cdist += $xl[$nl];
       $ttd += $tl[$nl];
       #printf "%6.2f %6.2f %6.3f\n",$ang[$nl]*180/$pi,$xl[$nl],$tl[$nl];
       }
   if ($verb > 0) {printf "%6.2f %4.0f %7.3f\n",$cdist,$cdepth,$ttd;}
   $tof = 180 - $tup;
   if ($verb > 0) {printf "Calculated Distance[km] %6.1f\n",$cdist;}
   }

printf "DIR Distance[km] %6.1f -- Phase %2s -- TravelTime[s] %7.3f -- Takeoff[Deg] %5.1f -- Depth[km] %5.1f\n",$dist,$phas,$ttd,$tof,$depth;

close (FILEIN);
