#! /opt/local/bin/perl
#
# calculate refracted travel time for 1-D velocity model at local/regional distance for flat earth
#
# requires input velocity file ; distance ; source depth ; phase [P|S] ; and verbosity
# tt_refracted.pl velfile 100 10 P 0
#
# velocity file: layer thickness, vp, [vs] for arbitrary number of layers and last being half space
#
# verbose = 0: only fastest refracted time (and interface it travels) is given
# verbose > 0: all possible refracted paths and times are given
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
# read velocity model ...
$i=1;
open (FILEIN, $filein) || die "Can't open in-file $!";
while (<FILEIN>)
   {
   chomp;
   ($D[$i],$vp[$i],$vs[$i]) = split(' ',$_);
   if ($vs[$i] == 0) { $vs[$i] = $vp[$i]/1.732 ; }
   $p[$i] = 1/$vp[$i];
   #print "$i $D[$i] $vp[$i] $vs[$i] $p[$i]\n";
   $i++;
   }
$i--;
if ($verb > 0) {print "Number of layers: $i\n";}

# Calculate S travel time?
if (substr($phas,0,1) eq "S" || substr($phas,0,1) eq "s")
   {
   if ($verb > 0) {print "phase is S\n";}
   $j=1;
   while ($j <= $i) { $vp[$j] = $vs[$j] ; $p[$j] = 1./$vs[$j] ; $j++ ; }
   }

# depth = 0: source in first layer, down- and up-parts are the same
#
if ($depth == 0)
   {
   $ttref = 999;     # travel time of refracted wave
   $vmax = $vp[1];
   # calculate travel time for each interface and select shortest time as headwave traveltime
   $j=2;
   while ($j <= $i) {
      if ($verb >0) {print "\nhead wave along interface $j with velocity $vp[$j]:\n";}

      if ($vp[$j] < $vmax) {if ($verb > 0) {print " No head wave in layer $j\n" ;} $j++ ; next ; }
      $vmax = $vp[$j];

      $k=1;
      $xcrit[$j]=0;
      while ($k < $j) {
         # calculate take off and refracted wave angles along the ray path
         $theta[$k] = asin($p[$j]/$p[$k]);                 #units of rad!
         if ($k == 1) { $tof[$j] = $theta[$k]/$pi*180; }   #units of degree!
         # critical distance xcrit
         $xcrit[$j] += 2 * $D[$k]*tan($theta[$k]);
         #print "$theta[$k] $tof[$j] $xcrit[$j]\n";
         $k++;
         }
      if ($verb > 0) {printf " take off angle %4.1f and critical distance %5.1f\n",$tof[$j],$xcrit[$j];}

      # calculate travel time only if distance exceeds the critical distance!
      if ($dist > $xcrit[$j]) {
         $tt[$j] = $dist * $p[$j];
         $k=1;
         while ($k < $j) {
            $tt[$j] += 2 * $D[$k] * $p[$k] * sqrt(1 - ($p[$j]/$p[$k])*($p[$j]/$p[$k]));
            $k++;
            }
         }
      else {
         $tt[$j] = 999;
         }
      if ($ttref > $tt[$j]) { $ttref = $tt[$j]; $tofa = $tof[$j]; $interf = $j ; $interf_v = $vp[$j] ; $interf_p = $p[$j]; }
      if ($verb > 0) {printf " travel time %6.2f take-off %5.2f min_tt: %6.2f\n",$tt[$j],$tof[$j],$ttref;}
      $j++;
      }
   if ($verb > 0) {printf "\n%6.2f %4.0f %7.3f %6.2f %6.2f %6.3f %2d %4.2f\n",0,$depth,0,$tofa,$dist,$ttref,$interf,$interf_v;}
   if ($verb > 0) {printf "%6.2f %4.0f %7.3f\n",$dist,0,$ttref;}

   # now calculate path and time spent in each layer for fastest headwave
   $k=1;
   if ($verb > 0) {print "\nPath (total and horizontal) and time for fastest headwave traveling in layer $interf:\n";}
   while ($k < $interf) {
      $length[$k] = 2 * $D[$k] / (sqrt(1 - ($p[$interf]/$p[$k])*($p[$interf]/$p[$k])));
      $horlen[$k] = $length[$k] * $interf_p / $p[$k];
      $timlen[$k] = $length[$k] * $p[$k];
      #printf "path horizontal_path and time spent in layer %2d:\n",$k;
      if ($verb >0) {printf "layer %2d %6.2f %6.2f %6.3f\n",$k,$length[$k],$horlen[$k],$timlen[$k];}
      $totl += $length[$k];
      $toth += $horlen[$k];
      $tott += $timlen[$k];
      $k++;
      }
   $length[$interf] = $dist - $toth;
   $timlen[$interf] = $ttref - $tott;
   #printf "path horizontal_path and time spent in refracted layer %2d:\n",$interf;
   if ($verb >0) {printf "layer %2d %6.2f %6.2f %6.3f\n",$interf,$length[$interf],$length[$interf],$timlen[$interf];}

   }

# depth > 0: source not at surface
#
elsif ($depth > 0)
   {
   #print "$depth    num-layers: $i\n";
   $k=1;
   $sl=0;
   SL:
   while ($k < $i) {
      $TD += $D[$k];
      if ($depth < $TD) { $sl = $k ; }
      last SL if ($sl>0) ;
      $k++;
      }
   if ($sl == 0) { $sl = $i ; if ($verb > 0) {print "Source is in half-space. No refracted wave.\n"; } }
   else {
   $DD[$sl] = $TD - $depth;   #vertical distance from source to first interface below
   if ($verb >0 ) {print "Source is in layer $sl\n";}
   #print "source-depth: $depth total-depth: $TD layer-thickness: $D[$k] source-2-interface: $DD[$sl]\n";

   $ttref = 999;     # travel time of refracted wave
   $vmax = $vp[$sl];
   # calculate travel time for each interface and select shortest time as headwave traveltime
   if ($DD[$sl] == $D[$sl]) {$j=$sl;}
   else {$j=$sl+1;}
   $tt[$j]=999;
   while ($j <= $i) {
      if ($verb >0) {print "\nhead wave along interface $j with velocity $vp[$j]:\n";}
      if ($vp[$j] < $vmax) {if ($verb >0) {print " No head wave in layer $j\n" ;} $j++ ; next ; }
      $vmax = $vp[$j];
      $k=1;
      #$xcrit[$j]=0;
      while ($k < $j) {
         # calculate take off and refracted wave angles along the ray path
         $theta[$k] = asin($p[$j]/$p[$k]);                 #units of rad!
         if ($k == $sl) { $tof[$j] = $theta[$sl]/$pi*180; }   #units of degree!
         $k++;
         }
      if (($DD[$sl] == $D[$sl]) && ($sl == $j)) {$tof[$j] = 90;}
      if ($verb >0) {printf " take off angle %4.1f\n",$tof[$j];}
      # calculate travel time and horizontal distance $toth for down- and up-going part. Is $toth < total-distance $dist?
      # only then does a refracted wave exist ...
      # downgoing part:
      $dum=$sl;
      $totl=$toth=$tott=0;
      while ($dum < $j) {
         if ($dum > $sl) { $DD[$dum] = $D[$dum]; }
         $dlength[$dum] = $DD[$dum] / (sqrt(1 - ($p[$j]/$p[$dum])*($p[$j]/$p[$dum])));
         $dhorlen[$dum] = $dlength[$dum] * $p[$j] / $p[$dum];
         $dtimlen[$dum] = $dlength[$dum] * $p[$dum];
         #printf "path horizontal_path and time spent in layer %2d:\n",$dum;
         #printf "D layer %2d %6.2f %6.2f %6.3f %6.3f\n",$dum,$dlength[$dum],$dhorlen[$dum],$dtimlen[$dum],$DD[$dum];
         #$totl += $length[$dum];
         #$toth += $horlen[$dum];
         #$tott += $timlen[$dum];
         $dum++;
         }
      # upgoing part:
      $k=1;
      while ($k < $j) {
         $ulength[$k] = $D[$k] / (sqrt(1 - ($p[$j]/$p[$k])*($p[$j]/$p[$k])));
         $uhorlen[$k] = $ulength[$k] * $p[$j] / $p[$k];
         $utimlen[$k] = $ulength[$k] * $p[$k];
         #printf "path horizontal_path and time spent in layer %2d:\n",$k;
         #printf "U layer %2d %6.2f %6.2f %6.3f %6.3f\n",$k,$ulength[$k],$uhorlen[$k],$utimlen[$k],$D[$k];
         $totl += ($dlength[$k] + $ulength[$k]);
         $toth += ($dhorlen[$k] + $uhorlen[$k]);
         $tott += ($dtimlen[$k] + $utimlen[$k]);
         $k++;
         }
      # does refracted wave exist?
      $exist = 1;
      if ($toth > $dist) { $tt[$j] = 999; $exist = 0; if ($verb >0) {print " distance less than critical distance - no head wave\n"; } }
      else {$tott += ($dist - $toth)*$p[$j] ; $tt[$j] = $tott; } 

      if ($ttref > $tt[$j]) { $ttref = $tt[$j]; $tofa = $tof[$j]; $interf = $j ; $interf_v = $vp[$j] ; $interf_p = $p[$j]; }
      if ($verb >0) {printf " travel time %6.2f take-off %5.2f min_tt: %6.2f\n",$tt[$j],$tof[$j],$ttref;}

      $j++;
      }

   if ($verb > 0) {printf "\n%6.2f %4.0f %7.3f %6.2f %6.2f %6.3f %2d %4.2f\n",0,$depth,0,$tofa,$dist,$ttref,$interf,$interf_v;}
   if ($verb > 0) {printf "%6.2f %4.0f %7.3f\n",$dist,$depth,$ttref;}
   
   # path and time spent in each layer for fastest headwave - along interface $interf
   if ($verb >0) {print "\nPath (total and horizontal) and time in each layer for fastest headwave traveling in layer $interf:\n";}
   $dum=$sl;
   $totl=$toth=$tott=0;
   while ($dum < $interf) {
      if ($dum > $sl) { $DD[$dum] = $D[$dum]; }
      $dlength[$dum] = $DD[$dum] / (sqrt(1 - ($p[$interf]/$p[$dum])*($p[$interf]/$p[$dum])));
      $dhorlen[$dum] = $dlength[$dum] * $p[$interf] / $p[$dum];
      $dtimlen[$dum] = $dlength[$dum] * $p[$dum];
      $dum++;
      }
   $k=1;
   while ($k < $interf) {
      $ulength[$k] = $D[$k] / (sqrt(1 - ($p[$interf]/$p[$k])*($p[$interf]/$p[$k])));
      $uhorlen[$k] = $ulength[$k] * $p[$interf] / $p[$k];
      $utimlen[$k] = $ulength[$k] * $p[$k];
      if ($verb > 0) {printf "layer %2d %6.2f %6.2f %6.3f %5.2f\n",$k,$dlength[$k]+$ulength[$k],$dhorlen[$k]+$uhorlen[$k],$dtimlen[$k]+$utimlen[$k],$vp[$k];}
      $totl += ($dlength[$k] + $ulength[$k]);
      $toth += ($dhorlen[$k] + $uhorlen[$k]);
      $tott += ($dtimlen[$k] + $utimlen[$k]);
      $k++;
      }
   $length[$interf] = $dist - $toth;
   $timlen[$interf] = $ttref - $tott;
   if ($verb > 0) {printf "layer %2d %6.2f %6.2f %6.3f %5.2f\n",$interf,$length[$interf],$length[$interf],$timlen[$interf],$vp[$interf];}

   }
   }

if ($verb > 0) {printf "\n";}
if ($ttref == 999)
   { printf "NO HEAD WAVE EXISTS at Distance[km] %6.1f -- Phase %2s -- Depth[km] %5.1f\n",$dist,$phas,$depth; }
else { printf "FRA Distance[km] %6.1f -- Phase %2s -- TravelTime[s] %7.3f -- Takeoff[Deg] %5.1f -- Depth[km] %5.1f -- Layer %2d\n",$dist,$phas,$ttref,$tofa,$depth,$interf; }

close (FILEIN);
