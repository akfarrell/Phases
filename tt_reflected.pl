#! /opt/local/bin/perl -w
#
# calculate reflected travel time for 1-D velocity model with reflection from interface N (with N-1 layers above interface)
#
# requires input velocity file ; distance ; source depth ; phase [P|S] ; interface ; and verbosity
# tt_reflected.pl velfile 100 10 P 1 0
#
# velocity file: layer thickness, vp, [vs] for arbitrary number of layers and last being half space
#
# interface: interface at which ray reflects ; if [M|m] calculates Moho reflected phase
#
# no output results for verbosity = 0 when no reflected phase exists ...
#
# version 1.0 01/08/13 jb
#
use Math::Trig;
use Math::Complex;
use List::Util qw(max);  #use: $a = max 1 , 2 ... then $a=2
$pi=3.14159265;

# verbose: 0 no ; >0 yes
$verb=0;
# number of trials to find correct take off angle for given distance
$maxcnt=100;

# input: velocity model, distance and P or S phase to calculate travel time for.
# Should be parsed from a main routine ...
$filein = $ARGV[0];
$dist = $ARGV[1];
$depth = $ARGV[2];
$phas = $ARGV[3];
$inte = $ARGV[4];
$verb = $ARGV[5];
#print "$filein $dist $depth $phas\n";
# read velocity model ...
$i=1;
open (FILEIN, $filein) || die "Can't open in-file $!";
while (<FILEIN>)
   { chomp;
   ($D[$i],$vp[$i],$vs[$i]) = split(' ',$_);
   if ($vs[$i] == 0) { $vs[$i] = $vp[$i]/1.732 ; }
   #print "$i $D[$i] $vp[$i] $vs[$i]\n";
   $i++; }
$i--;
# calculate S travel time?
if (substr($phas,0,1) eq "S" || substr($phas,0,1) eq "s")
   {
   if ($verb > 0) {print "phase is S\n";}
   $j=1;
   while ($j <= $i) { $vp[$j] = $vs[$j] ; $j++ ; }
   }
# which interface? answer: $inte ; $temp is 1 for Moho reflection
$temp = 0;
if ($inte eq "M" || $inte eq "m")
   { $temp=1;
   if ($verb > 0) {print "Moho reflection - assume mantle is where vp>=7.5 and vs>=4.2 km/s\n";}
   $inte=1;
   while (($vp[$inte+1] < 7.5) && ($vs[$inte+1] < 4.2)) {$inte++;}
   }
elsif (($vp[$inte+1] >= 7.5) && ($vs[$inte+1] >= 4.2) && ($vp[$inte] < 7.5) && ($vs[$inte] < 4.2) )
   { $temp=1; }
elsif ($inte >= $i)
   { print "Reflection from interface $inte impossible because only $i layers, including half space, exist\n"; exit; }
if ($verb > 0) {print "Reflection from interface $inte\n";}
# in which layer is source? answer: $sl
$k=1;
$sl=0;
SL:
while ($k < $i) {
   $TD += $D[$k];
   if ($depth < $TD) { $sl = $k ; }
   last SL if ($sl>0) ;
   $k++;
   }
if ($sl == 0) { $sl = $i ; if ($verb > 0) {print "Source is in half-space. No reflected wave.\n"} ; exit; }
elsif ($sl > $inte) {if ($verb > 0) {print "Source is in layer $sl below interface $inte specified. No reflected wave.\n"; } exit; }
else { $DD[$sl] = $TD - $depth;   #vertical distance from source to next interface below
   if ($verb > 0) {print "Source is in layer $sl\n";}
   }
#print "DD $DD[$sl]\n";

# Find correct take off and refraction angles along ray ... $tup is take off angle at source in layer $sl
$tup=0;          # in degree
$step=4;
$tol=0.00001;    # relative tolerance 0.01 is within 1%
$xtest=0;        # horizontal distance for a given take off angle $tup
$cnt=1;
while (abs ($xtest-$dist) > ($dist*$tol)) {
   if ($cnt == $maxcnt) { last } ;
   $tupr = $tup / 180. * $pi;   #convert to radians
   #print "$tupr $tup\n";
   $pcr=0;
   for ($nl = 1; $nl <= $inte; $nl++) {
      if ((sin($tupr) / $vp[$sl] * $vp[$nl]) > 1) { $pcr=1 ; last ; } # refracted angle would be post-critical, jump out of for-loop
      $ang[$nl] = asin ( sin($tupr) / $vp[$sl] * $vp[$nl]) ;
      #$ang[$nl] = Re (asin ( sin($tupr) / $vp[$sl] * $vp[$nl])) ;
      #print "cnt $cnt  layer $nl  angle $ang[$nl]\n";
      }

   if ($pcr == 1) { 
      # if statement: should not be reachable since ray is vetical for $cnt==1 ...?
      if ($cnt == 1) { print "NO REFLECTION from interface $inte does not exist because it goes post-critical along the path!\n" ; exit ; }
      # else: the step is too large; reduce step size, reduce take-off angle, and go back to try with the new take-off angle ...
      else { if ($verb > 0) { print "Refraction is becoming postcritical; reduce take-off angle, try again ... \n" } ; $step = $step/1.99 ; $tup -= $step ; $cnt++ ; next ; }
      }

   # we should get here only for pre-critical angles, so sqrt() is a real number ...
   $xsl = ($DD[$sl] * sin ($ang[1]) * $vp[$sl]/$vp[1]) / ( sqrt (1 - ( $vp[$sl]*$vp[$sl]*sin($ang[1])*sin($ang[1])/($vp[1]*$vp[1]))));
   $xtest = $xsl ;
   for ($nl = 1; $nl <= $sl; $nl++) { $x[$nl] = ($D[$nl] * sin ($ang[1]) * $vp[$nl]/$vp[1]) / ( sqrt (1-($vp[$nl]*$vp[$nl]*sin($ang[1])*sin($ang[1])/($vp[1]*$vp[1])))); $xtest += $x[$nl] ; }
   for ($nl = ($sl+1) ; $nl <= $inte; $nl++) { $x[$nl] = ($D[$nl] * sin ($ang[1]) * $vp[$nl]/$vp[1]) / ( sqrt (1-($vp[$nl]*$vp[$nl]*sin($ang[1])*sin($ang[1])/($vp[1]*$vp[1])))); $xtest += (2*$x[$nl]) ; }
   #$xtest = abs ($xtest);  #xtest though is real ...
   #print "$cnt $xtest $dist $ang[$inte]\n";

   $prev_tup = $tup;
   if ($xtest > $dist) {$tup -= $step ; $step = $step/1.99 ; }
   $tup += $step;
   
   $cnt++;
   }

if ($cnt == $maxcnt) { printf "WARNING correct distance %6.2f [km] was not found within $maxcnt trials! Distance used is %6.2f [km]\n",$dist,$xtest; }
# xtest is within specs of dist -- but tup has been updated afterwards, revert to previous value!
$tup = $prev_tup;
$tupr = $tup / 180. * $pi;   #convert to radians
#print "$tupr $tup $xtest $dist\n";

# Calculate $x[$nl] for correct take off angle ...
for ($nl = 1; $nl <= $inte; $nl++)
   {
   $ang[$nl] = asin ( sin($tupr) / $vp[$sl] * $vp[$nl]);          # should all be real now!
   #$ang[$nl] = abs (asin ( sin($tupr) / $vp[$sl] * $vp[$nl])) ;  # abs only necessary if asin() is complex ...
   if ($ang[$nl] >= 1.57079) { print "NO REFLECTION POSSIBLE from interface $inte because post-critical angle of refraction along reflection-path!\n" ; exit; }
   }

$PCR = "" ;
if ( (sin($ang[$inte])*$vp[$inte+1]/$vp[$inte]) > 1 ) {$PCR = "PC";}

$ttref = 0;
if ($verb > 0) {print "\nHorizontal path and time in each layer for phase reflected at interface $inte:\n";}
#for ($nl = 1; $nl <= $sl; $nl++)
#   {
#   $x[$nl] = ($D[$nl] * sin ($ang[1]) * $vp[$nl]/$vp[1]) / ( sqrt (1-($vp[$nl]*$vp[$nl]*sin($ang[1])*sin($ang[1])/($vp[1]*$vp[1]))));
#   $t[$nl] = sqrt ($D[$nl]*$D[$nl] + $x[$nl]*$x[$nl]) / $vp[$nl] ;
#   if ($verb > 0) {printf "layer %2d horizontal path %5.2f -- travel time %5.2f -- above source, traversed once upwards\n",$nl,$x[$nl],$t[$nl];}
#   $ttref += $t[$nl];
#   }

$xsl = ($DD[$sl] * sin ($ang[1]) * $vp[$sl]/$vp[1]) / ( sqrt (1 - ( $vp[$sl]*$vp[$sl]*sin($ang[1])*sin($ang[1])/($vp[1]*$vp[1]))));
$tsl = sqrt ($DD[$sl]*$DD[$sl] + $xsl*$xsl) / $vp[$sl] ;
if ($verb > 0) {printf "layer %2d horizontal path %6.2f -- travel time %6.2f -- source layer downgoing\n",$sl,$xsl,$tsl;}
$ttref += $tsl ;
#print "$xsl $tsl $ttref\n";

for ($nl = ($sl+1) ; $nl <= $inte; $nl++)
   {
   $x[$nl] = ($D[$nl] * sin ($ang[1]) * $vp[$nl]/$vp[1]) / ( sqrt (1-($vp[$nl]*$vp[$nl]*sin($ang[1])*sin($ang[1])/($vp[1]*$vp[1])))); 
   $t[$nl] = sqrt ($D[$nl]*$D[$nl] + $x[$nl]*$x[$nl]) / $vp[$nl] ;
   if ($verb > 0) {printf "layer %2d horizontal path %6.2f -- travel time %6.2f -- below source, traversed twice\n",$nl,(2*$x[$nl]),(2*$t[$nl]);}
   $ttref += (2*$t[$nl]);
   }

for ($nl = $sl; $nl >= 1; $nl--)
   {
   $x[$nl] = ($D[$nl] * sin ($ang[1]) * $vp[$nl]/$vp[1]) / ( sqrt (1-($vp[$nl]*$vp[$nl]*sin($ang[1])*sin($ang[1])/($vp[1]*$vp[1]))));
   $t[$nl] = sqrt ($D[$nl]*$D[$nl] + $x[$nl]*$x[$nl]) / $vp[$nl] ;
   if ($verb > 0) {printf "layer %2d horizontal path %6.2f -- travel time %6.2f -- above source, traversed once upwards\n",$nl,$x[$nl],$t[$nl];}
   $ttref += $t[$nl];
   }

$tofa = $tup;
if ($temp == 1) {$inte = "M" ;}
$inte = $inte.$PCR ;
#print "Total travel time of reflected phase is $ttref\n";
printf "FLE Distance[km] %6.1f -- Phase %2s -- TravelTime[s] %7.3f -- Takeoff[Deg] %5.1f -- Depth[km] %5.1f -- Layer %2d -- Interface %3s\n",$xtest,$phas,$ttref,$tofa,$depth,$sl,$inte;
close (FILEIN);
