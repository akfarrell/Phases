% if sec == 1
%     val_try = P_ind; indie = 1; ph_used = 'P';
% elseif sec == 3
%     val_try = S_ind; ph_used = 'S'; indie = 3; %CHANGE!!!!
% elseif sec == 4
%     val_try = S_ind; ph_used = 'S'; indie = 4;
% elseif sec == 2
%     try
%         vee = 1;
%         val_try = ph_indz(vee);
%         indie = vee+1;
%         if vee == 1
%             ph_used = 'P1';
%         elseif vee == 2
%             ph_used = 'P2';
%         elseif vee == 3
%             ph_used = 'P3';
%         end
%     catch
%     end
% elseif sec == 5
%     try
%         vee = 2;
%         val_try = ph_indz(vee);
%         indie = vee+1;
%         if vee == 1
%             ph_used = 'P1';
%         elseif vee == 2
%             ph_used = 'P2';
%         elseif vee == 3
%             ph_used = 'P3';
%         end
%     catch
%     end
% end
%  eqs_interest = [2501 2494 2251 2311 2312 2436 1786 1951 1864 1870 2000 2042 2043 2046 2049 2050 2054 2055 ...
%     2056 2057 2119 2152 2153 2161 2165 2166 2177 2210 2217 2223 2247 2249 2267 2272 2273 2276 2295 2313 ...
%     2314 2324 2347 2382 2386 2390];
tic
sec = [1 2 3];
sec2 = [1 2 5 4];
% for count = 2%1:numel(sec)
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2501,'PLMN',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2494,'PLBR',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2311,'PLMN',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2312,'PLAR',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2436,'PLSS',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(1786,'PLLL',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(1951,'PLSQ',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(1864,'PLLL',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(1870,'PLLL',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2000,'PLLL',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2042,'PLAR',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2043,'PLLA',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2046,'PLSS',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2049,'PLSS',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2050,'PLSS',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2054,'PLSM',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2055,'PLSM',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2056,'PLSM',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2057,'PLSM',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2119,'PL03',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2152,'PLMD',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2153,'PLMD',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2161,'PLSS',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2165,'PL03',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2165,'PLCM',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2165,'PLMD',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2166,'PLBR',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2177,'PLAR',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2177,'PLMD',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2210,'PLLA',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2217,'PLSS',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2223,'PLSS',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2247,'PLSS',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2249,'PLLL',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2267,'PLBR',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2272,'PLBR',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2273,'PLBR',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2276,'PLLL',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2295,'PL03',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2313,'PLLL',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2314,'PLBR',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2324,'PLLA',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2347,'PLJR',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2382,'PLLL',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2386,'PLSM',sec(count))
%     phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2390,'PLMN',sec(count))
% end

for count = 3%1:numel(sec2)
    phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(2251,'PLLL',sec2(count))
end
toc