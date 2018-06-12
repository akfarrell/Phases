tic
addpath('/raid/home/a/akfarrell/')
addpath('/raid/home/a/akfarrell/Uturuncu')
%clear;
clear all; close all; clc;

SECSPERDAY = 60 * 60 * 24;
siteStruct = siteStruct_for_Heather();
stationz = {'LZ2','LZ3','LZAC','LZAE','LZAZ','LZBB','LZLA','PLLC'};

ds = datasource('antelope', '/raid/data/antelope/databases/PLUTONS/dbmerged');
scnl = scnlobject(stationz, 'HH*');
stime = datenum(2011,11,12,3,40,20);
etime = datenum(2011,11,12,3,41,35);

% stime = datenum(2011,11,12,3,40,15);
% etime = datenum(2011,11,12,3,41,05);

w = waveform(ds, scnl, stime, etime);
w_clean = waveform_clean(w,filterobject('b', [3 20], 2));
plot_panels(w_clean)
%E, N, Z
%%
w_clean_E = w_clean(1:3:end);
w_clean_N = w_clean(2:3:end);
w_clean_Z = w_clean(3:3:end);
%get(w_clean_Z,'channel')
%%
ml = 3;
% oridStruct.eq_1.lon = 171.522;
% oridStruct.eq_1.lat = 51.869;
% oridStruct.eq_1.depth = 48.9;
% oridStruct.eq_1.ml = ml;
% 
% oridStruct.eq_2.lat = 35.545;
% oridStruct.eq_2.lon = 96.730;
% oridStruct.eq_2.depth = 3.9;
% oridStruct.eq_2.ml = ml;
% 
% oridStruct.eq_3.lat = 35.532;
% oridStruct.eq_3.lon = 96.753;
% oridStruct.eq_3.depth = 6.2;
% oridStruct.eq_3.ml = ml;
% 
% oridStruct.eq_4.lon = 65.572;
% oridStruct.eq_4.lat = 18.633;
% oridStruct.eq_4.depth = 82.1;
% oridStruct.eq_4.ml = ml;

%%
% oridStruct.eq_5.lon = -68.507;
% oridStruct.eq_5.lat = -25.168;
% oridStruct.eq_5.depth = -5.706;
% oridStruct.eq_5.ml = ml;
% 
% oridStruct.eq_6.lon = -68.521;
% oridStruct.eq_6.lat = -25.336;
% oridStruct.eq_6.depth = -5.467;
% oridStruct.eq_6.ml = ml;
%%
oridStruct.eq_7.lon = 122.261;
oridStruct.eq_7.lat = 23.989;
oridStruct.eq_7.depth = 34.6;
oridStruct.eq_7.ml = 4.5;

% w_AK = prep_for_rs(w_clean_E,'eq_1',siteStruct,oridStruct)
% w_OK1 = prep_for_rs(w,'eq_2',siteStruct,oridStruct)
% w_OK2 = prep_for_rs(w,'eq_3',siteStruct,oridStruct)
% w_PR = prep_for_rs(w,'eq_4',siteStruct,oridStruct)
% 
% w_lastarriaE = prep_for_rs(w_clean_E,'eq_5',siteStruct,oridStruct);
% w_lastarriaN = prep_for_rs(w_clean_N,'eq_5',siteStruct,oridStruct);
% w_lastarriaZ = prep_for_rs(w_clean_Z,'eq_5',siteStruct,oridStruct);
% 
% w_azufreE = prep_for_rs(w_clean_E,'eq_6',siteStruct,oridStruct);
% w_azufreN = prep_for_rs(w_clean_N,'eq_6',siteStruct,oridStruct);
% w_azufreZ = prep_for_rs(w_clean_Z,'eq_6',siteStruct,oridStruct);

w_TWE = prep_for_rs(w_clean_E,'eq_7',siteStruct,oridStruct);
w_TWN = prep_for_rs(w_clean_N,'eq_7',siteStruct,oridStruct);
w_TWZ = prep_for_rs(w_clean_Z,'eq_7',siteStruct,oridStruct);

%%
% plotw_rs_test(w_AK)
% plotw_rs_test(w_OK1)
% plotw_rs_test(w_OK2)
% plotw_rs_test(w_PR)


% plotw_rs_test(w_lastarriaE)
% plotw_rs_test(w_lastarriaN)
% plotw_rs_test(w_lastarriaZ)
% plotw_rs_test(w_azufreE)
% plotw_rs_test(w_azufreN)
% plotw_rs_test(w_azufreZ)

plotw_rs_test(w_TWE)
plotw_rs_test(w_TWN)
plotw_rs_test(w_TWZ)

%toc
%lat lon depth ml