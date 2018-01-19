% function high_corrs = phase_hunt_changingNeedle_mins(allorids,oridStruct,P,fil,cutoff_val)
% find(strcmp(fieldnames(oridStruct),'eq_whatever')) %to find eq
%orid2166 is origin 203
%orid2034 is origin 128
%orid2277 is origin 277
tic
addpath('/raid/home/a/akfarrell/')
addpath('/raid/home/a/akfarrell/Uturuncu')
addpath('/raid/home/a/akfarrell/Polarizemic-master/functions/')
%clear;
clc;
fil=[2 25];
%fil=[10 20];
[oridStruct, allorids] = get_eq_info(); 
load('siteStruct.mat');
%stylie = 'min'; %%%%% CHANGE!!!---------0----min,max,abz
stylie = {'max','min'};

%ch = {'HHE','HHN','HHZ'};
ch2 = {'HHZ','HHN','HHE'};
ch = {'HHZ','HHR','HHT'};
cutoff_val = 0.6;
S_pad = 25;
SECSPERDAY = 60 * 60 * 24;
num_samps = 31;
%close all

eq = 2066;
stasz = 'PLLL';
val_try = 499;%

%%
count=find(allorids == eq);
close all
directory = '/home/a/akfarrell/Uturuncu/Phase/wf_objs';
filename = sprintf('wf_%d.mat',allorids(count));
filename_wPath = fullfile(directory,filename);
%if exist(filename_wPath,'file')
%    load(filename_wPath)
%else
    %create and clean waveform object
    [w_raw,OrS,stations_inEq] = get_wf(allorids(count),oridStruct);
    w_clean_filt = waveform_clean(w_raw, filterobject('b', fil, 2));
    w_clean = waveform_clean(w_raw);
    %save(filename_wPath,'w_clean', 'OrS', 'stations_inEq');
%end
stationz = get(w_clean,'station');

%%
count2 = min(find(strcmp(stasz,stationz)));%:3:numel(w_clean) %43-45 is PLMN 43:3:43
%HHE = count, HHN = count+1, HHZ = count+2

fname = sprintf('corr_%d_%s_%s.mat',allorids(count),stationz{count2},stylie{1});
directory = '/home/a/akfarrell/Uturuncu/Phase/wf_objs';
directory2 = sprintf('/home/a/akfarrell/Uturuncu/Phase/corrs/%d',allorids(count));
filename_wPath2 = fullfile(directory2,fname);
load(filename_wPath2)
c_orig=c;
%P_pad = 15
dnum = zeros(1,numel(get(w_clean(count2),'data')));
dnum(1) = datenum(get(w_clean(count2),'start'));
freq = get(w_clean(count2), 'freq');
for l = 2:numel(get(w_clean(count2),'data'))
    dnum(l) = datenum((l/freq)/SECSPERDAY+dnum(1));
end
numz = find(strcmp(siteStruct.sta,stasz));
az = azimuth(siteStruct.lat(numz), siteStruct.lon(numz),oridStruct.(sprintf('eq_%d',eq)).lat(1),oridStruct.(sprintf('eq_%d',eq)).lon(1));
starttime = oridStruct.(sprintf('eq_%d',eq)).time_phase(intersect(find(strcmp(oridStruct.(sprintf('eq_%d',eq)).phase,'P')),...
    find(strcmp(oridStruct.(sprintf('eq_%d',eq)).sta,stasz))))-datenum(0,0,0,0,0,1);

%% Find best filter using periodogram and filter to that
fs = get(w_clean(count2),'freq')
w_clean_copy = w_clean;
w_clean_copy = w_clean_copy(count2:count2+2);
compnt = get(w_clean_copy,'channel');
data = get(w_clean_copy,'data');
close all
errorz = zeros(3,2);
for countz = 1:numel(data)
    sig_data.(compnt{countz}) = data{countz}(val_try-5:val_try+36); 
    %added 5 samples before val_try to allow for a padding for fft
    %also added 5 samples after 32 samples, to allow for padding for fft
    leng = numel(sig_data.(compnt{countz}));
    nfft = leng;
    %[pdgrm2.(compnt{countz}),f_pdgrm2.(compnt{countz})] = periodogram(sig_data.(compnt{countz}),[],nfft,fs);
    [pdgrm1.(compnt{countz}),f_pdgrm1.(compnt{countz})] = periodogram(sig_data.(compnt{countz}),hamming(leng),nfft,fs);
    [pdgrm.(compnt{countz}),f_pdgrm.(compnt{countz})] = periodogram(detrend(sig_data.(compnt{countz})),hamming(leng),nfft,fs);
    % TRY THIS WITH DETERMINE PEAKS OF BOTH!!!!!
    [k.(compnt{countz}),w.(compnt{countz})] = det_pks(pdgrm.(compnt{countz}),f_pdgrm.(compnt{countz}),3,'pdgrm')
    %[k1.(compnt{countz}),w1.(compnt{countz})] = det_pks(pdgrm1.(compnt{countz}),f_pdgrm1.(compnt{countz}),3,'pdgrm')
    %w1.(compnt{countz})=w1.(compnt{countz})/2;
    w.(compnt{countz})=w.(compnt{countz})/2;
    errorz(countz,:) = [k.(compnt{countz})-w.(compnt{countz}),k.(compnt{countz})+w.(compnt{countz})];
    %errorz1.(compnt{countz}) = [k1.(compnt{countz})+w1.(compnt{countz}),k1.(compnt{countz})-w1.(compnt{countz})];
   
end
%%
intersect_vals = RangesIntersect(errorz(1,:),errorz(2,:) ,errorz(3,:) )
%make sine waves of frequency picked up by det_pks



%% Subset data to around the P- and S-wave picks

try
    endtime = oridStruct.(sprintf('eq_%d',eq)).time_phase(intersect(find(strcmp(oridStruct.(sprintf('eq_%d',eq)).phase,'S')),...
        find(strcmp(oridStruct.(sprintf('eq_%d',eq)).sta,stasz))))+datenum(0,0,0,0,0,1);
    w_clean_e = extract(w_clean,'TIME',starttime,endtime);
catch
    w_clean_e = extract(w_clean,'TIME',starttime);
end





%%

for count3 = 1:numel(ch2)
    W(1,3-count3+1) = w_clean_e(count2+count3-1);
    W_all(1,3-count3+1) = w_clean(count2+count3-1);
end
clear count3
%% rotate waveforms
TC=threecomp(W_all,az);
TCR = TC.rotate();
TCRW = TCR.waveform();
TCRP = TCR.particlemotion();
%TCRP.plot2()

for count3 = 1:numel(ch)
    Haystack_data.(ch{count3}) = get(TCRW(count3),'data');
    get(TCRW(count3),'channel')
    needle.(ch{count3}) = Haystack_data.(ch{count3})(P_ind:P_ind+num_samps-1);

    %%
    lngX = length(Haystack_data.(ch{count3}));
    lngY = length(needle.(ch{count3}));
    assert(lngX >= lngY);
    lags = 0:(lngX-lngY);
    %%
    for valz = lags
        c.(ch{count3})(valz+1) = xcorr(Haystack_data.(ch{count3})(valz+1:valz+lngY) -...
            mean(Haystack_data.(ch{count3})(valz+1:valz+lngY)), needle.(ch{count3}) - mean(needle.(ch{count3})),0,'coeff');
    end
end

%% Check and plot directionality and polarity of waveforms
dtac = [Haystack_data.HHZ'; Haystack_data.HHT'; Haystack_data.HHR'];
wndo = 20;%14;%31;
[azim, incd, ellip] = polar_coherency(dtac, wndo);
[azim_cov, incd_cov, ellip_cov] = polar_covariance(dtac, wndo);
phase_plots(Haystack_data,P_ind, S_ind, azim,incd,ellip,eq,stasz)
%phase_plots(Haystack_data,P_ind, S_ind, azim_cov,incd_cov,ellip_cov,eq,stasz)

%%
toc
%end