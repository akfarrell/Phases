% function phase_hunt_determineRealphase_rotated_narrowLOWFreq_loop(eq,stasz,sec)
tic
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase/reexternalre/')

addpath(genpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase'))
addpath('/Users/alexandrafarrell/Desktop/akfarrell/')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Polarizemic-master/functions/')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/fillPage/')
%clear;


clearvars -except eq stasz sec


close all; clc;
%fil=[2 25];
%fil=[10 20];
[evidStruct, allevids] = get_eq_info();
load('siteStruct.mat');
%stylie = 'min'; %%%%% CHANGE!!!---------0----min,max,abz
stylie = {'max','min'};
warning('off')
%ch = {'HHE','HHN','HHZ'};
ch = {'HHZ','HHN','HHE'};
ch2 = {'HHZ','HHR','HHT'};
cutoff_val = 0.6;
S_pad = 25;
SECSPERDAY = 60 * 60 * 24;
num_samps = 31;
%close all
% eq = 2066;
% stasz = 'PLLL';
% eq = 2249;
% stasz = 'PLLL';

% 
eq = 1951;
stasz = 'PLSQ';
sec=2;

% eq = 2501;
% stasz = 'PLMN';


Pvel = 4.2; %km/s
Svel = 2.35; %km/s
%val_try = 569; % If commented out, find val_try (is P_ind or S_ind) ~line
%110
count=find(allevids == eq);
if isempty(count)
    clear evidStruct allevids count
    [evidStruct, allevids] = get_eq_info('/Users/alexandrafarrell/Desktop/dbplutons_alex','true'); %not relocated, with errors
    count = find(allevids == eq);
end

directory = '/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase/wf_objs';
filename = sprintf('wf_%d.mat',allevids(count));
filename_wPath = fullfile(directory,filename);
%if exist(filename_wPath,'file')
%    load(filename_wPath)
%else
%create and clean waveform object

[w_raw,EvS,stations_inEq] = get_wf(allevids(count),evidStruct);
if isempty(find(strcmp(get(w_raw,'channel'),'HHE')))
    error('Stations not aligned correctly')
end
%w_clean_filt = waveform_clean(w_raw, filterobject('b', fil, 2));
w_clean = waveform_clean(w_raw);
%save(filename_wPath,'w_clean', 'EvS', 'stations_inEq');
%end
stationz = get(w_clean,'station');


count2 = min(find(strcmp(stasz,stationz)));%:3:numel(w_clean) %43-45 is PLMN 43:3:43
chan_temp = get(w_clean,'channel');
if isempty(find(strcmp(chan_temp(count2:count2+2),'HHE')))
    error('Stations not aligned correctly')
end
%HHE = count, HHN = count+1, HHZ = count+2

ind_P = intersect(find(strcmp(evidStruct.(EvS).sta,stationz{count2})),find(strcmp(evidStruct.(EvS).phase,'P')));
time_Parr = evidStruct.(EvS).time_phase(ind_P);

try
    ind_S = intersect(find(strcmp(evidStruct.(EvS).sta,stationz{count2})),find(strcmp(evidStruct.(EvS).phase,'S')));
    time_Sarr = evidStruct.(EvS).time_phase(ind_S);
catch
    ind_S = numel(get(w_clean(1),'data'));
end
dnum = zeros(1,numel(get(w_clean(count2),'data')));
dnum(1) = datenum(get(w_clean(count2),'start'));
freq = get(w_clean(count2), 'freq');
for l = 2:numel(get(w_clean(count2),'data'))
    dnum(l) = datenum((l/freq)/SECSPERDAY+dnum(1));
end
P_ind = find(dnum>=time_Parr, 1,'first');
if isempty(P_ind)
    P_ind = 0;
end
S_ind = find(dnum>=time_Sarr, 1,'first');
if isempty(S_ind)
    S_ind = numel(get(w_clean(1),'data'))-num_samps;
end
%%
% fname = sprintf('corr_%d_%s_%s.mat',allevids(count),stationz{count2},stylie{1});
% directory = '/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase/wf_objs';
% directory2 = sprintf('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase/corrs/%d',allevids(count));
% filename_wPath2 = fullfile(directory2,fname);
% load(filename_wPath2)
% c_orig=c;
%P_pad = 15


numz = find(strcmp(siteStruct.sta,stasz));
% find azimuth, distance, and incidence angle
az = azimuth(siteStruct.lat(numz), siteStruct.lon(numz),evidStruct.(sprintf('eq_%d',eq)).lat(1),evidStruct.(sprintf('eq_%d',eq)).lon(1));
starttime = evidStruct.(sprintf('eq_%d',eq)).time_phase(intersect(find(strcmp(evidStruct.(sprintf('eq_%d',eq)).phase,'P')),...
    find(strcmp(evidStruct.(sprintf('eq_%d',eq)).sta,stasz))))-datenum(0,0,0,0,0,1);
dist = distance(siteStruct.lat(numz), siteStruct.lon(numz),evidStruct.(sprintf('eq_%d',eq)).lat(1),evidStruct.(sprintf('eq_%d',eq)).lon(1))*111.12;
inc = atand((siteStruct.elev(numz)+evidStruct.(sprintf('eq_%d',eq)).depth(1))/dist);

az2 = az;
if az > 180
    az = az-180;
else
    az = az+180;
end

load('phaseStruct_backup/pS1.mat')

testvar = phaseStruct;
phase_times = phaseStruct.(EvS).(stasz);
ph_indz = zeros(size(phase_times));
for cnt = 1:numel(phase_times)
    ph_indz(cnt) = find(dnum>=phase_times(cnt), 1,'first');
end
if sec == 1
    val_try = P_ind; indie = 1; ph_used = 'P';
elseif sec == 3
    val_try = S_ind; ph_used = 'S'; indie = 3; %CHANGE!!!!
elseif sec == 4
    val_try = S_ind; ph_used = 'S'; indie = 4;
elseif sec == 2
    try
        vee = 1;
        val_try = ph_indz(vee);
        indie = vee+1;
        if vee == 1
            ph_used = 'P1';
        elseif vee == 2
            ph_used = 'P2';
        elseif vee == 3
            ph_used = 'P3';
        end
    catch
    end
elseif sec == 5
    try
        vee = 2;
        val_try = ph_indz(vee);
        indie = vee+1;
        if vee == 1
            ph_used = 'P1';
        elseif vee == 2
            ph_used = 'P2';
        elseif vee == 3
            ph_used = 'P3';
        end
    catch
    end
end

endtime = P_ind+ceil((dist/Svel-dist/Pvel)*100)+100;
% Find best filter using periodogram and filter to that
fs = get(w_clean(count2),'freq');
w_clean_copy = w_clean;
w_clean_copy = w_clean_copy(count2:count2+2);
compnt = get(w_clean_copy,'channel');
data = get(w_clean_copy,'data');
%close all
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
    
    %CHECK WITH NORMALIZING AND ADDING TOGETHER SPECTRA
    max_p = max(pdgrm.(compnt{countz}));
    pdgrm_norm.(compnt{countz}) = pdgrm.(compnt{countz})./max_p;
    clear max_p
    
    % TRY THIS WITH DETERMINE PEAKS OF BOTH!!!!!
    try
        [k.(compnt{countz}),w.(compnt{countz})] = det_pks(pdgrm.(compnt{countz}),f_pdgrm.(compnt{countz}),3,'pdgrm','low-freq')
    catch
        [k.(compnt{countz}),w.(compnt{countz})] = det_pks(pdgrm.(compnt{countz}),f_pdgrm.(compnt{countz}),3,'pdgrm')
        disp('NOT LOW-FREQ')
    end
    %[k1.(compnt{countz}),w1.(compnt{countz})] = det_pks(pdgrm1.(compnt{countz}),f_pdgrm1.(compnt{countz}),3,'pdgrm')
    %w1.(compnt{countz})=w1.(compnt{countz})/2;
    w.(compnt{countz})=w.(compnt{countz})/2;
    errorz(countz,:) = [k.(compnt{countz})-w.(compnt{countz}),k.(compnt{countz})+w.(compnt{countz})];
    %errorz1.(compnt{countz}) = [k1.(compnt{countz})+w1.(compnt{countz}),k1.(compnt{countz})-w1.(compnt{countz})];
    
end
% Determine peak frequency of normalized, composite spectra

pdgrm_tot = pdgrm_norm.HHZ+pdgrm_norm.HHE+pdgrm_norm.HHN;
[pk_tot,w_tot,which_pk] = det_pks(pdgrm_tot,f_pdgrm.HHZ,3,'pdgrm','low-freq')
if strcmp(stasz,'PLCM') && eq==2165 && indie == 1
    warning('on','all')
    warning('PLCM for eq 2165 chooses 40 Hz - forcing to 14.29')
    pk_tot = 14.29; %PLCM picks dominant frequency at 40 Hz, because the first arrival is shitty
    which_pk = '-';
    warning('off','all')
end

%% Check into values for best fitting frequencies and see how they compare to waveforms
intersect_vals = RangesIntersect(errorz(1,:), errorz(2,:), errorz(3,:));
a = linspace(errorz(1,1),errorz(1,2),50);
b = linspace(errorz(2,1),errorz(2,2),50);
c = linspace(errorz(3,1),errorz(3,2),50);
d = mean(horzcat(a,b,c));

%make sine waves of frequency picked up by det_pks
sig_len = numel(sig_data.HHZ);
x = 0:sig_len-1; %128(127) okay, 256(255) is best
fs = 100;
y = 100*sin(2*pi*d/fs*x);
y2 = 100*sin(2*pi*pk_tot/fs*x);
figure() %figure to show what the proposed frequency looks like comared to the waveforms
plot(x,y)
hold on
plot(x,y2*1000,'b--')
plot(x,sig_data.HHZ,'r')
plot(x,sig_data.HHN,'k')
plot(x,sig_data.HHE,'g')
title(sprintf('d - pkTot = %2.2f',d-pk_tot))
hold off
%%
fprintf('d - pkTot = %2.2f\n',d-pk_tot)
disp('using pk_tot')

%% Use the filter above on the data
fil = [0.8 min(max(25,pk_tot+10),40)];
% fil = [max(pk_tot-10,0.5) pk_tot+10];%[max(pk_tot-5,0.5) pk_tot+5];%[0.5 25];%[pk_tot-2 pk_tot+2];
w_clean_filt = waveform_clean(w_clean, filterobject('b', fil, 2));
w_clean_nofilt = w_clean;
% warning('on','all')
% warning('USING UNFILTERED DATA IN W_CLEAN_FILT')
% warning('off','all')
% w_clean_filt = w_clean_nofilt;
% %This shows the ellipticity, azimuth, and incidence of the ray at that
%filter

%% Subset data to around the P- and S-wave picks

% try
%     endtime = evidStruct.(sprintf('eq_%d',eq)).time_phase(intersect(find(strcmp(evidStruct.(sprintf('eq_%d',eq)).phase,'S')),...
%         find(strcmp(evidStruct.(sprintf('eq_%d',eq)).sta,stasz))))+datenum(0,0,0,0,0,1);
%     w_clean_e = extract(w_clean_filt,'TIME',starttime,endtime);
% catch
%     w_clean_e = extract(w_clean_filt,'TIME',starttime);
% end





%%
for count3 = 1:numel(ch)
    %W(1,3-count3+1) = w_clean_e(count2+count3-1);
    W_all(1,3-count3+1) = w_clean_filt(count2+count3-1);
    %organizing and finding only waveforms for event of interest
    W_all_unfilt(1,3-count3+1) = w_clean_nofilt(count2+count3-1);
    
end
try
    W_all_short = extract(W_all,'TIME', time_Parr-datenum(0,0,0,0,0,2), time_Sarr+datenum(0,0,0,0,0,2));
    W_all_unfilt_short = extract(W_all_unfilt,'TIME', time_Parr-datenum(0,0,0,0,0,2), time_Sarr+datenum(0,0,0,0,0,2));
catch
    W_all_short = extract(W_all,'TIME', time_Parr-datenum(0,0,0,0,0,3), dnum(endtime));
    W_all_unfilt_short = extract(W_all_unfilt,'TIME', time_Parr-datenum(0,0,0,0,0,3), dnum(endtime));
end
clear count3
%% rotate waveforms
wndo = 32;%30;%20;%14;%31;
TC=threecomp(W_all,az2);
TCS = threecomp(W_all_short,az2);
TCR = TC.rotate();
TCSR = TCS.rotate();
TCW = TCR.waveform();
TCRW = TCSR.waveform();


TCUS = threecomp(W_all_unfilt_short,az2); %az2
TCUSR = TCUS.rotate();





wind_use = max(wndo,ceil(200/pk_tot));
fprintf('Using window length %2.4f\n', wind_use);
[TCP,vals,vecs] = TCR.particlemotion(0.01,wind_use/100);
[TCP_overall,vals_overall,vecs_overall] = TCR.particlemotion();
%[TCP,vals,vecs] = TCR.particlemotion();
[lambda, I]=sort(vals,1,'descend');
X = vecs(:,I(1),:); %get greatest eigenvector
TCSRP = TCSR.particlemotion(0.01,wind_use/100);
TCSRPL = TCSR.particlemotion();
TCUSRP = TCUSR.particlemotion(0.01,wind_use/100);
TCUSRPL = TCUSR.particlemotion();

TC_unfilt = threecomp(W_all_unfilt,az);
%TCR_unfilt = TC_unfilt.rotate();
TCW_unfilt = TC_unfilt.waveform();
[TCP_unfilt, vals_unfilt, vecs_unfilt] = TC_unfilt.particlemotion(0.01,wind_use/100);

clear c c2
cutoff_val = 0.8;
data_for_xcorr = get(TCW,'data');

for count3 = 1:numel(ch)
    Haystack_data.(ch{count3}) = get(W_all(count3),'data');
    Haystack_data_unfilt.(ch{count3}) = get(W_all_unfilt(count3),'data');
    Haystack_data_rotate.(ch2{count3}) = data_for_xcorr(1,count3);
    Haystack_data_rotate.(ch2{count3}) = Haystack_data_rotate.(ch2{count3}){1};
    get(W_all(count3),'channel')
    needle.(ch{count3}) = Haystack_data.(ch{count3})(P_ind:P_ind+num_samps-1);
    needle_r.(ch2{count3}) = Haystack_data_rotate.(ch2{count3})(P_ind:P_ind+num_samps-1);
    needle_r2.(ch2{count3}) = Haystack_data_rotate.(ch2{count3})(S_ind:S_ind+num_samps-1);
    
    lngX = length(Haystack_data_rotate.(ch2{count3}));
    lngY = length(needle_r.(ch2{count3}));
    assert(lngX >= lngY);
    lags = 0:(lngX-lngY);
    for valz = lags
        c.(ch2{count3})(valz+1) = xcorr(Haystack_data_rotate.(ch2{count3})(valz+1:valz+lngY) - ...
            mean(Haystack_data_rotate.(ch2{count3})(valz+1:valz+lngY)), needle_r.(ch2{count3}) - mean(needle_r.(ch2{count3})),0,'coeff');
        c2.(ch2{count3})(valz+1) = xcorr(Haystack_data_rotate.(ch2{count3})(valz+1:valz+lngY) - ...
            mean(Haystack_data_rotate.(ch2{count3})(valz+1:valz+lngY)), needle_r2.(ch2{count3}) - mean(needle_r2.(ch2{count3})),0,'coeff');
    end
    
end
i.HHZ = val_try; i.HHR = val_try; i.HHT = val_try;
plot_xcorrs(lags, ch2, Haystack_data_rotate, needle_r, c, stationz{count2}, cutoff_val, i, allevids(count), 4,P_ind,S_ind)
plot_xcorrs(lags, ch2, Haystack_data_rotate, needle_r2, c2, stationz{count2}, cutoff_val, i, allevids(count), 4,P_ind,S_ind)

%% Check and plot directionality and polarity of waveforms
dtac_i = [Haystack_data.(ch{1})'; Haystack_data.(ch{3})'; Haystack_data.(ch{2})'];
dtac_unfilt = [Haystack_data_unfilt.(ch{1})'; Haystack_data_unfilt.(ch{2})'; Haystack_data_unfilt.(ch{3})'];
dtac = dtac_i(:,max(P_ind-100,1):S_ind+100);
% dtac = dtac_i;
sps = 100;
% window length in (s)
windt = wndo/100;
% fraction of window for subwindows used in cross-spectrum estimation 
fsw = 2/wndo;
% fraction of window for overlap between time bins
fov = 0.9;

% time-frequency polarization analysis
[azim, incd, ellip, rpol, ppol, tvecC_zs, F, Cxypzz, Cxypnn, Cxypee, largest_eig] = ...
    polartf_cross_spectrum(dtac,windt,sps,fsw,fov);
largest_eig;


[aa bb] = size(Cxypzz);
pctr = 0; %0.002;
%mxc = max(max(Cxypzz));
mxc = max(max([Cxypzz Cxypee Cxypnn]));
for ii=1:aa
    for jj=1:bb
        if ((Cxypzz(ii,jj)> pctr*mxc) | (Cxypnn(ii,jj)> pctr*mxc) | (Cxypee(ii,jj)> pctr*mxc))
            incdm(ii,jj) = incd(ii,jj);
            azimm(ii,jj) = azim(ii,jj);
            ellipm(ii,jj) = ellip(ii,jj);
        else
            incdm(ii,jj) = NaN;
            azimm(ii,jj) = NaN;
            ellipm(ii,jj) = NaN;
        end        
    end
end

% plot polarization attributes
figure
fsize = 12;
% colormap to show NaN as blank
mycolormap = [ ones(1,3); jet(300)];
tvecC_z = tvecC_zs+(P_ind-100)/100;
[~,min_ind] = min(abs(tvecC_z.*100-val_try));
largest_eig(:,min_ind)
%
% subpanel 1
%
subplot(5,1,2)
imagesc(tvecC_z,F,incdm,[-1 180]); axis xy; axis([tvecC_z(1) tvecC_z(end) 0 30]); colormap(mycolormap)
ylabel('Frequency (Hz)')
hh = colorbar('EastOutside','FontSize',fsize,'Ytick',[0:90:180],'Ylim',[0 180]);
label = sprintf(' Incidence (deg) ');
set(get(hh,'YLabel'),'String',label,'FontSize',fsize)
set(gca,'Fontsize',fsize);
%
% subpanel 2
%
subplot(5,1,1)
imagesc(tvecC_z,F,azimm,[-1 180]); axis xy; axis([tvecC_z(1) tvecC_z(end) 0 30]); colormap(mycolormap)
ylabel('Frequency (Hz)')
hh = colorbar('EastOutside','FontSize',fsize,'Ytick',[0:90:180],'Ylim',[0 180]);
label = sprintf(' Azimuth (deg) ');
set(get(hh,'YLabel'),'String',label,'FontSize',fsize)
set(gca,'Fontsize',fsize);
%
% subpanel 3
%
subplot(5,1,3)
imagesc(tvecC_z,F,ellipm,[-.01 1]); axis xy; axis([tvecC_z(1) tvecC_z(end) 0 30]); colormap(mycolormap)
ylabel('Frequency (Hz)')
%ylabel('Frequency (Hz)')
hh = colorbar('EastOutside','FontSize',fsize,'Ytick',[0:0.5:1],'Ylim',[0 1]);
label = sprintf(' Ellipticity ');
set(get(hh,'YLabel'),'String',label,'FontSize',fsize)
set(gca,'Fontsize',fsize);
%
% subpanel 4
%
subplot(5,1,4)
plot(([1:length(dtac(1,:))]*(1/sps))+(max(P_ind-100,1))/100,dtac(1,:)*(10^6),'b-'); ...
    axis([tvecC_z(1) tvecC_z(end) ...
    -max(abs(dtac(1,:)*(10^6))) max(abs(dtac(1,:)*(10^6)))]); 
colorbar; colormap(mycolormap)
colorbar('Visible','off');
ylabel('Particle Velocity (\mum/s)')
%xlabel('Time (s)')
set(gca,'Fontsize',fsize);
%
%
subplot(5,1,5)
hurr = floor(min(min(10*log10(abs(Cxypzz/(10^-18))+eps))));
durr = ceil(max(max(10*log10(abs(Cxypzz/(10^-18))+eps))));
imagesc(tvecC_z,F,10*log10(abs(Cxypzz/(10^-18))+eps),[hurr durr]); axis xy; axis([tvecC_z(1) tvecC_z(end) 0 30]); colormap(mycolormap)
%caxis([30 90])
ylabel('Frequency (Hz)')
xlabel('Time (s)')
%ylabel('Frequency (Hz)')
hh = colorbar('EastOutside','FontSize',fsize,'Ytick',[hurr:20:durr],'Ylim',[hurr durr]);
label = sprintf(' dB rel. 1 nm/s ');
set(get(hh,'YLabel'),'String',label,'FontSize',fsize)
set(gca,'Fontsize',fsize);

orient portrait
%print(gcf,'-dpsc','-fillpage','-r300','-opengl','cleve.ps');
print(gcf,'-dpdf','-fillpage','-r300','-opengl','Utu.pdf');

% [azim, incd, ellip] = polar_coherency(dtac, wndo);
% 
% [azim_u, incd_u, ellip_u] = polar_coherency(dtac_unfilt,wndo);
% for cz = 1:numel(azim)
%     if azim(cz)<0
%         azim(cz) = -azim(cz)+180;
%     end
%     if azim_u(cz)<0
%         azim_u(cz) = -azim_u(cz)+180;
%     end
%     %     if incd(cz)>90
%     %         incd(cz)=90-(incd(cz)-90);
%     %     end
% end



%% Get phase information if availabe
%load('phaseStruct.mat');

% try
%     phase_plots(Haystack_data_unfilt,P_ind, S_ind, azim_adj,incd,ellip,eq,stasz,val_try,endtime,az_adj,inc,ph_indz)
%     phase_plots(Haystack_data,P_ind, S_ind, azim_adj,incd,ellip,eq,stasz,val_try,endtime,az_adj,inc,ph_indz)
% catch
%     phase_plots(Haystack_data_unfilt,P_ind, S_ind, azim_adj,incd,ellip,eq,stasz,val_try,endtime,az_adj,inc)
%     phase_plots(Haystack_data,P_ind, S_ind, azim_adj,incd,ellip,eq,stasz,val_try,endtime,az_adj,inc)
% end

%% Plot rectilinearity, planarity, azimuth, energy, and inclination
tcrp_az = get(TCP,'azimuth');
tcrp_inc = get(TCP,'inclination');
tcrp_rec = get(TCP,'rectilinearity');
tcrp_plan = get(TCP,'planarity');
tcrp_en = get(TCP,'energy');

tcrp_az_u = get(TCP_unfilt,'azimuth');
tcrp_inc_u = get(TCP_unfilt,'inclination');
tcrp_rec_u = get(TCP_unfilt,'rectilinearity');
tcrp_plan_u = get(TCP_unfilt,'planarity');
tcrp_en_u = get(TCP_unfilt,'energy');
%phase_plots(Haystack_data, P_ind, S_ind, get(tcrp_az,'data'),get(tcrp_inc,'data'),get(tcrp_rec,'data'),eq,stasz);

%% Plots using 0-360 for azimuth
% plot_threecomp_directionality(endtime,val_try,Haystack_data_unfilt,P_ind,S_ind,tcrp_az,tcrp_inc, tcrp_rec, tcrp_plan, tcrp_en,az,inc)
% plot_threecomp_directionality(endtime,val_try,Haystack_data,P_ind,S_ind,tcrp_az,tcrp_inc, tcrp_rec, tcrp_plan, tcrp_en,az,inc)
%% Plots using 240-500 for azimuth
tcrp_az_data = get(tcrp_az,'data');
% tcrp_az_data_adj = zeros(size(tcrp_az_data));
% for ct = 1:numel(tcrp_az_data)
%     if tcrp_az_data(ct) < min_num
%         tcrp_az_data_adj(ct) = tcrp_az_data(ct)+360;
%     else
%         tcrp_az_data_adj(ct) = tcrp_az_data(ct);
%     end
% end

% try
%With phase arrivals in red
% plot_threecomp_directionality(endtime,val_try,Haystack_data_unfilt,P_ind,...
%     S_ind,tcrp_az,tcrp_inc, tcrp_rec, tcrp_plan, tcrp_en,az,inc,tcrp_az_data_adj,ph_indz)
plot_threecomp_directionality(endtime,val_try,Haystack_data,P_ind,S_ind,...
    tcrp_az,tcrp_inc, tcrp_rec, tcrp_plan, tcrp_en,az,inc,tcrp_az_data,ph_indz)
%h = plotpm(TCP,'true')
%f = plotpm(TCP_overall,'true')

f = plotpm(TCSRPL,'true');
h = plotpm(TCSRP,'true');

% if ~exist('short_window', 'dir')
%     mkdir short_window
% end
% if ~exist('long_window','dir')
%     mkdir long_window
% end
for countiesz=1:2
    if countiesz == 1
        directory = sprintf('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase/short_window/unfilt');
        filename = sprintf('%d_%s_%d_%s.pdf',eq,stasz,wind_use,ph_used);
        filename_wPath = fullfile(directory,filename);
        %hgexport(h, filename_wPath, hgexport('factorystyle'), 'Format', 'pdf');
        fp = fillPage(h, 'margins', [0 0 0 0], 'papersize', [11 8.5]);
        print(h, '-dpdf', '-r400',filename_wPath);
        clear directory filename filename_wPath fp
    elseif countiesz == 2 && (val_try == P_ind || val_try == S_ind)
        directory = sprintf('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase/long_window/unfilt');
        filename = sprintf('%d_%s_%d_%s.pdf',eq,stasz,wind_use,ph_used);
        filename_wPath = fullfile(directory,filename);
        %hgexport(h, filename_wPath, hgexport('factorystyle'), 'Format', 'pdf');
        fp = fillPage(f, 'margins', [0 0 0 0], 'papersize', [11 8.5]);
        print(f, '-dpdf', '-r400',filename_wPath);
        clear directory filename filename_wPath fp
    end
end

%% make comparison figure%% make phase plot data range from 240-500, then plot
%[azim_cov, incd_cov, ellip_cov] = polar_covariance(dtac, wndo);
% ###make whole figure with all params and wfs
[min_diff,id_diff] = min(abs(F-pk_tot));
min_num = 240;
for count = 1:numel(tcrp_az_data)
    if tcrp_az_data(count)>180
        tcrp_az_data_compressed(count) = tcrp_az_data(count)-180;
    else
         tcrp_az_data_compressed(count) = tcrp_az_data(count);
    end
end
tcrp_az_use = tcrp_az_data_compressed(max(1,P_ind-100):S_ind+100);
figure()
% for count2 = 3:9
%     
%     plot(tvecC_z,incdm(count2,:))
%     hold on
% end
plot(tvecC_z,azim(id_diff,:))
hold on
xs = linspace(tvecC_z(1),tvecC_z(end),numel(tcrp_az_use));
plot(xs,tcrp_az_use,'r')
%%
plot_haney_directionality(endtime,val_try,Haystack_data,P_ind,S_ind,...
     azim(id_diff,:),incd(id_diff,:), ellip(id_diff,:), rpol(id_diff,:), ppol(id_diff,:),az,inc,needle,tvecC_z)
 
plot_both_directionality(endtime,val_try, Haystack_data,P_ind,S_ind,...
    azim(id_diff,:),incd(id_diff,:), ellip(id_diff,:), rpol(id_diff,:), ppol(id_diff,:),az, tvecC_z,inc,needle,tcrp_az,tcrp_inc,tcrp_rec, tcrp_plan,TC)

%%
% if az < min_num
%     az_adj = az+360;
% else
%     az_adj = az;
% end
% 
% azim_adj = zeros(size(azim));
% for ct = 1:numel(azim)
%     if azim(ct) < min_num
%         azim_adj(ct) = azim(ct)+360;
%     else
%         azim_adj(ct) = azim(ct);
%     end
% end



rec_vals = get(tcrp_rec,'data');
plan_vals = get(tcrp_plan,'data');
[valrec,indrec] = max(rec_vals(val_try:val_try+num_samps));
if plan_vals(indrec+val_try-1)<0.5
    disp('Low degree of planarity!!')
end
tcrp_inc_data = get(tcrp_inc,'data');
tcrp_inc_data = 90-tcrp_inc_data;

az_calc = tcrp_az_data(indrec+val_try-1);
inc_calc = tcrp_inc_data(indrec+val_try-1);
clear phaseStruct
load('phaseStruct.mat')
phaseStruct.(sprintf('eq_%d',eq)).(stasz).incWest(indie) = mean(tcrp_inc_data(val_try:val_try+6));
phaseStruct.(sprintf('eq_%d',eq)).(stasz).azWest(indie) = mean(tcrp_az_data(val_try:val_try+6));
phaseStruct.(sprintf('eq_%d',eq)).(stasz).peak{indie} = which_pk;
phaseStruct.(sprintf('eq_%d',eq)).(stasz).planWest(indie) = mean(plan_vals(val_try:val_try+6));
phaseStruct.(sprintf('eq_%d',eq)).(stasz).recWest(indie) = mean(rec_vals(val_try:val_try+6));
phaseStruct.(sprintf('eq_%d',eq)).(stasz).ind(indie) = val_try;
phaseStruct.(sprintf('eq_%d',eq)).(stasz).xcorrP_z(indie) = c.HHZ(val_try);
phaseStruct.(sprintf('eq_%d',eq)).(stasz).xcorrP_r(indie) = c.HHR(val_try);
phaseStruct.(sprintf('eq_%d',eq)).(stasz).xcorrP_t(indie) = c.HHT(val_try);

phaseStruct.(sprintf('eq_%d',eq)).(stasz).xcorrS_z(indie) = c2.HHZ(val_try);
phaseStruct.(sprintf('eq_%d',eq)).(stasz).xcorrS_r(indie) = c2.HHR(val_try);
phaseStruct.(sprintf('eq_%d',eq)).(stasz).xcorrS_t(indie) = c2.HHT(val_try);

if sec == 2
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).xcorr_vals_P_z = c.HHZ(P_ind:S_ind);
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).xcorr_vals_P_r = c.HHR(P_ind:S_ind);
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).xcorr_vals_P_t = c.HHZ(P_ind:S_ind);
    
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).xcorr_vals_S_z = c2.HHZ(P_ind:S_ind);
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).xcorr_vals_S_r = c2.HHR(P_ind:S_ind);
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).xcorr_vals_S_t = c2.HHT(P_ind:S_ind);
end


if sec==1
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).largest_eigs = largest_eig(:,min_ind);
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).haney_frequencies = F;
end

phaseStruct.(sprintf('eq_%d',eq)).(stasz).std_incWest(indie) = std(tcrp_inc_data(val_try:val_try+6));
phaseStruct.(sprintf('eq_%d',eq)).(stasz).std_azWest(indie) = std(tcrp_az_data(val_try:val_try+6));
phaseStruct.(sprintf('eq_%d',eq)).(stasz).std_planWest(indie) = std(plan_vals(val_try:val_try+6));
phaseStruct.(sprintf('eq_%d',eq)).(stasz).std_recWest(indie) = std(rec_vals(val_try:val_try+6));

try
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).haney_az(indie) = mean(azim(id_diff,((val_try):(val_try+6))./100));
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).haney_inc(indie) = mean(incd(id_diff,((val_try):(val_try+6))./100));
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).haney_ellip(indie) = mean(ellip(id_diff,((val_try):(val_try+6))./100));
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).haney_rect(indie) = mean(rpol(id_diff,((val_try):(val_try+6))./100));
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).haney_plan(indie) = mean(ppol(id_diff,((val_try):(val_try+6))./100));
    
    
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).std_haney_az(indie) = std(azim(id_diff,((val_try):(val_try+6))./100));
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).std_haney_inc(indie) = std(incd(id_diff,((val_try):(val_try+6))./100));
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).std_haney_ellip(indie) = std(ellip(id_diff,((val_try):(val_try+6))./100));
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).std_haney_rect(indie) = std(rpol(id_diff,((val_try):(val_try+6))./100));
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).std_haney_plan(indie) = std(ppol(id_diff,((val_try):(val_try+6))./100));
    
catch
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).haney_az(indie) = mean(interp1(tvecC_z,azim(id_diff,:),((val_try):(val_try+6))./100));
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).haney_inc(indie) = mean(interp1(tvecC_z,incd(id_diff,:),((val_try):(val_try+6))./100));
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).haney_ellip(indie) = mean(interp1(tvecC_z,ellip(id_diff,:),((val_try):(val_try+6))./100));
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).haney_rect(indie) = mean(interp1(tvecC_z,rpol(id_diff,:),((val_try):(val_try+6))./100));
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).haney_plan(indie) = mean(interp1(tvecC_z,ppol(id_diff,:),((val_try):(val_try+6))./100));
    
    
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).std_haney_az(indie) = std(interp1(tvecC_z,azim(id_diff,:),((val_try):(val_try+6))./100));
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).std_haney_inc(indie) = std(interp1(tvecC_z,incd(id_diff,:),((val_try):(val_try+6))./100));
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).std_haney_ellip(indie) = std(interp1(tvecC_z,ellip(id_diff,:),((val_try):(val_try+6))./100));
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).std_haney_rect(indie) = std(interp1(tvecC_z,rpol(id_diff,:),((val_try):(val_try+6))./100));
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).std_haney_plan(indie) = std(interp1(tvecC_z,ppol(id_diff,:),((val_try):(val_try+6))./100));
end


% phaseStruct.(sprintf('eq_%d',eq)).(stasz).az_W_calc(indie) = tcrp_az_data(indrec+val_try-1);
% phaseStruct.(sprintf('eq_%d',eq)).(stasz).inc_W_calc(indie) = inc_calc;
% phaseStruct.(sprintf('eq_%d',eq)).(stasz).ind_calc(indie) = indrec+val_try-1;
% phaseStruct.(sprintf('eq_%d',eq)).(stasz).rec_W_calc(indie) = valrec;
% phaseStruct.(sprintf('eq_%d',eq)).(stasz).plan_W_calc(indie) = plan_vals(indrec+val_try-1);
phaseStruct.(sprintf('eq_%d',eq)).(stasz).freq(indie) = pk_tot;
phaseStruct.(sprintf('eq_%d',eq)).(stasz).number_of_cycles(indie) = pk_tot*(wind_use/100);
phaseStruct.(sprintf('eq_%d',eq)).(stasz).window(indie) = wind_use;

%%
if ~isfield(phaseStruct.(sprintf('eq_%d',eq)).(stasz),'phase')
    if numel(testvar.(sprintf('eq_%d',eq)).(stasz)) == 3 && exist('time_Parr','var')
        if exist('time_Sarr','var') && ~isempty(time_Sarr)
            phaseStruct.(sprintf('eq_%d',eq)).(stasz).phase = {'P','P1','P2','P3','S'};
        else
            phaseStruct.(sprintf('eq_%d',eq)).(stasz).phase = {'P','P1','P2','P3'};
        end
    elseif numel(testvar.(sprintf('eq_%d',eq)).(stasz)) == 3 && ~exist('time_Parr','var')
        if exist('time_Sarr','var') && ~isempty(time_Sarr)
            phaseStruct.(sprintf('eq_%d',eq)).(stasz).phase = {'P1','P2','P3','S'};
        else
            phaseStruct.(sprintf('eq_%d',eq)).(stasz).phase = {'P1','P2','P3'};
        end
    elseif numel(testvar.(sprintf('eq_%d',eq)).(stasz)) == 2 && exist('time_Parr','var')
        if exist('time_Sarr','var') && ~isempty(time_Sarr)
            phaseStruct.(sprintf('eq_%d',eq)).(stasz).phase = {'P','P1','P2','S'};
        else
            phaseStruct.(sprintf('eq_%d',eq)).(stasz).phase = {'P','P1','P2'};
        end
    elseif numel(testvar.(sprintf('eq_%d',eq)).(stasz)) == 2 && ~exist('time_Parr','var')
        if exist('time_Sarr','var') && ~isempty(time_Sarr)
            phaseStruct.(sprintf('eq_%d',eq)).(stasz).phase = {'P1','P2','S'};
        else
            phaseStruct.(sprintf('eq_%d',eq)).(stasz).phase = {'P1','P2'};
        end
    elseif numel(testvar.(sprintf('eq_%d',eq)).(stasz)) == 1 && exist('time_Parr','var')
        if exist('time_Sarr','var') && ~isempty(time_Sarr)
            phaseStruct.(sprintf('eq_%d',eq)).(stasz).phase = {'P','P1','S'};
        else
            phaseStruct.(sprintf('eq_%d',eq)).(stasz).phase = {'P','P1'};
        end
    elseif numel(testvar.(sprintf('eq_%d',eq)).(stasz)) == 1 && ~exist('time_Parr','var')
        if exist('time_Sarr','var') && ~isempty(time_Sarr)
            phaseStruct.(sprintf('eq_%d',eq)).(stasz).phase = {'P1','S'};
        else
            phaseStruct.(sprintf('eq_%d',eq)).(stasz).phase = {'P1'};
        end
    end
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).time = evidStruct.eq_2165.time_origin(1);
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).az_exp_adj = az;
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).inc_exp = inc;
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).P_arr = time_Parr;
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).phase_arr = testvar.(sprintf('eq_%d',eq)).(stasz);
    if exist('time_Sarr','var') && ~isempty(time_Sarr)
        phaseStruct.(sprintf('eq_%d',eq)).(stasz).S_arr = time_Sarr;
    end
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).dist = dist;
    phaseStruct.(sprintf('eq_%d',eq)).(stasz).az = az2;
end
save('phaseStruct.mat','phaseStruct')
fil


% catch
%     plot_threecomp_directionality(endtime,val_try,Haystack_data_unfilt,P_ind,...
%         S_ind,tcrp_az,tcrp_inc, tcrp_rec, tcrp_plan, tcrp_en,az_adj,inc,tcrp_az_data_adj)
%     plot_threecomp_directionality(endtime,val_try,Haystack_data,P_ind,S_ind,...
%         tcrp_az,tcrp_inc, tcrp_rec, tcrp_plan, tcrp_en,az_adj,inc,tcrp_az_data_adj)
% end

%%
%plot_threecomp_directionality(endtime,val_try,Haystack_data_unfilt,P_ind,S_ind,tcrp_az_u,tcrp_inc_u, tcrp_rec_u, tcrp_plan_u, tcrp_en_u,az,inc)
% %% Plot record section
% w_clean_HHE= w_clean(1:3:end);
% w_clean_HHN= w_clean(2:3:end);
% w_clean_HHZ= w_clean(3:3:end);
%
% w_prep = prep_for_rs(w_clean,sprintf('eq_%d',eq),siteStruct,evidStruct);
%
% w_prep_E = prep_for_rs(w_clean_HHE,sprintf('eq_%d',eq),siteStruct,evidStruct);
% w_prep_N = prep_for_rs(w_clean_HHN,sprintf('eq_%d',eq),siteStruct,evidStruct);
% w_prep_Z = prep_for_rs(w_clean_HHZ,sprintf('eq_%d',eq),siteStruct,evidStruct);
%
% plotw_rs_test(w_prep)
%
% plotw_rs_test(w_prep_E)
% plotw_rs_test(w_prep_N)
% plotw_rs_test(w_prep_Z)

% %% Plot az, inc, and polarity/ellipticity together on same plot
% figure()
% pad = 100;
% subplot(2,1,1)
% plot(tcrp_az_data_adj,'r')
% hold on
% plot(azim_adj,'b')
% title('Azimuth (degrees clockwise from North)')
% line([0 numel(tcrp_az_data_adj)],[az_adj az_adj],'Color','g','LineStyle',':','LineWidth',1)
% if S_ind > endtime
%     xlim([P_ind-pad, xmax])
% else
%     xlim([P_ind-pad, S_ind+pad])
% end
% ylim([min(min(tcrp_az_data_adj,azim_adj')), max(max(tcrp_az_data_adj,azim_adj'))])
% line([P_ind P_ind], [min(min(tcrp_az_data_adj,azim_adj')), max(max(tcrp_az_data_adj,azim_adj'))],...
%         'Color','m','LineStyle',':','LineWidth',2)
% try
%     line([S_ind S_ind], [min(min(tcrp_az_data_adj,azim_adj')), max(max(tcrp_az_data_adj,azim_adj'))],...
%         'Color','m','LineStyle',':','LineWidth',2)
% catch
% end
%
% subplot(2,1,2)
% plot(get(tcrp_inc,'data'),'r')
% hold on
% plot(incd,'b')
% tcrp_inc_data = get(tcrp_inc,'data');
% title('Incidence Angle (degrees from vertical)')
% line([0 numel(tcrp_inc_data)],[inc inc],'Color','g','LineStyle',':','LineWidth',1)
% if S_ind > endtime
%     xlim([P_ind-pad, xmax])
% else
%     xlim([P_ind-pad, S_ind+pad])
% end
% ylim([min(min(tcrp_inc_data,incd')), max(max(tcrp_inc_data,incd'))])
% line([P_ind P_ind], [min(min(tcrp_inc_data,incd')), max(max(tcrp_inc_data,incd'))],...
%         'Color','m','LineStyle',':','LineWidth',2)
% try
%     line([S_ind S_ind], [min(min(tcrp_inc_data,incd')), max(max(tcrp_inc_data,incd'))],...
%         'Color','m','LineStyle',':','LineWidth',2)
% catch
% end
%
% %% Plot all 3 parts of particle motion on same plot
% figure()
% pad = 100;
%
% %plot(tcrp_az_data_adj,'r')
% tcrp_plan_data = get(tcrp_plan,'data');
% plot(tcrp_plan_data*100,'b')
% hold on
% [hAx,hLine1,hLine2] = plotyy(1:numel(tcrp_az_data_adj),tcrp_az_data_adj,1:numel(tcrp_inc_data),tcrp_inc_data);
%
% % plot(tcrp_inc_data,'g')
% %
%
% set(hLine1,'Color','r')
% set(hLine2,'Color','g')
% legend('Plan','Az','Inc','Location','Best')
% set(hAx(1),'Ylim',[min(min(tcrp_az_data_adj,tcrp_plan_data*100)),max(max(tcrp_az_data_adj,tcrp_plan_data*100))])
% set(hAx(2),'Ylim',[min(tcrp_inc_data),max(tcrp_inc_data)])
% line([0 numel(tcrp_az_data_adj)],[az_adj az_adj],'Color','r','LineStyle',':','LineWidth',1)
% if S_ind > endtime
%     set(hAx(2),'Xlim', [P_ind-pad, xmax]);
%     xlim([P_ind-pad, xmax])
% else
%     set(hAx(2),'Xlim', [P_ind-pad, S_ind+pad]);
%     xlim([P_ind-pad, S_ind+pad])
% end
%
% line([P_ind P_ind], [min(tcrp_plan_data), max(tcrp_az_data_adj)],...
%         'Color','m','LineStyle',':','LineWidth',2)
% try
%     line([S_ind S_ind], [min(tcrp_plan_data), max(tcrp_az_data_adj)],...
%         'Color','m','LineStyle',':','LineWidth',2)
% catch
% end
% hold(hAx(2),'on')
% plot(hAx(2),[0 numel(tcrp_inc_data)],[inc inc],'Color','g','LineStyle',':','LineWidth',1)
% %legend('Plan','Az','Inc','Location','Best','BF Az','P','S','BF Inc')
%%
plot3_edit(TCP,val_try/100,(val_try+15)/100)
plot3_edit(TCP,P_ind/100,(P_ind+15)/100)
plot3_edit(TCP,S_ind/100,(S_ind+15)/100)

azi = mean(tcrp_az_data(val_try:val_try+6));
if azi < 180
    azi = azi+180;
end
azp = phaseStruct.(sprintf('eq_%d',eq)).(stasz).azWest(1);
azs = phaseStruct.(sprintf('eq_%d',eq)).(stasz).azWest(2);
if azp<180
    azp = azp+180;
end
if azs<180
    azs = azs+180;
end
TC_edit=threecomp(W_all,azi);
TCR_edit = TC_edit.rotate();

TC_P=threecomp(W_all,azp);
TCR_P = TC_P.rotate();

TC_S=threecomp(W_all,azs);
TCR_S = TC_S.rotate();


plot3_edit(TCR_edit,val_try/100,(val_try+15)/100)
plot3_edit(TCR_P,P_ind/100,(P_ind+15)/100)
plot3_edit(TCR_S,S_ind/100,(S_ind+15)/100)
%%
toc
%end