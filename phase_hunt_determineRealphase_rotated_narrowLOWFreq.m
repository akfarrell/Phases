% function high_corrs = phase_hunt_changingNeedle_mins(allorids,oridStruct,P,fil,cutoff_val)
% find(strcmp(fieldnames(oridStruct),'eq_whatever')) %to find eq
%orid2166 is origin 203
%orid2034 is origin 128
%orid2277 is origin 277
tic
addpath(genpath('/home/a/akfarrell/Uturuncu/Phase'))
addpath('/raid/home/a/akfarrell/')
addpath('/raid/home/a/akfarrell/Uturuncu')
addpath('/raid/home/a/akfarrell/Polarizemic-master/functions/')
%clear;
clear all; close all; clc;
%fil=[2 25];
%fil=[10 20];
[oridStruct, allorids] = get_eq_info(); 
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
eq = 2166;
stasz = 'PLMD';
Pvel = 4.2; %km/s
Svel = 2.35; %km/s
%val_try = 569; % If commented out, find val_try (is P_ind or S_ind) ~line
%110


count=find(allorids == eq);
directory = '/home/a/akfarrell/Uturuncu/Phase/wf_objs';
filename = sprintf('wf_%d.mat',allorids(count));
filename_wPath = fullfile(directory,filename);
%if exist(filename_wPath,'file')
%    load(filename_wPath)
%else
    %create and clean waveform object
    [w_raw,OrS,stations_inEq] = get_wf(allorids(count),oridStruct);
    %w_clean_filt = waveform_clean(w_raw, filterobject('b', fil, 2));
    w_clean = waveform_clean(w_raw);
    %save(filename_wPath,'w_clean', 'OrS', 'stations_inEq');
%end
stationz = get(w_clean,'station');


count2 = min(find(strcmp(stasz,stationz)));%:3:numel(w_clean) %43-45 is PLMN 43:3:43
%HHE = count, HHN = count+1, HHZ = count+2

ind_P = intersect(find(strcmp(oridStruct.(OrS).sta,stationz{count2})),find(strcmp(oridStruct.(OrS).phase,'P')));
time_Parr = oridStruct.(OrS).time_phase(ind_P);

try
    ind_S = intersect(find(strcmp(oridStruct.(OrS).sta,stationz{count2})),find(strcmp(oridStruct.(OrS).phase,'S')));
    time_Sarr = oridStruct.(OrS).time_phase(ind_S);
catch
    ind_S = numel(get(w_clean(1),'data'));
end
dnum = zeros(1,numel(get(w_clean(count2),'data')));
dnum(1) = datenum(get(w_clean(count2),'start'));
freq = get(w_clean(count2), 'freq');
for l = 2:numel(get(w_clean(count2),'data'))
    dnum(l) = datenum((l/freq)/SECSPERDAY+dnum(1));
end
try
    P_ind = find(dnum>=time_Parr, 1,'first');
catch
    P_ind = 0;
end
try
    S_ind = find(dnum>=time_Sarr, 1,'first');
catch
    S_ind = numel(get(w_clean(1),'data'))-num_samps;
end

% fname = sprintf('corr_%d_%s_%s.mat',allorids(count),stationz{count2},stylie{1});
% directory = '/home/a/akfarrell/Uturuncu/Phase/wf_objs';
% directory2 = sprintf('/home/a/akfarrell/Uturuncu/Phase/corrs/%d',allorids(count));
% filename_wPath2 = fullfile(directory2,fname);
% load(filename_wPath2)
% c_orig=c;
%P_pad = 15


numz = find(strcmp(siteStruct.sta,stasz));
% find azimuth, distance, and incidence angle
az = azimuth(siteStruct.lat(numz), siteStruct.lon(numz),oridStruct.(sprintf('eq_%d',eq)).lat(1),oridStruct.(sprintf('eq_%d',eq)).lon(1));
starttime = oridStruct.(sprintf('eq_%d',eq)).time_phase(intersect(find(strcmp(oridStruct.(sprintf('eq_%d',eq)).phase,'P')),...
    find(strcmp(oridStruct.(sprintf('eq_%d',eq)).sta,stasz))))-datenum(0,0,0,0,0,1);
dist = distance(siteStruct.lat(numz), siteStruct.lon(numz),oridStruct.(sprintf('eq_%d',eq)).lat(1),oridStruct.(sprintf('eq_%d',eq)).lon(1))*111.12;
inc = atand((siteStruct.elev(numz)+oridStruct.(sprintf('eq_%d',eq)).depth(1))/dist);

az2 = az;
if az > 180
    az = az-180;
else
    az = az+180;
end
load('phaseStruct_backup/pS1.mat')
phase_times = phaseStruct.eq_2166.(stasz);
ph_indz = zeros(size(phase_times));
for cnt = 1:numel(phase_times)
    ph_indz(cnt) = find(dnum>=phase_times(cnt), 1,'first');
end
% val_try = P_ind; indie = 1;
% val_try = S_ind; indie = 4; %CHANGE!!!!
try
    vee = 1;
    val_try = ph_indz(vee);
    indie = vee+1;
catch
end
endtime = P_ind+ceil((dist/Svel-dist/Pvel)*100)+100;
% Find best filter using periodogram and filter to that
fs = get(w_clean(count2),'freq');
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
[pk_tot,w_tot] = det_pks(pdgrm_tot,f_pdgrm.HHZ,3,'pdgrm','low-freq')

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
plot(x,y2,'b--')
plot(x,sig_data.HHZ,'r')
plot(x,sig_data.HHN,'k')
plot(x,sig_data.HHE,'g')
title(sprintf('d - pkTot = %2.2f',d-pk_tot))
hold off
%%
fprintf('d - pkTot = %2.2f\n',d-pk_tot)
disp('using pk_tot')

%% Use the filter above on the data
fil = [pk_tot-2 pk_tot+2];
w_clean_filt = waveform_clean(w_clean, filterobject('b', fil, 2));
w_clean_nofilt = w_clean;
%This shows the ellipticity, azimuth, and incidence of the ray at that
%filter

%% Subset data to around the P- and S-wave picks

% try
%     endtime = oridStruct.(sprintf('eq_%d',eq)).time_phase(intersect(find(strcmp(oridStruct.(sprintf('eq_%d',eq)).phase,'S')),...
%         find(strcmp(oridStruct.(sprintf('eq_%d',eq)).sta,stasz))))+datenum(0,0,0,0,0,1);
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
clear count3
%% rotate waveforms
wndo = 20;%30;%20;%14;%31;
TC=threecomp(W_all,az);
%TCR = TC.rotate();
TCW = TC.waveform();
[TCP,vals,vecs] = TC.particlemotion(0.01,wndo/100);
[lambda, I]=sort(vals,1,'descend');
X = vecs(:,I(1),:); %get greatest eigenvector



TC_unfilt = threecomp(W_all_unfilt,az);
%TCR_unfilt = TC_unfilt.rotate();
TCW_unfilt = TC_unfilt.waveform();
[TCP_unfilt, vals_unfilt, vecs_unfilt] = TC_unfilt.particlemotion(0.01,wndo/100);


for count3 = 1:numel(ch)
    Haystack_data.(ch{count3}) = get(W_all(count3),'data');
    Haystack_data_unfilt.(ch{count3}) = get(W_all_unfilt(count3),'data');
    get(W_all(count3),'channel')
    needle.(ch{count3}) = Haystack_data.(ch{count3})(P_ind:P_ind+num_samps-1);
    lngX = length(Haystack_data.(ch{count3}));
    lngY = length(needle.(ch{count3}));
    assert(lngX >= lngY);
    lags = 0:(lngX-lngY);
    for valz = lags
        c.(ch{count3})(valz+1) = xcorr(Haystack_data.(ch{count3})(valz+1:valz+lngY) -...
            mean(Haystack_data.(ch{count3})(valz+1:valz+lngY)), needle.(ch{count3}) - mean(needle.(ch{count3})),0,'coeff');
    end
end

%% Check and plot directionality and polarity of waveforms
dtac = [Haystack_data.(ch{1})'; Haystack_data.(ch{3})'; Haystack_data.(ch{2})'];
dtac_unfilt = [Haystack_data_unfilt.(ch{1})'; Haystack_data_unfilt.(ch{3})'; Haystack_data_unfilt.(ch{2})'];

[azim, incd, ellip] = polar_coherency(dtac, wndo);

[azim_u, incd_u, ellip_u] = polar_coherency(dtac_unfilt,wndo);
for cz = 1:numel(azim)
    if azim(cz)<0
        azim(cz) = -azim(cz)+180;
    end
    if azim_u(cz)<0
        azim_u(cz) = -azim_u(cz)+180;
    end
%     if incd(cz)>90
%         incd(cz)=90-(incd(cz)-90);
%     end
end

%% make phase plot data range from 240-500, then plot 
%[azim_cov, incd_cov, ellip_cov] = polar_covariance(dtac, wndo);

min_num = 240;
if az < min_num
    az_adj = az+360;
else
    az_adj = az;
end

azim_adj = zeros(size(azim));
for ct = 1:numel(azim)
    if azim(ct) < min_num
        azim_adj(ct) = azim(ct)+360;
    else
        azim_adj(ct) = azim(ct);
    end
end


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
tcrp_az_data_adj = zeros(size(tcrp_az_data));
for ct = 1:numel(tcrp_az_data)
    if tcrp_az_data(ct) < min_num
        tcrp_az_data_adj(ct) = tcrp_az_data(ct)+360;
    else
        tcrp_az_data_adj(ct) = tcrp_az_data(ct);
    end
end

% try
%With phase arrivals in red
plot_threecomp_directionality(endtime,val_try,Haystack_data_unfilt,P_ind,...
    S_ind,tcrp_az,tcrp_inc, tcrp_rec, tcrp_plan, tcrp_en,az_adj,inc,tcrp_az_data_adj,ph_indz)
plot_threecomp_directionality(endtime,val_try,Haystack_data,P_ind,S_ind,...
    tcrp_az,tcrp_inc, tcrp_rec, tcrp_plan, tcrp_en,az_adj,inc,tcrp_az_data_adj,ph_indz)

%%
rec_vals = get(tcrp_rec,'data');
plan_vals = get(tcrp_plan,'data');
[valrec,indrec] = max(rec_vals(val_try:val_try+num_samps));
if plan_vals(indrec+val_try-1)<0.5
    disp('Low degree of planarity!!')
end
tcrp_inc_data = get(tcrp_inc,'data');

az_calc = tcrp_az_data_adj(indrec+val_try-1)
inc_calc = tcrp_inc_data(indrec+val_try-1)
clear phaseStruct
load('phaseStruct.mat')
phaseStruct.(sprintf('eq_%d',eq)).(stasz).az_W_calc(indie) = az_calc;
phaseStruct.(sprintf('eq_%d',eq)).(stasz).inc_W_calc(indie) = inc_calc;
phaseStruct.(sprintf('eq_%d',eq)).(stasz).ind_calc(indie) = indrec+val_try-1;
phaseStruct.(sprintf('eq_%d',eq)).(stasz).rec_W_calc(indie) = valrec;
phaseStruct.(sprintf('eq_%d',eq)).(stasz).plan_W_calc(indie) = plan_vals(indrec+val_try-1);
save('phaseStruct.mat','phaseStruct')


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
% w_prep = prep_for_rs(w_clean,sprintf('eq_%d',eq),siteStruct,oridStruct);
% 
% w_prep_E = prep_for_rs(w_clean_HHE,sprintf('eq_%d',eq),siteStruct,oridStruct);
% w_prep_N = prep_for_rs(w_clean_HHN,sprintf('eq_%d',eq),siteStruct,oridStruct);
% w_prep_Z = prep_for_rs(w_clean_HHZ,sprintf('eq_%d',eq),siteStruct,oridStruct);
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
toc
%end