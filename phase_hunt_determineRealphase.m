% function high_corrs = phase_hunt_changingNeedle_mins(allorids,oridStruct,P,fil,cutoff_val)
% find(strcmp(fieldnames(oridStruct),'eq_whatever')) %to find eq
%orid2166 is origin 203
%orid2034 is origin 128
%orid2277 is origin 277
tic
addpath('/Users/alexandrafarrell/Desktop/akfarrell/')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Polarizemic-master/functions/')

%clear;
clear all; close all; clc;
%fil=[2 25];
%fil=[10 20];
% [evidStruct, allevids] = get_eq_info();
[evidStruct_error, allevids_error] = get_eq_info('/Users/alexandrafarrell/Desktop/akfarrell/heather/dbmerged','true'); %relocated, with errors (excludes all m<0.5

[evidStruct_unchecked_error, allevids_unchecked_error] = get_eq_info('/Users/alexandrafarrell/Desktop/dbplutons_alex','true'); %not relocated, with errors
load('siteStruct.mat');
%stylie = 'min'; %%%%% CHANGE!!!---------0----min,max,abz
stylie = {'max','min'};
% load('phaseStruct.mat')
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
eq = 2501;
stasz = 'PLHS';
Pvel = 4.2; %km/s
Svel = 2.35; %km/s
%val_try = 297; % If commented out, find val_try (is P_ind or S_ind)
val_try = 482 %datenum(2011,04,29,19,15,59.686)

%%
if intersect(allevids_error,eq)
    count=find(allevids_error == eq);
    evidStruct = evidStruct_error;
    allevids = allevids_error;
elseif intersect(allevids_unchecked_error,eq)
    count=find(allevids_unchecked_error == eq);
    evidStruct = evidStruct_unchecked_error;
    allevids = allevids_unchecked_error;
end
    
directory = '/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase/wf_objs';
filename = sprintf('wf_%d.mat',allevids(count));
filename_wPath = fullfile(directory,filename);
%if exist(filename_wPath,'file')
%    load(filename_wPath)
%else
%create and clean waveform object
fil = [0.8 26];
[w_raw,EvS,stations_inEq] = get_wf(allevids(count),evidStruct);
w_clean = waveform_clean(w_raw, filterobject('b', fil, 2));
% w_clean = waveform_clean(w_raw);
%save(filename_wPath,'w_clean', 'OrS', 'stations_inEq');
%end
stationz = get(w_clean,'station');
P_pad = 12;
S_pad = 25;
%%
count2 = min(find(strcmp(stasz,stationz)));%:3:numel(w_clean) %43-45 is PLMN 43:3:43
%HHE = count, HHN = count+1, HHZ = count+2

fname = sprintf('corr_%d_%s_%s.mat',allevids(count),stationz{count2},stylie{1});
directory = '/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase/wf_objs';
directory2 = sprintf('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase/corrs/%d',allevids(count));
filename_wPath2 = fullfile(directory2,fname);
% load(filename_wPath2)

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

if val_try >2000
    val_try = find(dnum>=val_try, 1,'first'); %find index of value to try, if in dnum
end



for count3 = 1:numel(ch)
    Haystack_data.(ch{count3}) = get(w_clean(count2+count3-1),'data');
    needle.(ch{count3}) = Haystack_data.(ch{count3})(P_ind:P_ind+num_samps-1);
    %fix indexing issues
    lngX = length(Haystack_data.(ch{count3}));
    lngY = length(needle.(ch{count3}));
    assert(lngX >= lngY);
    lags = 0:(lngX-lngY);
    for valz = lags
        c.(ch{count3})(valz+1) = xcorr(Haystack_data.(ch{count3})(valz+1:valz+lngY) -...
            mean(Haystack_data.(ch{count3})(valz+1:valz+lngY)), needle.(ch{count3}) - mean(needle.(ch{count3})),0,'coeff');
    end
    
    c_backup.(ch{count3}) = c.(ch{count3}); %USE THIS VALUE TO ANNOTATE THE 3-COMPONENT PLOTS!!!!!!!!!!!-------------
    c.(ch{count3}) = abs(c.(ch{count3}));
    [m,i.(ch{count3})]=max(c.(ch{count3})(P_ind+P_pad:S_ind-S_pad)); %Pad P_ind and S_ind in range for max of abs val

    i.(ch{count3}) = i.(ch{count3})+P_ind+P_pad;
    %-------------------------UNCOMMENT PLOT IF WANT TO PLOT - 2 lines
    %plot_xcorrs(lags, ch{count3}, Haystack_data.(ch{count3}), needle.(ch{count3}), c.(ch{count3}), stationz{count2}, ...
    %cutoff_val, i.(ch{count3}), allevids(count), 1,P_ind+P_pad,S_ind-S_pad) %Make sure ind padding is same as line 64!!!
    if m>= cutoff_val %%%only half setup for multiple returns that are greater than cutoff_val - need to finish if looking at this
        values = find(c.(ch{count3})>=cutoff_val);
        %---------------------UNCOMMENT PLOT IF WANT TO PLOT - 2 lines
        %plot_xcorrs(lags, ch{count3}, Haystack_data.(ch{count3}), needle.(ch{count3}), c.(ch{count3}), stationz{count2}, ...
        %cutoff_val, i.(ch{count3}), allevids(count), 2, P_ind,S_ind)
        stations_w_highcorrs.(sprintf('eq_%d',allevids(count))).(stationz{count2}).(ch{count3}).val = m;
        stations_w_highcorrs.(sprintf('eq_%d',allevids(count))).(stationz{count2}).(ch{count3}).index = i.(ch{count3});
        %                 for indie = 1:numel(c.(ch{count3})>=cutoff_val)
        %                     stations_w_highcorrs.(sprintf('eq_%d',allevids(count))).(stationz{count2}).(ch{count3}).val(indie) = m;
        %                     stations_w_highcorrs.(sprintf('eq_%d',allevids(count))).(stationz{count2}).(ch{count3}).index(indie) = i.(ch{count3});
        %                 end
    end
    clear m
end




% c_orig=c;
c_total = c.(ch{1})+c.(ch{2})+c.(ch{3});
c_abs_total = abs(c.(ch{1}))+abs(c.(ch{2}))+abs(c.(ch{3}));
%P_pad = 15
dnum = zeros(1,numel(get(w_clean(count2),'data')));
dnum(1) = datenum(get(w_clean(count2),'start'));
freq = get(w_clean(count2), 'freq');
for l = 2:numel(get(w_clean(count2),'data'))
    dnum(l) = datenum((l/freq)/SECSPERDAY+dnum(1));
end
numz = find(strcmp(siteStruct.sta,stasz));
az = azimuth(siteStruct.lat(numz), siteStruct.lon(numz),evidStruct.(sprintf('eq_%d',eq)).lat(1),evidStruct.(sprintf('eq_%d',eq)).lon(1));
starttime = evidStruct.(sprintf('eq_%d',eq)).time_phase(intersect(find(strcmp(evidStruct.(sprintf('eq_%d',eq)).phase,'P')),...
    find(strcmp(evidStruct.(sprintf('eq_%d',eq)).sta,stasz))))-datenum(0,0,0,0,0,1);
try
    endtime = evidStruct.(sprintf('eq_%d',eq)).time_phase(intersect(find(strcmp(evidStruct.(sprintf('eq_%d',eq)).phase,'S')),...
        find(strcmp(evidStruct.(sprintf('eq_%d',eq)).sta,stasz))))+datenum(0,0,0,0,0,1);
    w_clean_e = extract(w_clean,'TIME',starttime,endtime);
catch
    endtime = dnum(end-num_samps);
    w_clean_e = extract(w_clean,'TIME',starttime,endtime);
end
for count3 = 1:numel(ch)
    
    Haystack_data.(ch{count3}) = get(w_clean(count2+count3-1),'data');
    W(1,3-count3+1) = w_clean_e(count2+count3-1);
    needle.(ch{count3}) = Haystack_data.(ch{count3})(P_ind:P_ind+num_samps-1);
    %%
    c_t.(ch{count3}) = c_total;
    c_a_t.(ch{count3}) = c_abs_total;
    for durr = 1:numel(stylie)
        if strcmp(stylie{durr},'max')
            c.HHE(val_try) =1;
            c.HHN(val_try)=1;
            c.HHZ(val_try)=1;
        elseif strcmp(stylie{durr},'min')
            c.HHE(val_try) =-1;
            c.HHN(val_try)=-1;
            c.HHZ(val_try)=-1;
        end
        if strcmp(stylie{durr}, 'min')
            [m,i2.(ch{count3})]=min(c.(ch{count3})(P_ind+P_pad:S_ind-S_pad)); %Pad P_ind and S_ind in range for min
            [m,itmin.(ch{count3})]=min(c_total(P_ind+P_pad:S_ind-S_pad));
            itmin.(ch{count3}) = itmin.(ch{count3})+P_ind+P_pad;
        elseif strcmp(stylie{durr}, 'max')
            [m,i2.(ch{count3})]=max(c.(ch{count3})(P_ind+P_pad:S_ind-S_pad)); %Pad P_ind and S_ind in range for max
            [m,itmax.(ch{count3})]=max(c_total(P_ind+P_pad:S_ind-S_pad));
            [m,iatmax.(ch{count3})]=max(c_abs_total(P_ind+P_pad:S_ind-S_pad));
            itmax.(ch{count3}) = itmax.(ch{count3})+P_ind+P_pad;
            iatmax.(ch{count3}) = iatmax.(ch{count3})+P_ind+P_pad;
        elseif strcmp(stylie{durr}, 'abz')
            c_backup.(ch{count3}) = c.(ch{count3}); %USE THIS VALUE TO ANNOTATE THE 3-COMPONENT PLOTS!!!!!!!!!!!-------------
            c.(ch{count3}) = abs(c.(ch{count3}));
            [m,i2.(ch{count3})]=max(c.(ch{count3})(P_ind+P_pad:S_ind-S_pad)); %Pad P_ind and S_ind in range for max of abs val
        end
        i2.(ch{count3}) = i2.(ch{count3})+P_ind+P_pad;
        
        clear m
        
    end
end
% 
% %% Check and plot directionality and polarity of waveforms
% dtac = [Haystack_data.HHZ'; Haystack_data.HHE'; Haystack_data.HHN'];
% wndo = 20;%14;%31;
% [azim, incd, ellip] = polar_coherency(dtac, wndo);
% [azim_cov, incd_cov, ellip_cov] = polar_covariance(dtac, wndo);
% 
% 
% phase_plots(Haystack_data,P_ind, S_ind, azim,incd,ellip,eq,stasz,val_try,xmax, az, inc)
%phase_plots(Haystack_data,P_ind, S_ind, azim_cov,incd_cov,ellip_cov,eq,stasz)

%%

plot_xcorrs(lags, ch, Haystack_data, needle, c, stationz{count2}, cutoff_val, i2, allevids(count), 4,P_ind+P_pad,S_ind-S_pad)
plot_xcorrs(lags, ch, Haystack_data, needle, c_t, stationz{count2}, cutoff_val*3, itmin, allevids(count), 4,P_ind+P_pad,S_ind-S_pad)
plot_xcorrs(lags, ch, Haystack_data, needle, c_t, stationz{count2}, cutoff_val*3, itmax, allevids(count), 4,P_ind+P_pad,S_ind-S_pad)
% if exist('iatmax','var')
%     plot_xcorrs(lags, ch, Haystack_data, needle, c_a_t, stationz{count2}, cutoff_val*3, iatmax, allorids(count), 4,P_ind+P_pad,S_ind-S_pad)
% end

h = figure();
set(h, 'Position', [500 5 2000 1300])
for countz = 1:numel(ch)
    subplot(1,6,countz*2-1)
    plot(Haystack_data.(ch{countz}),1:numel(Haystack_data.(ch{countz})))
    if S_ind + 300 > numel(Haystack_data.(ch{countz}))
        ylim([P_ind-300, S_ind])
    else
        ylim([P_ind-300, S_ind+300])
    end
    xlim([min(Haystack_data.(ch{countz})), max(Haystack_data.(ch{countz}))])
    line([min(Haystack_data.(ch{countz})), max(Haystack_data.(ch{countz}))], [P_ind, P_ind],'Color','m','LineStyle',':','LineWidth',2)
    line([min(Haystack_data.(ch{countz})), max(Haystack_data.(ch{countz}))],[S_ind, S_ind],'Color','m','LineStyle',':','LineWidth',2)
    title(sprintf('%s',ch{countz}))
    %ylim([1 numel(get(w_clean(1),'data'))])
    hold on
    subplot(1,6,countz*2)
    spectrogram(Haystack_data.(ch{countz}),20,19,128,100) %make sure this is same as below
    s.(ch{countz})=spectrogram(Haystack_data.(ch{countz}),20,19,128,100); %make sure this is same as above
    if S_ind + 300 > numel(Haystack_data.(ch{countz}))
        ylim([(P_ind-300)/100, S_ind/100])
    else
        ylim([(P_ind-300)/100, (S_ind+300)/100])
    end
    hold on
    line([0, 50], [P_ind/100, P_ind/100],'Color','m','LineStyle',':','LineWidth',2)
    line([0, 50], [S_ind/100, S_ind/100],'Color','m','LineStyle',':','LineWidth',2)
    %% lines of stuff
    line([0, 50], [iatmax.(ch{countz})/100, iatmax.(ch{countz})/100],'Color','y','LineStyle',':','LineWidth',2)
    line([0, 50], [itmax.(ch{countz})/100, itmax.(ch{countz})/100],'Color','r','LineStyle',':','LineWidth',2)
    line([0, 50], [itmin.(ch{countz})/100, itmin.(ch{countz})/100],'Color','g','LineStyle',':','LineWidth',2)
    line([0, 50], [val_try/100, val_try/100],'Color','w','LineStyle',':','LineWidth',2)
    xlim([0,50])
    col = caxis;
    caxis([-15 col(2)])
    clear col
end
hold off

%% Trying things
g = figure();
set(g, 'Position', [1000 1000 1400 1200])
for countzs = 1:numel(ch)
    subplot(6,1,countzs*2-1)
    plot(abs(s.(ch{countzs})(1,:)),'k'); grid on;
    if S_ind + 300 > numel(Haystack_data.(ch{countzs}))
        xlim([P_ind-300, S_ind])
    else
        xlim([P_ind-300, S_ind+300])
    end
    ylim([0, max(abs(s.(ch{countzs})(1,:)))])
    line([iatmax.(ch{countzs}), iatmax.(ch{countzs})],[0, max(abs(s.(ch{countzs})(1,:)))],'Color','y','LineStyle',':','LineWidth',2)
    line([itmax.(ch{countzs}), itmax.(ch{countzs})],[0, max(abs(s.(ch{countzs})(1,:)))],'Color','r','LineStyle',':','LineWidth',2)
    line([itmin.(ch{countzs}), itmin.(ch{countzs})],[0, max(abs(s.(ch{countzs})(1,:)))],'Color','g','LineStyle',':','LineWidth',2)
    line([val_try, val_try],[0, max(abs(s.(ch{countzs})(1,:)))],'Color','b','LineStyle',':','LineWidth',2)
    hold on
    
    subplot(6,1,countzs*2)
    plot(1:numel(Haystack_data.(ch{countzs})),Haystack_data.(ch{countzs}))
    if S_ind + 300 > numel(Haystack_data.(ch{countzs}))
        xlim([P_ind-300, S_ind])
    else
        xlim([P_ind-300, S_ind+300])
    end
    hold on
    plot(val_try:val_try+30,Haystack_data.(ch{countzs})(val_try:val_try+30),'r')
    ylim([min(Haystack_data.(ch{countzs})), max(Haystack_data.(ch{countzs}))])
    line([P_ind, P_ind], [min(Haystack_data.(ch{countzs})), max(Haystack_data.(ch{countzs}))], 'Color','m','LineStyle',':','LineWidth',2)
    line([S_ind, S_ind], [min(Haystack_data.(ch{countzs})), max(Haystack_data.(ch{countzs}))],'Color','m','LineStyle',':','LineWidth',2)
    title(sprintf('%s',ch{countzs}))
end
hold off

% %% 
% f = figure();
% set(f, 'Position', [1000 1000 1400 1200])
% for countzs2 = 1:numel(ch)
%     subplot(6,1,countzs2*2-1)
%     b = abs(s.HHZ(1,:))+abs(s.HHN(1,:))+abs(s.HHE(1,:));
%     plot(abs(s.HHZ(1,:))+abs(s.HHN(1,:))+abs(s.HHE(1,:)),'k'); grid on;
%     if S_ind + 300 > numel(Haystack_data.(ch{countzs2}))
%         xlim([P_ind-300, S_ind])
%     else
%         xlim([P_ind-300, S_ind+300])
%     end
%     ylim([0, max(abs(s.HHZ(1,:))+abs(s.HHN(1,:))+abs(s.HHE(1,:)))])
%     line([iatmax.(ch{countzs2}), iatmax.(ch{countzs2})],[0, max(abs(s.HHZ(1,:))+abs(s.HHN(1,:))+abs(s.HHE(1,:)))],'Color','y','LineStyle',':','LineWidth',2)
%     line([itmax.(ch{countzs2}), itmax.(ch{countzs2})],[0, max(abs(s.HHZ(1,:))+abs(s.HHN(1,:))+abs(s.HHE(1,:)))],'Color','r','LineStyle',':','LineWidth',2)
%     line([itmin.(ch{countzs2}), itmin.(ch{countzs2})],[0, max(abs(s.HHZ(1,:))+abs(s.HHN(1,:))+abs(s.HHE(1,:)))],'Color','g','LineStyle',':','LineWidth',2)
%     line([val_try, val_try],[0, max(abs(s.HHZ(1,:))+abs(s.HHN(1,:))+abs(s.HHE(1,:)))],'Color','b','LineStyle',':','LineWidth',2)
%     hold on
%     
%     subplot(6,1,countzs2*2)
%     plot(1:numel(Haystack_data.(ch{countzs2})),Haystack_data.(ch{countzs2}))
%     if S_ind + 300 > numel(Haystack_data.(ch{countzs2}))
%         xlim([P_ind-300, S_ind])
%     else
%         xlim([P_ind-300, S_ind+300])
%     end
%     ylim([min(Haystack_data.(ch{countzs2})), max(Haystack_data.(ch{countzs2}))])
%     line([P_ind, P_ind], [min(Haystack_data.(ch{countzs2})), max(Haystack_data.(ch{countzs2}))], 'Color','m','LineStyle',':','LineWidth',2)
%     line([S_ind, S_ind], [min(Haystack_data.(ch{countzs2})), max(Haystack_data.(ch{countzs2}))],'Color','m','LineStyle',':','LineWidth',2)
%     title(sprintf('%s',ch{countzs2}))
% end
% hold off

%% Check and plot the above, but with GISMO
% TC=threecomp(W,az);
% TCR = TC.rotate();
% TCRP = TCR.particlemotion();
% TCRP.plot2()

%%
toc

%end