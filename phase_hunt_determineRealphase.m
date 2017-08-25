% function high_corrs = phase_hunt_changingNeedle_mins(allorids,oridStruct,P,fil,cutoff_val)
% find(strcmp(fieldnames(oridStruct),'eq_whatever')) %to find eq
%orid2166 is origin 203
%orid2034 is origin 128
%orid2277 is origin 277
tic
addpath('/raid/home/a/akfarrell/')
addpath('/raid/home/a/akfarrell/Uturuncu')
clear all;
clc;
fil=[2 25];
%fil=[10 20];
[oridStruct, allorids] = get_eq_info(); 
load('siteStruct.mat');
%stylie = 'min'; %%%%% CHANGE!!!---------0----min,max,abz
stylie = {'max','min'};

ch = {'HHE','HHN','HHZ'};
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
%close all
% directory = '/home/a/akfarrell/Uturuncu/Phase/wf_objs';
% filename = sprintf('wf_%d.mat',allorids(count));
% filename_wPath = fullfile(directory,filename);
% if exist(filename_wPath,'file')
%     load(filename_wPath)
% else
%     create and clean waveform object
    [w_raw,OrS,stations_inEq] = get_wf(allorids(count),oridStruct);
    w_clean = waveform_clean(w_raw, filterobject('b', fil, 2));
%     save(filename_wPath,'w_clean', 'OrS', 'stations_inEq');
% end
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
c_total = c.(ch{1})+c.(ch{2})+c.(ch{3});
c_abs_total = abs(c.(ch{1}))+abs(c.(ch{2}))+abs(c.(ch{3}));
%P_pad = 15
dnum = zeros(1,numel(get(w_clean(count2),'data')));
dnum(1) = datenum(get(w_clean(count2),'start'));
freq = get(w_clean(count2), 'freq');
for l = 2:numel(get(w_clean(count2),'data'))
    dnum(l) = datenum((l/freq)/SECSPERDAY+dnum(1));
end

for count3 = 1:numel(ch)
    Haystack_data.(ch{count3}) = get(w_clean(count2+count3-1),'data');
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

%% Check and plot directionality and polarity of waveforms
dtac = [Haystack_data.HHZ'; Haystack_data.HHE'; Haystack_data.HHN'];
wndo = 20;%31;
[azim, incd, ellip] = polar_coherency(dtac, wndo);
[azim_cov, incd_cov, ellip_cov] = polar_covariance(dtac, wndo);
%phase_plots(Haystack_data,P_ind, S_ind, azim,incd,ellip,eq,stasz)
phase_plots(Haystack_data,P_ind, S_ind, azim_cov,incd_cov,ellip_cov,eq,stasz)
numz = find(strcmp(siteStruct.sta,stasz));
az = azimuth(siteStruct.lat(numz), siteStruct.lon(numz),oridStruct.(sprintf('eq_%d',eq)).lat(1),oridStruct.(sprintf('eq_%d',eq)).lon(1));

%%

% plot_xcorrs(lags, ch, Haystack_data, needle, c, stationz{count2}, cutoff_val, i2, allorids(count), 4,P_ind+P_pad,S_ind-S_pad)
% plot_xcorrs(lags, ch, Haystack_data, needle, c_t, stationz{count2}, cutoff_val*3, itmin, allorids(count), 4,P_ind+P_pad,S_ind-S_pad)
% plot_xcorrs(lags, ch, Haystack_data, needle, c_t, stationz{count2}, cutoff_val*3, itmax, allorids(count), 4,P_ind+P_pad,S_ind-S_pad)
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

%%
toc

%end