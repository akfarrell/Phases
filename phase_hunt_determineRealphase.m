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
[oridStruct, allorids] = get_eq_info();

%stylie = 'min'; %%%%% CHANGE!!!---------0----min,max,abz
stylie = {'max','min'};

ch = {'HHE','HHN','HHZ'};
cutoff_val = 0.8;
S_pad = 25;
SECSPERDAY = 60 * 60 * 24;
num_samps = 31;
close all

eq = 1982;
stasz = 'PLSE';
val_try = 672;

%%
count=find(allorids == eq);
close all
directory = '/home/a/akfarrell/Uturuncu/Phase/wf_objs';
filename = sprintf('wf_%d.mat',allorids(count));
filename_wPath = fullfile(directory,filename);
if exist(filename_wPath,'file')
    load(filename_wPath)
else
    create and clean waveform object
    [w_raw,OrS,stations_inEq] = get_wf(allorids(count),oridStruct);
    w_clean = waveform_clean(w_raw, filterobject('b', fil, 2));
    save(filename_wPath,'w_clean', 'OrS', 'stations_inEq');
end
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

plot_xcorrs(lags, ch, Haystack_data, needle, c, stationz{count2}, cutoff_val, i2, allorids(count), 4,P_ind+P_pad,S_ind-S_pad)
plot_xcorrs(lags, ch, Haystack_data, needle, c_t, stationz{count2}, cutoff_val*3, itmin, allorids(count), 4,P_ind+P_pad,S_ind-S_pad)
plot_xcorrs(lags, ch, Haystack_data, needle, c_t, stationz{count2}, cutoff_val*3, itmax, allorids(count), 4,P_ind+P_pad,S_ind-S_pad)
if exist('iatmax','var')
    plot_xcorrs(lags, ch, Haystack_data, needle, c_a_t, stationz{count2}, cutoff_val*3, iatmax, allorids(count), 4,P_ind+P_pad,S_ind-S_pad)
end

%%
% inds = find(c_total) >= overall_cutoff;
% inds_backup = inds;
% %remove inds that correspond to P-wave pick
% P_array = P_ind-P_pad:P_ind+P_pad;
% for indexy = 1:numel(P_array)
%     if intersect(inds,P_array(indexy))
%         inds = inds(inds~=intersect(inds,P_array(indexy)));
%     end
% end
% count4 = 1;
% %remove inds that are too close to each other
% %         while count4 <= numel(inds)
% %             count4;
% %             if numel(intersect(inds,inds(count4)-P_pad:inds(count4)+P_pad)) > 1
% %                 int_vals = intersect(inds,inds(count4)-P_pad:inds(count4)+P_pad);
% %                 [dumb,throw_val] = min(c_total(intersect(inds,int_vals)));
% %                 inds = inds(inds~=int_vals(throw_val));
% %                 count4 = count4+1;
% %              end
% %             count4 = count4+1;
% %         end
% 
% for counter = 1:numel(inds)
%     data_sig(counter*numel(needle.HHZ)-numel(needle.HHZ)+1:counter*numel(needle.HHZ)) = inds(counter):inds(counter)+numel(needle.HHZ)-1;
% end
%%
toc

%end