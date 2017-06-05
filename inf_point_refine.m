%% RUN APMB_schematic BEFORE THIS!!!
tic
close all
clear all
load('inf_point_2166and2034.mat')
[oridStruct, allorids] = get_eq_info();
%orid2166 is origin 203
%orid2034 is origin 128
for count=128%1:numel(allorids)
    siteSub = get_station_info(oridStruct, allorids(count));
    clear stationz
    OrS = sprintf('eq_%d',allorids(count));
    stationz = fieldnames(inf_point.(OrS));
    
    %% Stuff for Orid2034
    if allorids(count) == 2034
    site_subset_phase = {'PL03','PLMK'};
    site_subset_maybe = {'PLLL','PLQU','PLSS'};
    end
    
    %% Stuff for Orid2166
    if allorids(count) == 2166
    site_subset_phase = {'PL03','PLLL','PLMD','PLSP'};
    site_subset_3phase = {'PLDK','PLJR','PLMK','PLMN'};
    site_subset_maybe = {'PLSM'};
    %time_vals = {'', '0.51s or 1.23s', '0.47s or 0.79s', '0.45s or 0.72s'};
    end

    %% Type of phase
    for count2 = 1:numel(stationz)
        siteInd = find(strcmp(siteSub.sta, stationz{count2}));
        [latout,lonout] = reckon(oridStruct.(OrS).lat(1),oridStruct.(OrS).lon(1),inf_point.(OrS).(stationz{count2}).infl(1)/111.12,inf_point.(OrS).(stationz{count2}).infl(2));
        inf_point.(OrS).(stationz{count2}).infl(3) = latout;
        inf_point.(OrS).(stationz{count2}).infl(4) = lonout;
        if find(strcmp(site_subset_phase,stationz{count2}))
            inf_point.(OrS).(stationz{count2}).phase = 'y';
        elseif exist('site_subset_3phase','var') && numel(unique(strcmp(site_subset_3phase,stationz{count2}))) == 2
            inf_point.(OrS).(stationz{count2}).phase = 't';
        elseif exist('site_subset_maybe','var') && numel(unique(strcmp(site_subset_maybe,stationz{count2}))) == 2
            inf_point.(OrS).(stationz{count2}).phase = 'm';
        else
            inf_point.(OrS).(stationz{count2}).phase = 'n';
        end
        clear latout; clear lonout
    end
    clear siteSub
end
toc