%Load and sort origin information - depths >=7.5 km from sea level
%460 events in dbplutons_alex <=7.5 km
%clc
clear all
close all
addpath('/raid/home/a/akfarrell/')
addpath('/raid/home/a/akfarrell/Uturuncu')
depth_to_reflector = 7.5;
subset_exp = {sprintf('depth<=%1.1f',depth_to_reflector),'lat<=-21.8955 && lat>=-22.6698','lon<=-66.7865 && lon>=-67.6381'};

[PhaseStruct, phase, sta, chan] = loadArrivals_joinTables('/raid/home/a/akfarrell/Uturuncu/dbplutons_alex',subset_exp); %ordered by time
%429 events, 6679 rows
PhaseStruct.phase = phase;
PhaseStruct.sta = sta;
PhaseStruct.chan = chan;
PhaseStruct.time_phase = epoch2datenum(PhaseStruct.time_phase);
PhaseStruct.time_origin = epoch2datenum(PhaseStruct.time_origin);

allorids = unique(PhaseStruct.orid);
names=fields(PhaseStruct);
for count=1:numel(allorids)
    thisorid = allorids(count);
    oridS=sprintf('eq_%i',thisorid);
    indexes = find(PhaseStruct.orid == thisorid);
    for ncount=1:numel(names)
        oridStruct.(oridS).(names{ncount}) = PhaseStruct.(names{ncount})(indexes);
    end
end
%--------confidence for smajax,sminax,strike,sdepth,and stime is 0.683 ----
clear count
clear PhaseStruct
clear names
%%
%Create waveform objects
%Specify orid
Or = 2166
[w_raw,OrS,stations_inEq] = get_wf(Or,oridStruct)

%%
%clean waveforms and sort by distance from earthquake
%close all;
fil=[2 25];
w_clean = waveform_clean(w_raw, filterobject('b', fil, 2));
sta=unique(get(w_clean,'station'));
siteStruct = loadSiteTable('/raid/home/a/akfarrell/Uturuncu/dbplutons_alex','lat>=-23');
for count=1:numel(stations_inEq)
    siteSub.sta(count) = siteStruct.sta(strcmp(stations_inEq(count), siteStruct.sta));
    siteSub.lat(count) = siteStruct.lat(strcmp(stations_inEq(count), siteStruct.sta));
    siteSub.lon(count) = siteStruct.lon(strcmp(stations_inEq(count), siteStruct.sta));
    siteSub.elev(count) = siteStruct.elev(strcmp(stations_inEq(count), siteStruct.sta));
    siteSub.dist(count) = distance(oridStruct.(OrS).lat(1), oridStruct.(OrS).lon(1), ...
        siteStruct.lat(strcmp(stations_inEq(count), siteStruct.sta)), siteStruct.lon(strcmp(stations_inEq(count), siteStruct.sta)));
end
clear count

if isequal(sta,siteSub.sta)
    n=3;
    stas_site=reshape(repmat(siteSub.sta(:).',n,1),1,[]);
    dist_site=reshape(repmat(siteSub.dist(:).',n,1),1,[]);
else
    clear w_clean
end
[Y,order] = sort(dist_site);
w_clean_sort = w_clean(order);

len = numel(w_clean_sort);
max_vals=[];
min_vals=[];
for i=1:len
    max_vals(i) = max(w_clean_sort(i));
    min_vals(i) = min(w_clean_sort(i));
end


%Plot waveform objects
phase_mulplt(w_clean_sort,0,max_vals,min_vals)

%%
%run perl scripts to find the inflection point, travel time, etc.
%result=perl('tt_reflected.pl', 'velmodel_PLUTONS_forPerl.txt', '10', '8', 'P', '75', '1')