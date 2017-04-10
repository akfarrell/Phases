%Load and sort origin information - depths >=7.5 km from sea level
%460 events in dbplutons_alex <=7.5 km
%clc
clear all
close all

addpath('/raid/home/a/akfarrell/')
addpath('/raid/home/a/akfarrell/Uturuncu')
addpath('/raid/home/a/akfarrell/Uturuncu/Phase/rgb')
addpath('../scalebar_v3/scalebar')
addpath('../pathdist_v4')
addpath('../Inpaint_nans/Inpaint_nans')
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
% p=figure; hold on;
% set(p, 'Position', [1000 1000 1000 1000])
% latlim = [-22.75 -21.75]; %[southern_limit northern_limit] 
% lonlim = [-67.75 -66.5]; %[western_limit eastern_limit]
% worldmap(latlim, lonlim);
% borders = shaperead('../BOL_adm0.shp', 'UseGeoCoords', true);
% colormap(flipud(colormap))
% arg_borders = shaperead('../ARG_adm0.shp', 'UseGeoCoords', true);
% geoshow(borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
% geoshow(arg_borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
for count=1:numel(allorids)
    thisorid = allorids(count);
    oridS=sprintf('eq_%i',thisorid);
    indexes = find(PhaseStruct.orid == thisorid);
    for ncount=1:numel(names)
        oridStruct.(oridS).(names{ncount}) = PhaseStruct.(names{ncount})(indexes);
    end
%     scatterm(oridStruct.(oridS).lat,oridStruct.(oridS).lon,'k')
%     hold on
end
% scatterm(oridStruct.eq_2166.lat,oridStruct.eq_2166.lon,'filled','r')
% scatterm(-22.27, -67.18, 100,'^','filled','r')
%--------confidence for smajax,sminax,strike,sdepth,and stime is 0.683 ----
clear count
clear PhaseStruct
clear names

%Create waveform objects
%Specify orid
Or = 2166
%2166,2034
%%
[w_raw,OrS,stations_inEq] = get_wf(Or,oridStruct);

%%
%clean waveforms and sort by distance from earthquake
%close all;
fil=[2 25];
w_clean = waveform_clean(w_raw, filterobject('b', fil, 2));
%%
sta=unique(get(w_clean,'station'));
siteStruct = loadSiteTable('/raid/home/a/akfarrell/Uturuncu/dbplutons_alex','lat>=-23');
% scatterm(siteStruct.lat,siteStruct.lon,100,'*','m')
% textm(siteStruct.lat-0.04,siteStruct.lon-0.025,siteStruct.sta,'m')
% hold off
for count=1:numel(stations_inEq)
    siteSub.sta(count) = siteStruct.sta(strcmp(stations_inEq(count), siteStruct.sta));
    siteSub.lat(count) = siteStruct.lat(strcmp(stations_inEq(count), siteStruct.sta));
    siteSub.lon(count) = siteStruct.lon(strcmp(stations_inEq(count), siteStruct.sta));
    siteSub.elev(count) = siteStruct.elev(strcmp(stations_inEq(count), siteStruct.sta));
    siteSub.dist(count) = distance(oridStruct.(OrS).lat(1), oridStruct.(OrS).lon(1), ...
    siteStruct.lat(strcmp(stations_inEq(count), siteStruct.sta)), siteStruct.lon(strcmp(stations_inEq(count), siteStruct.sta)));
end
clear count
siteSub.distKM=siteSub.dist*111.12;

if isequal(sta,siteSub.sta)
    n=3;
    stas_site=reshape(repmat(siteSub.sta(:).',n,1),1,[]);
    dist_site=reshape(repmat(siteSub.dist(:).',n,1),1,[]);
else
    clear w_clean
end
[Y,order] = sort(dist_site);
w_clean_sort = w_clean(order);
clear w_clean w_raw 

len = numel(w_clean_sort);
max_vals=[];
min_vals=[];
for i=1:len
    max_vals(i) = max(w_clean_sort(i));
    min_vals(i) = min(w_clean_sort(i));
end
%stats n things
[stat_struct,most_phases_eq,num_phase] = find_info_shit(oridStruct);
%%
timez = automated_arrivals(siteSub,oridStruct.(OrS).depth(1)); %%---------currently set for smaller utu vel model!!!!!----------------------

%Find time in s between OT and picked phase arrival
park=datevec(oridStruct.(OrS).time_phase-oridStruct.(OrS).time_origin(1));
oridStruct.(OrS).time_from_ot = park(:,6);


%Plot waveform objects
style={'wf','linez','linez_shifted','plotz'}; %wf is regular waveforms, linez is with phases
%station = 'all';
station = {'PLMN','PLSM','PLMD','PL03','PLSP'};
%%
phase_mulplt(w_clean_sort,0,max_vals,min_vals,fil,OrS,timez,oridStruct,style{1},station)
%%
close all
phase_mulplt(w_clean_sort,0,max_vals,min_vals,fil,OrS,timez,oridStruct,style{2},station)
%%
phase_mulplt(w_clean_sort,0,max_vals,min_vals,fil,OrS,timez,oridStruct,style{3},station)

%%
%Plot error ellipse
%title with order/total eqs for each dimension of the error ellipse
h = figure;
latlim = [min(siteStruct.lat)-0.15, max(siteStruct.lat)+0.15]; %[southern_limit northern_limit] 
lonlim = [min(siteStruct.lon)-0.15, max(siteStruct.lon)+0.15]; %[western_limit eastern_limit]
set(h, 'Position', [1000 1000 1000 1000])
worldmap(latlim, lonlim);
borders = shaperead('BOL_adm0.shp', 'UseGeoCoords', true);
arg_borders = shaperead('ARG_adm0.shp', 'UseGeoCoords', true);
%geoshow(borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
%geoshow(arg_borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
hold on
scatterm(siteSub.lat, siteSub.lon, '^', 'k')
scatterm(oridStruct.(OrS).lat(1), oridStruct.(OrS).lon(1), 'm', 'filled')
northarrow('latitude', -22.7, 'longitude', -67.6, 'scaleratio', 1/20);
scalebar('location', 'se', 'FontSize', 11)
textm(siteSub.lat+0.01, siteSub.lon, siteSub.sta)
directory = '/home/a/akfarrell/Uturuncu/Phase';
filename = sprintf('%s_location_%1.4f_%1.4f.png',OrS,fil(1),fil(2));
filename_wPath = fullfile(directory,filename);
hgexport(h, filename_wPath, hgexport('factorystyle'), 'Format', 'png');
clear file*
clear directory

%%
%------- Show particle motion ----------
%find and cut data to area of interest
close all
%st_int = 'PLSP';
%st_int = 'PL03';
st_int = 'PLMN';
stations = get(w_clean_sort, 'station');
channels = get(w_clean_sort, 'channel');
st_ind = find(strcmp(stations,st_int));
chz = channels(st_ind);
SECSPERDAY = 60 * 60 * 24;
%datenum(2011,4,29,19,16,2.27) sttime(1) PLSP
%datenum(2011,4,29,19,16,2.40) ettime(1) PLSP
time_phase = oridStruct.(OrS).time_phase(find(strcmp(oridStruct.(OrS).sta,st_int)));
ph_time.PLMN.s = datenum(2011,4,29,19,15,56.03);
ph_time.PLMN.e = datenum(2011,4,29,19,15,56.33); %12,33,20
ph_time.PL03.s = datenum(2011,4,29,19,15,59.97);
ph_time.PL03.e = datenum(2011,4,29,19,16,0.09);
ph_time.PLSP.s = datenum(2011,4,29,19,16,2.27);
ph_time.PLSP.e = datenum(2011,4,29,19,16,2.40);
if numel(time_phase) > 1
%     if strcmp(st_int, 'PLMN')
%         time_phase(2) = time_phase(2)-datenum(0,0,0,0,0,0.045)
%     end
    sttime = [ph_time.(st_int).s,time_phase(1),time_phase(2)]; %(30-43)
    ettime = [ph_time.(st_int).e,time_phase(1) + datenum(0,0,0,0,0,0.33), time_phase(2)+datenum(0,0,0,0,0,0.2)]; %time_phase(1) + datenum(0,0,0,0,0,0.1)%for motion plots
else
    sttime = [ph_time.(st_int).s,time_phase(1)]; %(30-43)
    ettime = [ph_time.(st_int).e,time_phase(1) + datenum(0,0,0,0,0,0.33)]; %time_phase(1) + datenum(0,0,0,0,0,0.1)%for motion plots
end
for i = 1:numel(sttime)
    %fwave = subtime(w_clean_sort(st_ind),sttime(i),ettime(i));
    dataz = get(w_clean_sort(st_ind(1:3)), 'data');
    start = get(w_clean_sort(st_ind(1)),'start');
    dnum(1)=datenum(start);
    freq = get(w_clean_sort(st_ind(1)),'freq');
    for l = 2:numel(dataz{1})
        dnum(l) = datenum((l/freq)/SECSPERDAY+dnum(1));
    end
    [val1,minz] = min(abs(dnum-sttime(i)));
    [val2,maxz] = min(abs(ettime(i)-dnum));
    

    dataHHE = dataz{1}(minz:maxz); %get(fwave(1), 'data');
    dataHHN = dataz{2}(minz:maxz); %get(fwave(2),'data');
    dataHHZ = dataz{3}(minz:maxz); %get(fwave(3),'data');
    %if i == 1 %old way of getting needle - CHANGE WHEN CHANGING NEEDLE
    %if i is 1, looks at phase. If i is 2, looks at P wave
    if i == 2
        P.dataHHE = dataHHE;
        P.dataHHN = dataHHN;
        P.dataHHZ = dataHHZ;
        P.minz = minz;
        P.maxz = maxz;
        P.minTime = val1;
        P.maxTime = val2;
        P.start = start;
    end

    phase_mulplt(w_clean_sort(st_ind),0,max_vals(st_ind),min_vals(st_ind),fil,OrS,timez,oridStruct,style{4},{st_int},[minz maxz])
    a= floor(numel(dataHHE)/3);
    b= floor(numel(dataHHE)/1.5);

    h = figure();
    set(h, 'Position', [1000 1000 1250 1250])
    subplot(2,2,1)
%     plot3(dataHHE(1:a), dataHHN(1:a), dataHHZ(1:a), 'k') %Change so later in time is different color
%     hold on
%     plot3(dataHHE(a:b), dataHHN(a:b), dataHHZ(a:b), 'r') %Change so later in time is different color
%     plot3(dataHHE(b:end), dataHHN(b:end), dataHHZ(b:end), 'm') %Change so later in time is different color
    plot3(dataHHE,dataHHN,dataHHZ)
    hold on
    scatter3(dataHHE,dataHHN,dataHHZ,'o')
    c = [1:numel(dataHHE)];
    d = num2cell(c);
    %e = cellstr(d);
    dx = 10; dy = 10;dz=10;
    text(dataHHE+dx,dataHHN+dy,dataHHZ+dz,d);
    grid on
    box on
    xlabel '- W, + E'
    ylabel '- S, + N'
    zlabel '- down, + up'

    subplot(2,2,2)

    plot(dataHHE, dataHHN) %Change so later in time is different color
    hold on
    scatter(dataHHE, dataHHN, 'o')
    text(dataHHE+dx,dataHHN+dy,d);
    grid on
    box on
    xlabel '- W, + E'
    ylabel '- S, + N'

    subplot(2,2,3)
    plot(dataHHN, dataHHZ) %Change so later in time is different color
    hold on
    scatter(dataHHN, dataHHZ, 'o')
    text(dataHHN+dy,dataHHZ+dz,d);
    grid on
    box on
    xlabel '- S, + N'
    ylabel '- down, + up'


    subplot(2,2,4)
    plot(dataHHE, dataHHZ) %Change so later in time is different color
    hold on
    scatter(dataHHE, dataHHZ, 'o')
    text(dataHHE+dx,dataHHZ+dz,d);
    grid on
    box on
    xlabel '- W, + E'
    ylabel '- down, + up'
    hold off
end

%%
%phase_hunt()
%%
figure()
s=spectralobject(8,6,25,[40 80]);
iceweb.spectrogram_iceweb(s,w_clean_sort(st_ind),0.5)

%%
%run perl scripts to find the inflection point, travel time, etc.
%result=perl('tt_reflected.pl', 'velmodel_PLUTONS_forPerl.txt', '10', '8', 'P', '75', '1')

%%

