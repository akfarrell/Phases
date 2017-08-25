%% Loading in inflection point data
tic
%close all
clear inf_point %@!!!!!!!!!!-----make sure that inf_point is saved if there are changes-------!!!!!
%load('inf_point_2166and2034.mat')
load('g_oridStruct.mat')
load('good_orids.mat')
%load('oridStruct.mat')
load('inf_point.mat')
erqs = fieldnames(inf_point);
opt = 'zoom' %if I want to zoom in on a specific area

%% Plot 
p=figure; hold on;
set(p, 'Position', [1000 1000 1200 1200])
if strcmp(opt,'zoom')
    latlim = [-22.45 -22.1]; %[southern_limit northern_limit] 
    lonlim = [-67.42 -66.95]; %[western_limit eastern_limit]
else
    latlim = [-22.75 -21.85]; %[southern_limit northern_limit] 
    lonlim = [-67.75 -66.75]; %[western_limit eastern_limit]
end
worldmap(latlim, lonlim);
borders = shaperead('../BOL_adm0.shp', 'UseGeoCoords', true);
colormap(flipud(colormap))
arg_borders = shaperead('../ARG_adm0.shp', 'UseGeoCoords', true);
try
    geoshow(borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
catch
end
try
    geoshow(arg_borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
catch
end
hold on
scatterm(-22.27, -67.18, 100,'^','r')

%% Plotting reflection points
for count=41:66%numel(good_orids)
    stationz = fieldnames(inf_point.(erqs{count}));
    %scatterm(g_oridStruct.(erqs{count}).lat,g_oridStruct.(erqs{count}).lon,100,'d','filled','r')
    for count2 = 1:numel(stationz)
        latz = inf_point.(erqs{count}).(stationz{count2}).infl(3);
        lonz = inf_point.(erqs{count}).(stationz{count2}).infl(4);
        if find(strcmp(fieldnames(g_oridStruct.(erqs{count})),'refl')) && any(strcmp(fieldnames(g_oridStruct.(erqs{count}).refl),stationz{count2}))~=0
            if any(strcmp(g_oridStruct.(erqs{count}).refl.(stationz{count2}).status,'y'))~=0 && numel(g_oridStruct.(erqs{count}).refl.(stationz{count2}).time) == 1 %yes
                  scatterm(latz,lonz,30,'o','m','f') %only one, and it's a 'y'
            elseif any(strcmp(g_oridStruct.(erqs{count}).refl.(stationz{count2}).status,'y'))~=0 && numel(g_oridStruct.(erqs{count}).refl.(stationz{count2}).time) > 1 %three
                scatterm(latz,lonz,30,'o','r','f') %multiple, and only one has to be a 'y'
            elseif any(strcmp(g_oridStruct.(erqs{count}).refl.(stationz{count2}).status,'m'))~=0 %maybe
                scatterm(latz,lonz,30,'o','b','f')
            elseif any(strcmp(g_oridStruct.(erqs{count}).refl.(stationz{count2}).status,'n'))~=0 %something wrong with the waveform
            end
        else %no
            scatterm(latz,lonz,30,'o','k','f')
        end
        clear latz; clear lonz;
    end
end

%% Saving file
hold off
directory = '/home/a/akfarrell/Uturuncu/Phase/';
filename = sprintf('ALLinf_points_map_%deqs_%s.png',numel(erqs),opt);
filename_wPath = fullfile(directory,filename);
hgexport(p, filename_wPath, hgexport('factorystyle'), 'Format', 'png');
toc