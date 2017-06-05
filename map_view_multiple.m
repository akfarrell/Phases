%% Loading in inflection point data
tic
close all
clear inf_point %@!!!!!!!!!!-----make sure that inf_point is saved if there are changes-------!!!!!
load('inf_point_2166and2034.mat')
load('oridStruct.mat')
erqs = fieldnames(inf_point);

%% Plot 
p=figure; hold on;
set(p, 'Position', [1000 1000 1200 1200])
latlim = [-22.75 -21.85]; %[southern_limit northern_limit] 
lonlim = [-67.75 -66.75]; %[western_limit eastern_limit]
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
for count=1:numel(erqs)
    stationz = fieldnames(inf_point.(erqs{count}));
    scatterm(oridStruct.(erqs{count}).lat,oridStruct.(erqs{count}).lon,100,'d','filled','r')
    for count2 = 1:numel(stationz)
        latz = inf_point.(erqs{count}).(stationz{count2}).infl(3);
        lonz = inf_point.(erqs{count}).(stationz{count2}).infl(4);
        if strcmp(inf_point.(erqs{count}).(stationz{count2}).phase,'y') %yes
              scatterm(latz,lonz,75,'*','m')
        elseif strcmp(inf_point.(erqs{count}).(stationz{count2}).phase,'t') %three
            scatterm(latz,lonz,75,'*','r')
        elseif strcmp(inf_point.(erqs{count}).(stationz{count2}).phase,'m') %maybe
            scatterm(latz,lonz,75,'*','b')
        else %no
            scatterm(latz,lonz,75,'*','k')
        end
        clear latz; clear lonz;
    end
end

%% Saving file
hold off
directory = '/home/a/akfarrell/Uturuncu/Phase/';
filename = sprintf('ALLinf_points_map_%deqs.png',numel(erqs));
filename_wPath = fullfile(directory,filename);
hgexport(p, filename_wPath, hgexport('factorystyle'), 'Format', 'png');
toc