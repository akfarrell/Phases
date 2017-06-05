%% Loading in inflection point data
tic
close all
load('oridStruct.mat')
load('failz_final.mat')
load('siteStruct.mat')
load('good_orids.mat');

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
utuLat = -22.265;
utuLon = -67.185;

inflLat = -22.297;
inflLon = -67.203;%from Pritchard and Simons 2002 table 3
scatterm(utuLat, utuLon, 100,'^','r','f')
scatterm(inflLat, inflLon, 100, 'r','c','f')
scatterm(siteStruct.lat,siteStruct.lon,'*','r')
textm(siteStruct.lat,siteStruct.lon,siteStruct.sta)

%% Plot earthquakes
for count = 1:20%numel(good_orids)
    scatterm(oridStruct.(sprintf('eq_%d',good_orids(count))).lat(1), oridStruct.(sprintf('eq_%d',good_orids(count))).lon(1), 50,'s','k')
end

%% Saving file
hold off
directory = '/home/a/akfarrell/Uturuncu/Phase/';
filename = 'eqs_mapview_map.png';
filename_wPath = fullfile(directory,filename);
hgexport(p, filename_wPath, hgexport('factorystyle'), 'Format', 'png');

toc