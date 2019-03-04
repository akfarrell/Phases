%% Loading in inflection point data
tic
close all
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/scalebar_v3/scalebar')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/pathdist_v4')
[evidStruct, allevids] = get_eq_info('/Users/alexandrafarrell/Desktop/akfarrell/heather/dbmerged'); %just relocated, no errors
[evidStruct_error, allevids_error] = get_eq_info('/Users/alexandrafarrell/Desktop/akfarrell/heather/dbmerged','true'); %relocated, with errors (excludes all m<0.5

[evidStruct_unchecked_error, allevids_unchecked_error] = get_eq_info('/Users/alexandrafarrell/Desktop/dbplutons_alex','true'); %not relocated, with errors
% load('failz_final.mat')
load('siteStruct.mat')
% load('allevids_unchecked_error.mat');

%% Plot 
p=figure; hold on;


eqs_interest = [2120 2121 2122 2123 2124 2125 2126 2127 2128 2130 2131 2132];
% eqs_interest = [1992 1993 1994 1995 1996 1997];
% eqs_interest = [1847 1848 1849 1850 1851 1852 1853 1854 1855 1856 1857 1858 1859 1860 1861 1862 1863 1864 1866 1867 1868 1869 1870 1871 1872 1874 1875];
% eqs_interest = [2008 2009 2010 2011 2012 2013 2014 2015 2016 2018 2019 2020 2021 2022 2023 2024 2025 2026];
% eqs_interest = [2305 2306 2307 2308 2309 2310 2311 2312];
% eqs_interest = [2409 2410 2411 2412 2413];

lat = zeros(1,numel(eqs_interest));
lon = zeros(1,numel(eqs_interest));
for count = 1:numel(eqs_interest)
    inds = find(allevids_unchecked_error==eqs_interest(count));
    lat(count) = evidStruct_unchecked_error.(sprintf('eq_%d',allevids_unchecked_error(inds))).lat(1);
    lon(count) = evidStruct_unchecked_error.(sprintf('eq_%d',allevids_unchecked_error(inds))).lon(1);
end
clear count inds

set(p, 'Position', [1000 1000 1200 1200])
latlim = [-22.75 -21.85]; %[southern_limit northern_limit] 
lonlim = [-67.75 -66.75]; %[western_limit eastern_limit]

% latlim = [min(lat)-0.05 max(lat)+0.05];
% lonlim = [min(lon)-0.05 max(lon)+0.05];

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
% eqs_interest = [2062 2068 2070 2075];



% for count = 43:58%numel(allevids_unchecked_error)
%     %plotm([oridStruct.(sprintf('eq_%d',good_orids(count))).lat(1) utuLat], [oridStruct.(sprintf('eq_%d',good_orids(count))).lon(1) utuLon])
%     scatterm(evidStruct.(sprintf('eq_%d',good_orids(count))).lat(1), evidStruct.(sprintf('eq_%d',allevids_unchecked_error(count))).lon(1), 50,'s','k')
% end
smajax = zeros(1,numel(eqs_interest));
sminax = smajax;
strike = smajax;
for count = 1:numel(eqs_interest)
    inds = find(allevids_unchecked_error==eqs_interest(count));
%     lat = evidStruct.(sprintf('eq_%d',allevids_unchecked_error(inds))).lat(1);
%     lon = evidStruct.(sprintf('eq_%d',allevids_unchecked_error(inds))).lon(1);
    textm(lat(count)-0.001, lon(count)-0.001,sprintf('%d',eqs_interest(count)))
    scatterm(lat(count), lon(count), 50,'s','f','b')
    if intersect(eqs_interest(count),allevids_error)
        smajax(count) = evidStruct_error.(sprintf('eq_%d',allevids_unchecked_error(inds))).smajax(1);
        sminax(count) = evidStruct_error.(sprintf('eq_%d',allevids_unchecked_error(inds))).sminax(1);
        strike(count) = evidStruct_error.(sprintf('eq_%d',allevids_unchecked_error(inds))).strike(1);
    elseif intersect(eqs_interest(count),allevids_unchecked_error)
%         count
        smajax(count) = evidStruct_unchecked_error.(sprintf('eq_%d',allevids_unchecked_error(inds))).smajax(1);
        sminax(count) = evidStruct_unchecked_error.(sprintf('eq_%d',allevids_unchecked_error(inds))).sminax(1);
        strike(count) = evidStruct_unchecked_error.(sprintf('eq_%d',allevids_unchecked_error(inds))).strike(1);
    end
    ecc = axes2ecc(smajax(count)/111.12,sminax(count)/111.12);
    [elat,elon] = ellipse1(lat(count), lon(count), [smajax(count)/111.12, ecc],strike(count));
    plotm(elat,elon)
    clear elat elon ecc
end

northarrow('latitude', -21.9, 'longitude', -67.6, 'scaleratio', 1/20);
scalebar('location', 'sw', 'FontSize', 11)

%% Saving file
hold off
directory = '/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase/examples';
filename = 'eqs_mapview_map_2120-2132.png';
filename_wPath = fullfile(directory,filename);
hgexport(p, filename_wPath, hgexport('factorystyle'), 'Format', 'png');

toc