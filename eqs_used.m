%% Loading in inflection point data
tic
close all
[evidStruct, allevids] = get_eq_info();
load('failz_final.mat')
load('siteStruct.mat')
namez = fieldnames(evidStruct);
good_evids = zeros(1,1);
discard_evids = zeros(1,1);

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
scatterm(-22.27, -67.18, 100,'^','r','f')
scatterm(siteStruct.lat,siteStruct.lon,'*','r')
textm(siteStruct.lat,siteStruct.lon,siteStruct.sta)

%% Define outer polygon to see if eq's in network or not
stnz = {'PLHS','PLCL','PL03','PLTP','PLRR','PLWB','PL07','PLAR','PLSP','PLAN','PLCO','PLSQ','PLRV','PLLC','PLJR','PLLB','PLHS'};
for hurr = 1:numel(stnz)
    stnz_inds(hurr) = find(strcmp(siteStruct.sta,stnz{hurr}));
end
stn_lon = siteStruct.lon(stnz_inds);
stn_lat = siteStruct.lat(stnz_inds);
plotm(stn_lat,stn_lon);
%%
for count = 1:numel(namez)
    name_val = str2double(strrep(namez{count},'eq_',''));
%     if any(name_val==failz)
%     else
        if inpolygon(evidStruct.(namez{count}).lon(1),evidStruct.(namez{count}).lat(1),stn_lon,stn_lat)
            good_evids(numel(good_evids)+1) = name_val; %indexing is at +1
            scatterm(evidStruct.(namez{count}).lat(1), evidStruct.(namez{count}).lon(1), 50,'o','k')
        else
            discard_evids(numel(discard_evids)+1) = name_val;
            scatterm(evidStruct.(namez{count}).lat(1), evidStruct.(namez{count}).lon(1), 50,'o','c')
        end
%     end
end

%% Saving file
hold off
directory = '/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase/';
filename = 'ALLeqs_used_map.png';
filename_wPath = fullfile(directory,filename);
hgexport(p, filename_wPath, hgexport('factorystyle'), 'Format', 'png');

%% Refine which orids are used and make text file of discard orids
good_evids = good_evids(2:end);
discard_evids = discard_evids(2:end);
save('good_evids.mat','good_evids')
save('discard_evids.mat','discard_evids')
fileID = fopen('discard_evids.txt','w');
fprintf(fileID,'%d\n',discard_evids);
fclose(fileID);
toc