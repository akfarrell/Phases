%% Loading in inflection point data
tic
close all
load('oridStruct.mat')
load('failz_final.mat')
load('siteStruct.mat')
namez = fieldnames(oridStruct);
good_orids = zeros(1,1);
discard_orids = zeros(1,1);

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
stnz = {'PLHS','PLCL','PL03','PLTP','PLRR','PLWB','PLAR','PLSP','PLAN','PLCO','PLSQ','PLRV','PLLC','PLJR','PLLB','PLHS'};
for hurr = 1:numel(stnz)
    stnz_inds(hurr) = find(strcmp(siteStruct.sta,stnz{hurr}));
end
stn_lon = siteStruct.lon(stnz_inds);
stn_lat = siteStruct.lat(stnz_inds);
plotm(stn_lat,stn_lon);
%%
for count = 1:numel(namez)
    name_val = str2double(strrep(namez{count},'eq_',''));
    if any(name_val==failz)
    else
        if inpolygon(oridStruct.(namez{count}).lon(1),oridStruct.(namez{count}).lat(1),stn_lon,stn_lat)
            good_orids(numel(good_orids)+1) = name_val; %indexing is at +1
            scatterm(oridStruct.(namez{count}).lat(1), oridStruct.(namez{count}).lon(1), 50,'o','k')
        else
            discard_orids(numel(discard_orids)+1) = name_val;
            scatterm(oridStruct.(namez{count}).lat(1), oridStruct.(namez{count}).lon(1), 50,'o','c')
        end
    end
end

%% Saving file
hold off
directory = '/home/a/akfarrell/Uturuncu/Phase/';
filename = 'ALLeqs_used_map.png';
filename_wPath = fullfile(directory,filename);
hgexport(p, filename_wPath, hgexport('factorystyle'), 'Format', 'png');

%% Refine which orids are used and make text file of discard orids
good_orids = good_orids(2:end);
discard_orids = discard_orids(2:end);
save('good_orids.mat','good_orids')
save('discard_orids.mat','discard_orids')
fileID = fopen('discard_orids.txt','w');
fprintf(fileID,'%d\n',discard_orids);
fclose(fileID);
toc