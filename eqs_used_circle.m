%% Loading in inflection point data
tic
close all
load('oridStruct.mat')
load('failz_final.mat')
load('siteStruct.mat')
namez = fieldnames(oridStruct);
used_orids = zeros(1,1);
used_inds = zeros(1,1);

methods = 'circlez';
ctr = 1998;
r = 0.05;

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
th = 0:pi/50:2*pi;
xunit = r*cos(th)+oridStruct.(sprintf('eq_%d',ctr)).lat(1);
yunit = r*sin(th)+oridStruct.(sprintf('eq_%d',ctr)).lon(1);
plotm(xunit,yunit)

%%
for count = 1:numel(namez)
    name_val = str2double(strrep(namez{count},'eq_',''));
    if any(name_val==failz)
    else
        if inpolygon(oridStruct.(namez{count}).lon(1),oridStruct.(namez{count}).lat(1),yunit,xunit)
            used_orids(numel(used_orids)+1) = name_val; %indexing is at +1
            scatterm(oridStruct.(namez{count}).lat(1), oridStruct.(namez{count}).lon(1), 50,'o','k')
            used_inds(numel(used_inds)+1) = find(name_val==good_orids);
        else
            scatterm(oridStruct.(namez{count}).lat(1), oridStruct.(namez{count}).lon(1), 50,'o','c')
        end
    end
end

%% Saving file
hold off
directory = '/home/a/akfarrell/Uturuncu/Phase/eqs_used';
filename = sprintf('eqs_used_map_%d_%1.3f.png',ctr,r);
filename_wPath = fullfile(directory,filename);
hgexport(p, filename_wPath, hgexport('factorystyle'), 'Format', 'png');

%% Refine which orids are used and make text file of discard orids
used_orids = used_orids(2:end);
used_inds = used_inds(2:end);
save('/home/a/akfarrell/Uturuncu/Phase/eqs_used/used_orids.mat','used_orids')
save('/home/a/akfarrell/Uturuncu/Phase/eqs_used/used_inds.mat','used_inds')
toc