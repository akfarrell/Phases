% function maps_phasestuff(eq)
%%
eq = 'eq_2132';
close all;
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/scalebar_v3/scalebar')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/pathdist_v4')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu')
%run from phases_visual_attempt1.m
load('phaseStruct.mat')
load('siteStruct.mat')
[evidStruct, allevids] = get_eq_info();
h = figure;
latlim = [-22.75 -21.75]; %[southern_limit northern_limit] 
lonlim = [-67.75 -66.6]; %[western_limit eastern_limit]
set(h, 'Position', [1000 1000 1100 1000])

worldmap(latlim, lonlim);
% --- Define station characteristics ----%
stations = fieldnames(phaseStruct.(eq));
setm(gca, 'fontsize',11);
siteStaRaw = siteStruct.sta(1:numel(siteStruct.sta)); %remove Lazufre stations
siteLatRaw = siteStruct.lat(1:numel(siteStruct.sta));
siteLonRaw = siteStruct.lon(1:numel(siteStruct.sta));
for i=1:numel(siteStaRaw)
    for k = 1:numel(stations) 
        if strcmp(siteStaRaw{i}, stations{k})
            lat_sta(k) = siteLatRaw(i);
            lon_sta(k) = siteLonRaw(i);
            sta_name(k) = siteStaRaw(i);
        end
    end
end

eq_lat = evidStruct.(eq).lat(1);
eq_lon = evidStruct.(eq).lon(1);
grey = rgb('Grey');
thistle = rgb('Thistle');
colors = { 'm', 'b', 'c', 'g', 'r', 'y', grey, thistle}; %need to add more when I have more events

% ----- Make Plot ------- %            
borders = shaperead('BOL_adm0.shp', 'UseGeoCoords', true);
arg_borders = shaperead('ARG_adm0.shp', 'UseGeoCoords', true);
geoshow(borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
geoshow(arg_borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
northarrow('latitude', -21.9, 'longitude', -67.6, 'scaleratio', 1/20);
scalebar('location', 'sw', 'FontSize', 11)
hold on
scatterm(eq_lat, eq_lon,100,'b','d','f')
% colormap(jet)
% colormap(flipud(colormap))
% set(gca,'CLim',[min(),max()]);

for count = 1:numel(stations)
    %az = phaseStruct.(eq).(stations{count}).az_W_calc;
    az = phaseStruct.(eq).(stations{count}).azWest;
    hypotenuse = 0.05;
    inc = phaseStruct.(eq).(stations{count}).inc_W_calc;
    lat_az = lat_sta(count);
    lon_az = lon_sta(count);
    phase = phaseStruct.(eq).(stations{count}).phase;
    if phaseStruct.(eq).(stations{count}).az_exp_adj > 360
        az_exp = phaseStruct.(eq).(stations{count}).az_exp_adj-360
    else
        az_exp = phaseStruct.(eq).(stations{count}).az_exp_adj
    end
    if phaseStruct.(eq).(stations{count}).inc_exp > 360
        inc_exp = phaseStruct.(eq).(stations{count}).inc_exp-360;
    else
        inc_exp = phaseStruct.(eq).(stations{count}).inc_exp;
    end
%     for count2 = [1,5]%1:numel(az)%[1,5]%numel(az)
    for count2 = 1:numel(az)
        if strcmp(phase{count2},'S') || strcmp(phase{count2},'P') %Take out if don't want
            if az(count2) > 360
                az(count2) = az(count2)-360;
            end
            if az(count2)<=90
                u = hypotenuse*cosd(az(count2)); %vertical
                v = hypotenuse*sind(az(count2)); %horizontal
            elseif az(count2)>90 && az(count2)<=180
                u = -hypotenuse*sind(az(count2)-90); %vertical
                v = hypotenuse*cosd(az(count2)-90); %horizontal
            elseif az(count2)>180 && az(count2)<=270 
                u = -hypotenuse*cosd(az(count2)-180); %vertical
                v = -hypotenuse*sind(az(count2)-180); %horizontal
            elseif az(count2)>270 && az(count2)<=360
                u = hypotenuse*sind(az(count2)-270); %vertical
                v = -hypotenuse*cosd(az(count2)-270); %horizontal
            end
            if strcmp(phase{count2},'P')
                color = colors{1};
            elseif strcmp(phase{count2},'P1')
                color = colors{3};
            elseif strcmp(phase{count2},'S')
                color = colors{2};
            elseif strcmp(phase{count2},'P2')
                color = colors{4};
            elseif strcmp(phase{count2},'P3')
                color = colors{5};
            end
            quiverm(lat_az, lon_az,u, v, color) %color?
            clear u v color
            
        end 
    end
    if az_exp<=90
        u = hypotenuse*cosd(az_exp); %vertical
        v = hypotenuse*sind(az_exp); %horizontal
    elseif az_exp>90 && az_exp<=180
        u = -hypotenuse*sind(az_exp-90); %vertical
        v = hypotenuse*cosd(az_exp-90); %horizontal
    elseif az_exp>180 && az_exp<=270
        u = -hypotenuse*cosd(az_exp-180); %vertical
        v = -hypotenuse*sind(az_exp-180); %horizontal
    elseif az_exp>270 && az_exp<=360
        u = hypotenuse*sind(az_exp-270); %vertical
        v = -hypotenuse*cosd(az_exp-270); %horizontal
    end
    quiverm(lat_az, lon_az, u, v, 'k')
    textm(lat_az-0.01, lon_az+0.01,stations{count})
    legend('P exp','','P actual','','S actual')
    clear inc_exp az_exp lat_az lon_az inc az phase

end
%%
% makeTXTfiles(phaseStruct,eq,stations,'az_W_calc')
makeTXTfiles(phaseStruct,eq,stations,'azWest')
% makeTXTfiles(phaseStruct,eq,stations,'inc_W_calc')
makeTXTfiles(phaseStruct,eq,stations,'incWest')
disp('done')

%%

% tick_labels = (60:20:220)./(scale_factor*2);
% for i=1:numel(tick_labels)
%     tick_vals(i) = {sprintf('%1.2f',tick_labels(i))}
% end
% tick_vals
% % c=colorbar('eastoutside', 'TickLabels',tick_vals)
% % c.Label.String='Average Normalized Amplitude'
% % c.Label.FontSize=11;
% %title(sprintf('%s %d',id, num_eqs))
% directory = sprintf('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Synth_data/');
% if exist('num_averaged', 'var')
%     filename = sprintf('averaged_amps_%d_%s_nums.tiff',num_eqs, id);
% else
%     filename = sprintf('averaged_amps_%d_%s.tiff',num_eqs, id);
% end
% filename_wPath = fullfile(directory,filename);
% hold off
% %hgexport(h, filename_wPath, hgexport('factorystyle'), 'Format', 'pdf');
% fp = fillPage(h, 'margins', [0 0 0 0], 'papersize', [11 8.5]);
% print(h, '-dtiff', '-r400',filename_wPath);