%% Loading in inflection point data
tic
close all
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase/CurveIntersect')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/scalebar_v3/scalebar')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/pathdist_v4')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase/reexternalre/')

addpath(genpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase'))
addpath('/Users/alexandrafarrell/Desktop/akfarrell/')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Polarizemic-master/functions/')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/fillPage/')
[evidStruct, allevids] = get_eq_info();
[evidStruct_unchecked_error, allevids_unchecked_error] = get_eq_info('/Users/alexandrafarrell/Desktop/dbplutons_alex','true'); %not relocated, with errors
[evidStruct_error, allevids_error] = get_eq_info('/Users/alexandrafarrell/Desktop/akfarrell/heather/dbmerged','true'); %relocated, with errors (excludes all m<0.5

% load('failz_final.mat')
load('siteStruct.mat')
load('orientation_struct.mat')
% load('good_evids.mat');

%% Plot 
% p=figure; hold on;
% set(p, 'Position', [1000 1000 1200 1200])
% latlim = [-22.75 -21.85]; %[southern_limit northern_limit] 
% lonlim = [-67.75 -66.75]; %[western_limit eastern_limit]
% worldmap(latlim, lonlim);
% borders = shaperead('../BOL_adm0.shp', 'UseGeoCoords', true);
% colormap(flipud(colormap))
% arg_borders = shaperead('../ARG_adm0.shp', 'UseGeoCoords', true);
% try
%     geoshow(borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
% catch
% end
% try
%     geoshow(arg_borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
% catch
% end
% hold on
% utuLat = -22.265;
% utuLon = -67.185;
% 
% inflLat = -22.297;
% inflLon = -67.203;%from Pritchard and Simons 2002 table 3
% scatterm(utuLat, utuLon, 100,'^','r','f')
% scatterm(inflLat, inflLon, 100, 'r','c','f')
% scatterm(siteStruct.lat,siteStruct.lon,'*','r')
% textm(siteStruct.lat,siteStruct.lon,siteStruct.sta)

%% Plot earthquakes
% eqs_interest = [2062 2068 2070 2075];
eqs = fieldnames(orientation_struct);
for count = 1:numel(eqs)
    eqz{count} = strrep(eqs{count},'eq_','');
end
eqs_interest = str2double(eqz);
clear count

% eqs_interest = [2501 2494 2251 2311 2312 2436 1786 1951 1864 1870 2000 2042 2043 2046 2049 2050 2054 2055 ...
%     2056 2057 2119 2152 2153 2161 2165 2166 2177 2210 2217 2223 2247 2249 2267 2272 2273 2276 2295 2313 ...
%     2314 2324 2347 2382 2386 2390];
% eqs_interest = [2120 2121 2122 2123 2124 2125 2126 2127 2128 2130 2131 2132];

% for count = 43:58%numel(good_evids)
%     %plotm([oridStruct.(sprintf('eq_%d',good_orids(count))).lat(1) utuLat], [oridStruct.(sprintf('eq_%d',good_orids(count))).lon(1) utuLon])
%     scatterm(evidStruct.(sprintf('eq_%d',good_orids(count))).lat(1), evidStruct.(sprintf('eq_%d',good_evids(count))).lon(1), 50,'s','k')
% end
smajax = zeros(1,numel(eqs_interest));
sminax = smajax;
strike = smajax;
%%
for count = 1:numel(eqs_interest)
    if intersect(eqs_interest(count),allevids)
        inds = find(allevids==eqs_interest(count));
        lat = evidStruct.(sprintf('eq_%d',allevids(inds))).lat(1);
        lon = evidStruct.(sprintf('eq_%d',allevids(inds))).lon(1);
%         textm(lat-0.001, lon-0.001,sprintf('%d',eqs_interest(count)))
%         scatterm(lat, lon, 50,'s','f','b')
    elseif intersect(eqs_interest(count),allevids_unchecked_error)
        inds = find(allevids_unchecked_error==eqs_interest(count));
        lat = evidStruct_unchecked_error.(sprintf('eq_%d',allevids_unchecked_error(inds))).lat(1);
        lon = evidStruct_unchecked_error.(sprintf('eq_%d',allevids_unchecked_error(inds))).lon(1);
%         textm(lat-0.001, lon-0.001,sprintf('%d',eqs_interest(count)))
%         scatterm(lat, lon, 50,'s','f','b')
    else
        error('not plotting all eqs')
    end
    if intersect(eqs_interest(count),allevids_error)
        smajax(count) = evidStruct_error.(sprintf('eq_%d',allevids(inds))).smajax(1)/111.12;
        sminax(count) = evidStruct_error.(sprintf('eq_%d',allevids(inds))).sminax(1)/111.12;
        strike(count) = evidStruct_error.(sprintf('eq_%d',allevids(inds))).strike(1);
    elseif intersect(eqs_interest(count),allevids_unchecked_error)
        smajax(count) = evidStruct_unchecked_error.(sprintf('eq_%d',allevids_unchecked_error(inds))).smajax(1)/111.12;
        sminax(count) = evidStruct_unchecked_error.(sprintf('eq_%d',allevids_unchecked_error(inds))).sminax(1)/111.12;
        strike(count) = evidStruct_unchecked_error.(sprintf('eq_%d',allevids_unchecked_error(inds))).strike(1);
    end
    ecc = axes2ecc(smajax(count),sminax(count));
    [elat,elon] = ellipse1(lat, lon, [smajax(count), ecc],strike(count));
    stnz = fieldnames(orientation_struct.(eqs{count}));
%     plotm(elat,elon)
    for count2 = 1:numel(stnz)
        ind = find(strcmp(siteStruct.sta,stnz{count2}));
        siteLat = siteStruct.lat(ind);
        siteLon = siteStruct.lon(ind);
%         scatterm(siteLat,siteLon)
%         textm(siteLat, siteLon, siteStruct.sta{count2})
        
        x = elat';
        y = elon';
        x0 = siteLat; y0 = siteLon;                           % given point
        
        xe = [x0, x]; ye = [y0, y];               % put all points together
        
        cvx = convhull(xe,ye);                    % find convex hull
        
        x1 = xe(cvx(2)); y1=ye(cvx(2));           % one neighbor of (x0,y0)
        x2 = xe(cvx(end-1)); y2=ye(cvx(end-1));   % another neighbor
        s = linspace(0,2,300);
        xline1 = x0+s*(x1-x0);
        xline2 = x0+s*(x2-x0);
        yline1 = y0+s*(y1-y0);
        yline2 = y0+s*(y2-y0);
%         plotm(xline1, yline1, 'r')     % plot the lines
%         plotm(xline2, yline2, 'r')
        [x1,y1] = polyxpoly(xline1, yline1, elat, elon);
        [x2,y2] = polyxpoly(xline2, yline2, elat, elon);
%         if isempty(x1)
%             [xy1,distance1,t_a1] = distance2curve(durr,hurr,'spline')
% %             [x1,y1] = curveintersect(xline1,yline1,elat,elon)
%         end
%         durr = horzcat(elat,elon);
%         hurr = horzcat(xline1',yline1');
%         if isempty(x2)
%             [xy2,distance2,t_a2] = distance2curve(durr,hurr,'spline')
% %             [x2,y2] = curvefitintersect(xline2,yline2,elat,elon)
%         end
        if isempty(x1) || isempty(x2)
            fprintf('event %d %d sta %s %d\n',eqs_interest(count), count, stnz{count2}, count2)
            continue
        end
        
        
        az1 = azimuth(siteLat,siteLon, x1, y1);
        az2 = azimuth(siteLat,siteLon, x2, y2);
        orientation_struct.(eqs{count}).(stnz{count2}).azRange = [min(az1) max(az2)];
        save('orientation_struct.mat','orientation_struct')
        clear x1 x2 xline1 xline2 yline1 yline2 y1 y2 az1 az2 x y x0 xe cvx s ind siteLat siteLon
    end
        
    clear lat lon elat elon ecc stnz smajaz sminax strike inds lat lon
end

% northarrow('latitude', -21.9, 'longitude', -67.6, 'scaleratio', 1/20);
% scalebar('location', 'sw', 'FontSize', 11)

%% Saving file
% hold off
% directory = '/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase/examples';
% filename = 'eqs_mapview_map_useful.png';
% filename_wPath = fullfile(directory,filename);
% hgexport(p, filename_wPath, hgexport('factorystyle'), 'Format', 'png');


%% Look into!!!

% event 2097 84 sta PLRR 4
% event 2160 95 sta PLCO 8
% event 2166 97 sta PLSQ 15
% event 2307 117 sta PLCO 7
% event 2314 120 sta PLSM 14
% event 2321 121 sta PLJR 7
% event 2364 126 sta PL03 1
% event 2409 135 sta PL03 1
% event 2436 141 sta PL03 1
% event 2441 143 sta PLSM 8
% event 2562 158 sta PLDK 2

toc