%% Load topographic data and close previous figure
load('topo_data.mat');
[lon, lat] = meshgrid(topo_data.lon, topo_data.lat);
topo_data.z = double(topo_data.z)./1000;
close

%% Set up dat loop for all stations
for sta_ind = 1:numel(siteSub.sta)
    site = siteSub.sta{sta_ind};
    siteInd = find(strcmp(siteSub.sta,site));
    APMB_dep = -7.5;
    eq_depth = -oridStruct.eq_2166.depth(1);
    
    %% APMB Boundary
    h = figure();
    line([-2 siteSub.distKM(siteInd)+2],[APMB_dep,APMB_dep],'Color','r','LineWidth',2);

    %% Plot eq and station
    hold on
    scatter(0,eq_depth,200,'*','r')
    text(-0.7,eq_depth+0.75,'2166')
    scatter(siteSub.distKM(siteInd),siteSub.elev(siteInd),200,'*','k')
    text(siteSub.distKM(siteInd)-0.7,siteSub.elev(siteInd)+0.75,site)

    %% Interp some stuff
    startlat = min(siteSub.lat(siteInd),oridStruct.eq_2166.lat(1));
    stoplat = max(siteSub.lat(siteInd),oridStruct.eq_2166.lat(1));
    startlon = min(siteSub.lon(siteInd),oridStruct.eq_2166.lon(1));
    stoplon = max(siteSub.lon(siteInd),oridStruct.eq_2166.lon(1));

    lat_inds = linspace(oridStruct.eq_2166.lat(1),siteSub.lat(siteInd)); %100 evenly spaced samples between
    lon_inds = linspace(oridStruct.eq_2166.lon(1),siteSub.lon(siteInd)); %100 evenly spaced samples between

    elev_vals = zeros(1,numel(lon_inds)); %define sizes to keep box green :)
    distance_vals = elev_vals;
    azimuth_vals = elev_vals;

    if numel(lat_inds) == numel(lon_inds)
        for i = 1:numel(lat_inds)
            elev_vals(i) = interp2(lon,lat,topo_data.z,lon_inds(i),lat_inds(i));
            [distance_vals(i),azimuth_vals(i)] = distance(oridStruct.eq_2166.lat(1),oridStruct.eq_2166.lon(1),lat_inds(i),lon_inds(i));
        end
    end

    %% Plot interped stuff
    plot(distance_vals.*111.12,elev_vals,'k')
    xlim([-2 siteSub.distKM(siteInd)+2])
    axis equal

    %% Run to find reflection points from Jochen's reflection code
    result = perl('tt_reflected.pl', fullfile('velmodels2',sprintf('velmodel_%s.txt',site)), num2str(siteSub.distKM(siteInd)), num2str(siteSub.elev(siteInd)-eq_depth), 'P' ,'3','1');
    valz = strfind(result, 'path');
    travel_times=zeros(1,numel(valz)-1);
    trav_timez = strfind(result,' -- travel time');
    for ind = 1:numel(strfind(string(result), 'path'))-1 %path occurs layers+1 time, first is before results
        travel_times(ind) = str2double(result(valz(ind+1)+5:trav_timez(ind)-1));
    end
    inf_val = find(valz>=strfind(result, 'twice'),1)-2;
    clear ind
    inf_dist = sum(travel_times(1:inf_val-1))+travel_times(inf_val)/2;
    inf_point.(site) = [inf_dist mean(azimuth_vals(2:end))]; %distance AND azimuth

    %% Plot reflection points from Jochen's reflection code
    scatter(inf_point.(site)(1),APMB_dep,'^','b')
    text(inf_point.(site)(1)-0.75,APMB_dep+0.75,sprintf('%2.3f',inf_point.(site)(1)))
    line([0 inf_point.(site)(1)],[eq_depth APMB_dep])
    line([inf_point.(site)(1) siteSub.distKM(siteInd)],[APMB_dep siteSub.elev(siteInd)])

    %% Save figure
    hold off
    directory = '/home/a/akfarrell/Uturuncu/Phase/ref_pts/';
    filename = sprintf('%s_APMB_%d.png',site,APMB_dep);
    filename_wPath = fullfile(directory,filename);
    hgexport(h, filename_wPath, hgexport('factorystyle'), 'Format', 'png');

    %% Set up for loop
    clear h
end