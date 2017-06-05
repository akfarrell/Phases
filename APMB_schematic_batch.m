%% Load topographic data and close previous figure
load('topo_data.mat');
%load('inf_point_for2166.mat')
[lon, lat] = meshgrid(topo_data.lon, topo_data.lat);
topo_data.z = double(topo_data.z)./1000;
close
load('good_orids.mat')
for count = 1:14%numel(good_orids)
    eq = good_orids(count);
    OrS = sprintf('eq_%d',good_orids(count));
    OrS
    figz = 'no';
    siteSub = get_station_info(oridStruct, good_orids(count));

    %% Set up dat loop for all stations
    for sta_ind = 1:numel(siteSub.sta)
        if find(intersect(find(strcmp(oridStruct.(OrS).sta,siteSub.sta{sta_ind})),find(strcmp(oridStruct.(OrS).phase,'P'))))==1
            site = siteSub.sta{sta_ind};
            siteInd = find(strcmp(siteSub.sta,site));
            APMB_dep = -7.5;
            eq_depth = -oridStruct.(OrS).depth(1);

            %% APMB Boundary
            if strcmp(figz,'yes')
                h = figure();
                line([-2 siteSub.distKM(siteInd)+2],[APMB_dep,APMB_dep],'Color','r','LineWidth',2);

                %% Plot eq and station
                hold on
                scatter(0,eq_depth,200,'*','r')
                text(-0.7,eq_depth+0.75,strrep(OrS,'eq_',''))
                scatter(siteSub.distKM(siteInd),siteSub.elev(siteInd),200,'*','k')
                text(siteSub.distKM(siteInd)-0.7,siteSub.elev(siteInd)+0.75,site)
            end

            %% Interp some stuff
            lat_inds = linspace(oridStruct.(OrS).lat(1),siteSub.lat(siteInd)); %100 evenly spaced samples between
            lon_inds = linspace(oridStruct.(OrS).lon(1),siteSub.lon(siteInd)); %100 evenly spaced samples between

            elev_vals = zeros(1,numel(lon_inds)); %define sizes to keep box green :)
            distance_vals = elev_vals;
            azimuth_vals = elev_vals;

            if numel(lat_inds) == numel(lon_inds)
                for i = 1:numel(lat_inds)
                    if strcmp(figz,'yes')
                        elev_vals(i) = interp2(lon,lat,topo_data.z,lon_inds(i),lat_inds(i));
                    end
                    [distance_vals(i),azimuth_vals(i)] = distance(oridStruct.(OrS).lat(1),oridStruct.(OrS).lon(1),lat_inds(i),lon_inds(i));
                end
            end

                %% Plot interped stuff
            if strcmp(figz,'yes')
                plot(distance_vals.*111.12,elev_vals,'k')
                xlim([-2 siteSub.distKM(siteInd)+2])
                axis equal
            end

            %% Run to find reflection points from Jochen's reflection code
            result = perl('tt_reflected.pl', fullfile('velmodels2',sprintf('velmodel_%s.txt',site)), num2str(siteSub.distKM(siteInd)), num2str(siteSub.elev(siteInd)-eq_depth), 'P' ,'3','1');
            valz = strfind(result, 'path');
            travel_times=zeros(1,numel(valz)-1);
            trav_timez = strfind(result,' -- travel time');
            result
            for ind = 1:numel(strfind(string(result), 'path'))-1 %path occurs layers+1 time, first is before results
                travel_times(ind) = str2double(result(valz(ind+1)+5:trav_timez(ind)-1));
            end
            clear ind
            if strfind(result, 'Source is in layer 3') %Might need to change if velmodel is changed
                inf_dist = travel_times(1);
            else
                inf_val = find(valz>=strfind(result, 'twice'),1)-2; %REVISEEEE!!!
                inf_dist = sum(travel_times(1:inf_val-1))+travel_times(inf_val)/2;
            end
            inf_point.(OrS).(site).infl = [inf_dist mean(azimuth_vals(2:end))]; %distance AND azimuth
            [latout,lonout] = reckon(oridStruct.(OrS).lat(1),oridStruct.(OrS).lon(1),inf_point.(OrS).(site).infl(1)/111.12,inf_point.(OrS).(site).infl(2));
            inf_point.(OrS).(site).infl(3) = latout;
            inf_point.(OrS).(site).infl(4) = lonout;

            %% Plot reflection points from Jochen's reflection code
            if strcmp(figz,'yes')
                scatter(inf_point.(OrS).(site).infl(1),APMB_dep,'^','b')
                text(inf_point.(OrS).(site).infl(1)-0.75,APMB_dep+0.75,sprintf('%2.3f',inf_point.(OrS).(site).infl(1)))
                line([0 inf_point.(OrS).(site).infl(1)],[eq_depth APMB_dep])
                line([inf_point.(OrS).(site).infl(1) siteSub.distKM(siteInd)],[APMB_dep siteSub.elev(siteInd)])


                %% Save figure
                hold off
                directory = sprintf('/home/a/akfarrell/Uturuncu/Phase/ref_pts/%s',OrS);
                if ~exist(directory,'dir')
                    mkdir(directory)
                end
                filename = sprintf('%s_APMB_%d_%s.png',site,APMB_dep);
                filename_wPath = fullfile(directory,filename);
                hgexport(h, filename_wPath, hgexport('factorystyle'), 'Format', 'png');
            end

            %% Set up for loop
            clear h
        end
    end
end
save('inf_point.mat','inf_point')