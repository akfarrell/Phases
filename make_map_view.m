%orid2166 is origin 203
%% RUN APMB_schematic BEFORE THIS!!!
tic
close all
for count=203%1:numel(allorids)
    directory = '/home/a/akfarrell/Uturuncu/Phase/wf_objs';
    filename = sprintf('wf_%d.mat',allorids(count));
    filename_wPath = fullfile(directory,filename);
    if exist(filename_wPath,'file')
        load(filename_wPath)
    else
        %create and clean waveform object
        [w_raw,OrS,stations_inEq] = get_wf(allorids(count),oridStruct);
        w_clean = waveform_clean(w_raw, filterobject('b', fil, 2));
        save(filename_wPath,'w_clean', 'OrS', 'stations_inEq');
    end
    stationz = get(w_clean,'station');
    
    %% plot 
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
    scatterm(oridStruct.eq_2166.lat,oridStruct.eq_2166.lon,100,'d','filled','r')
    scatterm(-22.27, -67.18, 100,'^','r')
    site_subset_phase = {'PL03','PLLL','PLMD','PLSP'};
    site_subset_3phase = {'PLDK','PLJR','PLMK','PLMN'};
    site_subset_maybe = {'PLSM'};
    time_vals = {'', '0.51s or 1.23s', '0.47s or 0.79s', '0.45s or 0.72s'};
    %% some more plotting and such
    for count2 = 1:3:numel(w_clean) %43-45 is PLMN 43:3:43
       if find(intersect(find(strcmp(oridStruct.(OrS).sta,stationz{count2})),find(strcmp(oridStruct.(OrS).phase,'P'))))==1
            siteInd = find(strcmp(siteSub.sta, stationz{count2}));
            [latout,lonout] = reckon(oridStruct.eq_2166.lat(1),oridStruct.eq_2166.lon(1),inf_point.(stationz{count2})(1)/111.12,inf_point.(stationz{count2})(2));
            inf_point.(stationz{count2})(3) = latout;
            inf_point.(stationz{count2})(4) = lonout;
            if find(strcmp(site_subset_phase,stationz{count2}))
                scatterm(siteSub.lat(siteInd),siteSub.lon(siteInd),'o','filled','m')
                scatterm(latout,lonout,75,'*','m')
            elseif find(strcmp(site_subset_2phase,stationz{count2}))
                scatterm(siteSub.lat(siteInd),siteSub.lon(siteInd),'o','filled','r')
                dink = find(strcmp(site_subset_3phase,stationz{count2}));
                textm(siteSub.lat(siteInd)+0.015,siteSub.lon(siteInd)+0.01,time_vals{dink},'color','r')
                scatterm(latout,lonout,75,'*','r')
            elseif find(strcmp(site_subset_maybe,stationz{count2}))
                scatterm(siteSub.lat(siteInd),siteSub.lon(siteInd),'o','filled','b')
                scatterm(latout,lonout,75,'*','b')
            else
                scatterm(siteSub.lat(siteInd),siteSub.lon(siteInd),'o','filled','k')
                scatterm(latout,lonout,75,'*','k')
            end
            distElev.(stationz{count2}) = sqrt((siteSub.distKM(siteInd))^2 + (siteSub.elev(siteInd)+oridStruct.eq_2166.depth(1))^2);
            
            %% Changed to make texty map or less text
            distVal = sprintf('%s\n%2.2f\n%2.2f',siteSub.sta{siteInd},siteSub.distKM(siteInd),distElev.(stationz{count2}));
            %distVal = siteSub.sta{siteInd};
            if strcmp(stationz{count2},'PLCM')
                textm(siteSub.lat(siteInd)+0.03,siteSub.lon(siteInd)-0.025,distVal)
                %textm(siteSub.lat(siteInd)+0.015,siteSub.lon(siteInd)-0.025,distVal)
            else
                textm(siteSub.lat(siteInd)-0.03,siteSub.lon(siteInd)-0.025,distVal)
                %textm(siteSub.lat(siteInd)-0.015,siteSub.lon(siteInd)-0.025,distVal)
            end
            clear latout; clear lonout
        end

    end
    hold off
    directory = '/home/a/akfarrell/Uturuncu/Phase/';
    %filename = sprintf('inf_points_map_noDist_%d.png',allorids(count));
    filename = sprintf('inf_points_map_%d.png',allorids(count));
    filename_wPath = fullfile(directory,filename);
    hgexport(p, filename_wPath, hgexport('factorystyle'), 'Format', 'png');
end
toc