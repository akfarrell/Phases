function siteSub = get_station_info(oridStruct, Or)
    [w_raw,OrS,stations_inEq] = get_wf(Or,oridStruct);


    %clean waveforms and sort by distance from earthquake
    %close all;
    fil=[2 25];
    w_clean = waveform_clean(w_raw, filterobject('b', fil, 2));

    sta=unique(get(w_clean,'station'));
    siteStruct = loadSiteTable('/raid/home/a/akfarrell/Uturuncu/dbplutons_alex','lat>=-23');
    for count=1:numel(stations_inEq)
        siteSub.sta(count) = siteStruct.sta(strcmp(stations_inEq(count), siteStruct.sta));
        siteSub.lat(count) = siteStruct.lat(strcmp(stations_inEq(count), siteStruct.sta));
        siteSub.lon(count) = siteStruct.lon(strcmp(stations_inEq(count), siteStruct.sta));
        siteSub.elev(count) = siteStruct.elev(strcmp(stations_inEq(count), siteStruct.sta));
        siteSub.dist(count) = distance(oridStruct.(OrS).lat(1), oridStruct.(OrS).lon(1), ...
        siteStruct.lat(strcmp(stations_inEq(count), siteStruct.sta)), siteStruct.lon(strcmp(stations_inEq(count), siteStruct.sta)));
    end
    clear count
    siteSub.distKM=siteSub.dist*111.12;
end