function [oridStruct, allorids] = get_eq_info()
    depth_to_reflector = 7.5;
    subset_exp = {sprintf('depth<=%1.1f',depth_to_reflector),'lat<=-21.8955 && lat>=-22.6698','lon<=-66.7865 && lon>=-67.6381'};

    [PhaseStruct, phase, sta, chan] = loadArrivals_joinTables('/raid/home/a/akfarrell/Uturuncu/dbplutons_alex',subset_exp); %ordered by time
    %429 events, 6679 rows
    PhaseStruct.phase = phase;
    PhaseStruct.sta = sta;
    PhaseStruct.chan = chan;
    PhaseStruct.time_phase = epoch2datenum(PhaseStruct.time_phase);
    PhaseStruct.time_origin = epoch2datenum(PhaseStruct.time_origin);

    allorids = unique(PhaseStruct.orid);
    names=fields(PhaseStruct);
    for count=1:numel(allorids)
        thisorid = allorids(count);
        oridS=sprintf('eq_%i',thisorid);
        indexes = find(PhaseStruct.orid == thisorid);
        for ncount=1:numel(names)
            oridStruct.(oridS).(names{ncount}) = PhaseStruct.(names{ncount})(indexes);
        end
    end
end