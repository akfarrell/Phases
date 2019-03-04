function [evidStruct, allevids] = get_eq_info(dbpath,varargin)
if nargin <1
    dbpath = '/Users/alexandrafarrell/Desktop/akfarrell/heather/dbmerged'
end
    depth_to_reflector = 7.5;
    subset_exp = {sprintf('depth<=%1.1f',depth_to_reflector),'lat<=-21.8955 && lat>=-22.6698','lon<=-66.7865 && lon>=-67.6381'};

    %[PhaseStruct, phase, sta, chan] = loadArrivals_joinTables('/raid/home/a/akfarrell/Uturuncu/dbplutons_alex',subset_exp); %ordered by time
    %[PhaseStruct, phase, sta, chan] = loadArrivals_joinTables('/home/a/akfarrell/Uturuncu/db_heather/alex/dbmerged',subset_exp); %ordered by time
    if nargin < 1 || nargin == 1
        [PhaseStruct, phase, sta, chan, auth, fm] = loadArrivals_joinTables(dbpath,subset_exp); %ordered by time
    else
        [PhaseStruct, phase, sta, chan, auth, fm] = loadArrivals_joinTables(dbpath,subset_exp,'yes'); %ordered by time
    end
    %429 events, 6679 rows
    PhaseStruct.phase = phase;
    PhaseStruct.sta = sta;
    PhaseStruct.chan = chan;
    PhaseStruct.auth = auth;
    PhaseStruct.fm = fm;
    PhaseStruct.time_phase = epoch2datenum(PhaseStruct.time_phase);
    PhaseStruct.time_origin = epoch2datenum(PhaseStruct.time_origin);
    if isfield(PhaseStruct,'time_model')
    PhaseStruct.time_model = epoch2datenum(PhaseStruct.time_model);
    end

    allevids = unique(PhaseStruct.evid);
    names=fields(PhaseStruct);
    for count=1:numel(allevids)
        thisevid = allevids(count);
        evidS=sprintf('eq_%i',thisevid);
        indexes = find(PhaseStruct.evid == thisevid);
        for ncount=1:numel(names)
            evidStruct.(evidS).(names{ncount}) = PhaseStruct.(names{ncount})(indexes);
        end
    end
end