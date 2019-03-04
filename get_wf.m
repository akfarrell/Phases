function [w_raw,EvS,stations_inEq] = get_wf(Ev,evidStruct)
EvS= sprintf('eq_%i',Ev)
stations_inEq = unique(evidStruct.(EvS).sta)

%ds = datasource('antelope', '/raid/data/antelope/databases/PLUTONS/dbmerged');
%ds = datasource('antelope', '/home/a/akfarrell/Uturuncu/db_heather/alex/dbmerged');
%ds = datasource('antelope', '/raid/data/PLUTONS/antelope/heather/dbmerged');
ds = datasource('antelope', '/Users/alexandrafarrell/Desktop/akfarrell/heather/dbmerged');
scnl = scnlobject(stations_inEq, 'HH*', 'PL');
inds_phase_p=find(strcmp(evidStruct.(EvS).phase, 'P'))
stime = min(evidStruct.(EvS).time_phase(inds_phase_p))-datenum(0,0,0,0,0,1)%0.000005%(0,0,0,0,0,3)
if any(strcmp(evidStruct.(EvS).phase, 'S')) %(0,0,0,0,0,5) below
    etime = max((max(evidStruct.(EvS).time_phase(setdiff(1:end,inds_phase_p)))+datenum(0,0,0,0,0,1)),max(evidStruct.(EvS).time_phase(inds_phase_p))+datenum(0,0,0,0,0,15))
else
    etime = max(evidStruct.(EvS).time_phase(inds_phase_p))+datenum(0,0,0,0,0,10) %(0,0,0,0,0,15)
end
w_raw = waveform(ds, scnl, stime, etime);
end