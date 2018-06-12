function [w_raw,OrS,stations_inEq] = get_wf(Or,oridStruct)
OrS= sprintf('eq_%i',Or)
stations_inEq = unique(oridStruct.(OrS).sta)

ds = datasource('antelope', '/raid/data/antelope/databases/PLUTONS/dbmerged');
scnl = scnlobject(stations_inEq, 'HH*', 'PL');
inds_phase_p=find(strcmp(oridStruct.(OrS).phase, 'P'))
stime = min(oridStruct.(OrS).time_phase(inds_phase_p))-datenum(0,0,0,0,0,1)%0.000005%(0,0,0,0,0,3)
if any(strcmp(oridStruct.(OrS).phase, 'S')) %(0,0,0,0,0,5) below
    etime = max((max(oridStruct.(OrS).time_phase(setdiff(1:end,inds_phase_p)))+datenum(0,0,0,0,0,1)),max(oridStruct.(OrS).time_phase(inds_phase_p))+datenum(0,0,0,0,0,15))
else
    etime = max(oridStruct.(OrS).time_phase(inds_phase_p))+datenum(0,0,0,0,0,10) %(0,0,0,0,0,15)
end
w_raw = waveform(ds, scnl, stime, etime);
end