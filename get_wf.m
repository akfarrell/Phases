function [w_raw,OrS,stations_inEq] = get_wf(Or,oridStruct)
OrS= sprintf('eq_%i',Or)
stations_inEq = unique(oridStruct.(OrS).sta)

ds = datasource('antelope', '/raid/data/antelope/databases/PLUTONS/dbmerged');
scnl = scnlobject(stations_inEq, 'HH*', 'PL');
inds_phase_p=find(strcmp(oridStruct.(OrS).phase, 'P'))
stime = min(oridStruct.(OrS).time_phase(inds_phase_p))-0.00003
if any(strcmp(oridStruct.(OrS).phase, 'S'))
    etime = max(oridStruct.(OrS).time_phase(setdiff(1:end,inds_phase_p)))+0.0003
else
    etime = max(oridStruct.(OrS).time_phase(inds_phase_p))+0.0003
end
w_raw = waveform(ds, scnl, stime, etime);
end