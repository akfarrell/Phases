function w = prep_for_rs(w,eq_id,siteStruct,oridStruct)
%% Prep_for_rs.m - preps waveform object to work with plotw_rs_test.m
% Inputs:
%    w - waveform object
%    eq_id - earthquake ID, as a string (ex 'eq_2166')
%    siteStruct - structure of site information
%    oridStruct - structure of event information
% Output:
%    w - waveform object with added fields 'start','station','STLA','STLO',
%        'EVLA','EVLO','EVDP','KEVNM','mag','KNETWK'
%
% suggested usage: w = prep_for_rs(w,eq_id,siteStruct,oridStruct)
%    written by Alex Farrell 2/8/18

% w = w_clean;
% eq_id = 'eq_2166';
st_w = get(w, 'station');
stas = siteStruct.sta;
for count = 1:numel(w)
    ind = find(strcmp(stas,st_w(count)));
    w(count) = addfield(w(count), 'STLA', siteStruct.lat(ind));
    w(count) = addfield(w(count), 'STLO', siteStruct.lon(ind));
    w(count) = addfield(w(count), 'STEL', siteStruct.elev(ind));
    w(count) = addfield(w(count), 'EVLA', oridStruct.(eq_id).lat(1));
    w(count) = addfield(w(count), 'EVLO', oridStruct.(eq_id).lon(1));
    w(count) = addfield(w(count), 'EVDP', oridStruct.(eq_id).depth(1));
    w(count) = addfield(w(count), 'KEVNM', eq_id);
    w(count) = addfield(w(count), 'mag', oridStruct.(eq_id).ml(1));
    w(count) = addfield(w(count), 'KNETWK', 'PLUTONS');
end
% [starttime,sta,rlat,rlon,elat,elon,edep,eid,mag,netwk] = ...
%     getm(w,'start','station','STLA','STLO','EVLA','EVLO','EVDP','KEVNM','mag','KNETWK');