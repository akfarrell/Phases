%Define variables

a = 5.17; %P-wave velocity in km/s
b = 2.95; %S-wave velocity in km/s
eq_dep = 0; %Depth of earthquake
dist = 1; %Map distance from earthquake to station, in km
Phase = 'R'; %Phase looking for

%Output: P S PxP SxS PxS SxP
if strcmp(Phase, 'P')
elseif strcmp(Phase, 'S')
elseif strcmp(Phase, 'PxP')
elseif strcmp(Phase, 'SxS')
elseif strcmp(Phase, 'PxS')
elseif strcmp(Phase, 'SxP')
else
    fprintf('ERROR: %s not an acceptable phase\nChoose from P, S, PxP, SxS, PxS, or SxP\n', Phase)
end