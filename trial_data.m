clear all; close all; clc
p = load('trial_data.mat');

%% ---------- Explanation of variables in trial_data.mat ------------------
% dtac -> matrix input for polar_coherency.m
% ch -> cell array of channel names
% az -> expected azimuth given location of earthquake and station
% inc -> expected incidence angle given location of earthquake and station
% W_all -> waveform object for this station
% wndo -> length of window; input for polar_coherency.m

%---Just to compare with the values from this code, I included:----------
% Haystack_data -> a structure, indexed by channel, of seismic data
% azim -> azimuth output from when I ran polar_coherency.m on my computer
% incd -> incidence angle output from when I ran polar_coherency.m
% ellip -> ellipticity output from when I ran polar_coherency.m

%% -------------- Shows how I came up with dtac ------------------
% Get Haystack_data from waveform object
disp('Channel names are:')
disp(p.ch)

for count3 = 1:numel(p.ch)
    Haystack_data.(p.ch{count3}) = get(p.W_all(count3),'data');
end

% Get dtac from Haystack_data
dtac = [Haystack_data.(p.ch{1})'; Haystack_data.(p.ch{3})'; Haystack_data.(p.ch{2})'];


%% ---------------- Actually Running polar_coherency ----------------------
% Get azimuth, incidence angle, and ellipticity
[azim, incd, ellip] = polar_coherency(dtac, p.wndo);
for cz = 1:numel(azim)
    if azim(cz)<0
        azim(cz) = -azim(cz)+180;
    end
end