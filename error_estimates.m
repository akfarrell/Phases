clc; clear; close all;
evid = 2132;
%evid = 2120;
addpath('/Users/alexandrafarrell/Desktop/akfarrell/')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Polarizemic-master/functions/')
[evidStruct, allevids] = get_eq_info('/Users/alexandrafarrell/Desktop/akfarrell/heather/dbmerged'); %just relocated, no errors
[evidStruct_error, allevids_error] = get_eq_info('/Users/alexandrafarrell/Desktop/akfarrell/heather/dbmerged','true'); %relocated, with errors (excludes all m<0.5

[evidStruct_unchecked_error, allevids_unchecked_error] = get_eq_info('/Users/alexandrafarrell/Desktop/dbplutons_alex','true'); %not relocated, with errors
load('phaseStruct.mat')
time_diff = zeros(size(evidStruct_error.(sprintf('eq_%d',evid)).time_model));
model_arr = cell(size(time_diff));
phase_arr = model_arr;
for count = 1:numel(time_diff)
    time_diff(count) = etime(datevec(evidStruct_error.(sprintf('eq_%d',evid)).time_model(count)),datevec(evidStruct_error.(sprintf('eq_%d',evid)).time_phase(count)));
    model_arr{count} = datestr(evidStruct_error.(sprintf('eq_%d',evid)).time_model(count),'mmmm dd yyyy HH:MM:SS.FFF');
    phase_arr{count} = datestr(evidStruct_error.(sprintf('eq_%d',evid)).time_phase(count),'mmmm dd yyyy HH:MM:SS.FFF');
    evidStruct_error.(sprintf('eq_%d',evid)).time_diff(count) = time_diff(count);
end
%%
p_inds = find(strcmp(evidStruct_error.(sprintf('eq_%d',evid)).phase,'P'));
s_inds = find(strcmp(evidStruct_error.(sprintf('eq_%d',evid)).phase,'S'));
ps = time_diff(p_inds);
ss = time_diff(s_inds);
mnp = mean(abs(ps))
mns = mean(abs(ss))
mn_overall = mean(abs(time_diff))
PLLL_P = evidStruct_error.(sprintf('eq_%d',evid)).time_diff(p_inds(strcmp(evidStruct_error.(sprintf('eq_%d',evid)).sta(p_inds),'PLLL')))
PLSS_P = evidStruct_error.(sprintf('eq_%d',evid)).time_diff(p_inds(strcmp(evidStruct_error.(sprintf('eq_%d',evid)).sta(p_inds),'PLSS')))
PLSS_S = evidStruct_error.(sprintf('eq_%d',evid)).time_diff(s_inds(strcmp(evidStruct_error.(sprintf('eq_%d',evid)).sta(s_inds),'PLSS')))
PLLL_S = evidStruct_error.(sprintf('eq_%d',evid)).time_diff(s_inds(strcmp(evidStruct_error.(sprintf('eq_%d',evid)).sta(s_inds),'PLLL')))