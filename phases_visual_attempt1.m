load('phaseStruct.mat')
eq = 'eq_2166';
stations = fieldnames(phaseStruct.(eq));
%% load in P and S arrivals and distance if not already done
% for count = 1:numel(stations) 
%     phaseStruct.(eq).(stations{count}).P_arr = p_arree(count);
%     phaseStruct.(eq).(stations{count}).S_arr = s_arree(count); %no S is 1
%     phaseStruct.(eq).(stations{count}).dist = hypee(count); %diagonaldistance
% end

%% Get data and plot it
