%% Get data from figure generated in time_plots.m
%close all
r = openfig('velocity_vs_dist_circle');
ch=get(gca,'ch');
marker_shape = cell(1,numel(ch));
for i = 1:numel(ch)
    marker_shape{i} = ch(i).Marker;
end
plus_ind = find(strcmp('+',marker_shape)); %Direct P
square_ind = find(strcmp('square',marker_shape)); %maybe
diamond_ind = find(strcmp('diamond',marker_shape)); %definitely

plus_x = zeros(1,numel(plus_ind)); plus_y = plus_x;
for count = 1:numel(plus_ind)
    plus_x(count) = ch(plus_ind(count)).XData;
    plus_y(count) = ch(plus_ind(count)).YData;
end

square_x = zeros(1,numel(square_ind)); square_y = square_x;
for count2 = 1:numel(square_ind)
    square_x(count2) = ch(square_ind(count2)).XData;
    square_y(count2) = ch(square_ind(count2)).YData;
end

diamond_x = zeros(1,numel(diamond_ind)); diamond_y = diamond_x; %definitely
for count3 = 1:numel(diamond_ind)
    diamond_x(count3) = ch(diamond_ind(count3)).XData;
    diamond_y(count3) = ch(diamond_ind(count3)).YData;
end

%% Calculate velocities in km/s
vel_plus = plus_x./plus_y;
vel_square = square_x./square_y;
vel_diamond = diamond_x./diamond_y;
combo_vel = horzcat(vel_square,vel_diamond);

%% Calculate histograms for velocities
g = figure();hold on;
set(g, 'Position', [1000 1000 1200 1200])
% Direct P velocities
subplot(4,1,1)
[N,X] = hist(vel_plus,30);
bar(X,N)
title('Direct P-wave velocities (km/s)')
%maybe velocities
subplot(4,1,3)
[N2,X2] = hist(vel_square,30);
bar(X2,N2)
title('Maybe phase velocities (km/s)')
%definitely velocities
subplot(4,1,2)
[N3,X3] = hist(vel_diamond,30);
bar(X3,N3)
title('Definite phase velocities (km/s)')
%combo maybe and definitely velocities
subplot(4,1,4)
[N4,X4] = hist(combo_vel,30);
bar(X4,N4)
title('Definite and maybe phase velocities (km/s)')
samexaxis()