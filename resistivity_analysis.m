%% The files are xyzd ascii and you can open them with a text editor and plot them easily with matlab.
% The columns are lat/long coords, depth in metres, and apparent resistivity in ohm-m.
% Note the model space extends far beyond where we have resolution, and we can only trust 
%     what is within the station network and to depths of <80km; ie. only the central model area, where the cells are finely spaced.

res = dlmread('../resistivity_data_comeau/s96p18_1b_model.03_wsinv3dmt_all_layers_cent.dat'); %reads the larger of the two models
lat = res(:,1); % 147436 x 1 array
lon = res(:,2);
depth = res(:,3);
resistivity = res(:,4);