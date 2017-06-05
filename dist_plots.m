%% load stuff and determine distances
[oridStruct, allorids] = get_eq_info();
load('good_orids.mat')
utuLat = -22.265;
utuLon = -67.185;

inflLat = -22.297;
inflLon = -67.203;%from Pritchard and Simons 2002 table 3
dist = zeros(1,numel(good_orids));
az = dist; dist2 = dist; az2 = dist;
for count = 1:numel(good_orids)
    ind = find(good_orids(count)==allorids);
    lat = oridStruct.(sprintf('eq_%d',good_orids(count))).lat(1);
    lon = oridStruct.(sprintf('eq_%d',good_orids(count))).lon(1);
    [dist(count),az(count)] = distance(utuLat,utuLon,lat,lon);
    [dist2(count),az2(count)] = distance(inflLat,inflLon,lat,lon);
end
dist = dist.*111.12;
az = az.*(pi/180);

dist2 = dist2.*111.12;
az2 = az2.*(pi/180);

%% Form histogram for dist from volcano
close all
edgez = 0:2:38;
N = histc(dist,edgez);
a = figure();
bar(edgez,N,0.8,'histc')
xlim([0 38]);
xlabel('Distance from Volcano (km)')
ylabel('Number of Events')
xlabs={'1','3','5','7','9','11','13','15','17','19','21','23','25','27','29','31','33','35','37'};
%xlabs={'0','2','4','6','8','10','12','14','16','18','20','22','24','26','28','30','32','34','36','38'};
set(gca, 'xTickLabel', xlabs);
%set(gca, 'xTick', edgez);
set(gca, 'xTick', 1:2:37);
title('Histogram of distance from volcano')
hold off

directory = '/home/a/akfarrell/Uturuncu/Phase/';
filename = 'hist_dist_from_volc.png';
filename_wPath = fullfile(directory,filename);
hgexport(a, filename_wPath, hgexport('factorystyle'), 'Format', 'png');

%% Form histogram for dist from inflation center
edgez = 0:2:42;
N = histc(dist2,edgez);
b = figure();
bar(edgez,N,0.8,'histc')
xlim([0 42]);
xlabel('Distance from Center of Inflation (km)')
ylabel('Number of Events')
xlabs={'1','3','5','7','9','11','13','15','17','19','21','23','25','27','29','31','33','35','37','39','41'};
%xlabs={'0','2','4','6','8','10','12','14','16','18','20','22','24','26','28','30','32','34','36','38'};
set(gca, 'xTickLabel', xlabs);
%set(gca, 'xTick', edgez);
set(gca, 'xTick', 1:2:41);
title('Histogram of distance from center of inflation')
hold off

directory = '/home/a/akfarrell/Uturuncu/Phase/';
filename = 'hist_dist_from_infl.png';
filename_wPath = fullfile(directory,filename);
hgexport(b, filename_wPath, hgexport('factorystyle'), 'Format', 'png');

%% Form rose diagram for az from volcano
c = figure();
rose(az)
title('Azimuths from volcano to events')

directory = '/home/a/akfarrell/Uturuncu/Phase/';
filename = 'rose_az_from_volc.png';
filename_wPath = fullfile(directory,filename);
hgexport(c, filename_wPath, hgexport('factorystyle'), 'Format', 'png');


%% Form rose diagram for az from inflation center
d = figure();
rose(az2)
title('Azimuths from inflation center to events')

directory = '/home/a/akfarrell/Uturuncu/Phase/';
filename = 'rose_az_from_infl.png';
filename_wPath = fullfile(directory,filename);
hgexport(d, filename_wPath, hgexport('factorystyle'), 'Format', 'png');
