%% load stuff and determine distances
addpath('WindRose/')
addpath(genpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase'))
addpath('/Users/alexandrafarrell/Desktop/akfarrell/')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu')
[evidStruct, allorids] = get_eq_info();
load('good_evids.mat')
utuLat = -22.265;
utuLon = -67.185;

inflLat = -22.297;
inflLon = -67.203;%from Pritchard and Simons 2002 table 3
dist = zeros(1,numel(good_evids));
az = dist; dist2 = dist; az2 = dist;
for count = 1:numel(good_evids)
    ind = find(good_evids(count)==allorids);
    lat = evidStruct.(sprintf('eq_%d',good_evids(count))).lat(1);
    lon = evidStruct.(sprintf('eq_%d',good_evids(count))).lon(1);
    [dist(count),az(count)] = distance(utuLat,utuLon,lat,lon);
    [dist2(count),az2(count)] = distance(inflLat,inflLon,lat,lon);
end
dist = dist.*111.12;
%az = az.*(pi/180);

dist2 = dist2.*111.12;
%az2 = az2.*(pi/180);

%% Form histogram for dist from volcano
close all
%edgez = 0:2:38;
edgez = 0:3:42;
N = histc(dist,edgez);
a = figure();
burr = bar(edgez,N,0.8,'histc','k');
set(burr,'FaceColor','k');
xlim([0 43]);
%xlim([0 39]);
xlabel('Distance from Volcano (km)','FontSize',18)
ylabel('Number of Events','FontSize',18)
xlabs={'0','3','6','9','12','15','18','21','24','27','30','33','36','39','42'};
%xlabs={'1','3','5','7','9','11','13','15','17','19','21','23','25','27','29','31','33','35','37'};
%xlabs={'0','2','4','6','8','10','12','14','16','18','20','22','24','26','28','30','32','34','36','38'};
set(gca, 'xTickLabel', xlabs);
%set(gca, 'xTick', edgez);
%set(gca, 'xTick', 1:2:37);
set(gca, 'xTick', 0:3:42);
% title('Histogram of distance from volcano')
hold off

directory = '/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase/';
filename = 'hist_dist_from_volc.png';
filename_wPath = fullfile(directory,filename);
hgexport(a, filename_wPath, hgexport('factorystyle'), 'Format', 'png');

%% Form histogram for dist from inflation center
%edgez = 0:2:42;
edgez = 0:3:42;
N = histc(dist2,edgez);
b = figure();
burr2=bar(edgez,N,0.8,'histc','k');
set(burr2,'FaceColor','k');
xlim([0 43]);
xlabel('Distance from Center of Inflation (km)','FontSize',18)
ylabel('Number of Events','FontSize',18)
xlabs={'0','3','6','9','12','15','18','21','24','27','30','33','36','39','42'};
%xlabs={'1','3','5','7','9','11','13','15','17','19','21','23','25','27','29','31','33','35','37','39','41'};
%xlabs={'0','2','4','6','8','10','12','14','16','18','20','22','24','26','28','30','32','34','36','38'};
set(gca, 'xTickLabel', xlabs);
%set(gca, 'xTick', edgez);
%set(gca, 'xTick', 1:2:41);
set(gca, 'xTick', 0:3:42);
% title('Histogram of distance from center of inflation')
hold off

directory = '/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase/';
filename = 'hist_dist_from_infl.png';
filename_wPath = fullfile(directory,filename);
hgexport(b, filename_wPath, hgexport('factorystyle'), 'Format', 'png');

%% Form rose diagram for az from volcano
deg_labels = {sprintf('0%s',char(176)), sprintf('180%s',char(176)), sprintf('90%s',char(176)), sprintf('270%s',char(176))};
intensity = ones(1,numel(az));
intensity(1:end) = 0.2;
c= WindRose(az,intensity,18,'nspeeds',2,'colors',[0.6 0.6 0.6; 0.6 0.6 0.6],'titlestring','',...
    'labels',deg_labels,'ndirections',20,'freqlabelangle',45,'legendtype',0);
%set(gca,'LineWidth',2,'Color',[0 0 0],'FontSize',18)
%set(q,'LineWidth',2,'Color','k')
%g = patch(get(q,'XData'),get(q,'YData'),'blue');
% title('Azimuths from volcano to events')

directory = '/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase/';
filename = 'rose_az_from_volc.png';
filename_wPath = fullfile(directory,filename);
hgexport(c, filename_wPath, hgexport('factorystyle'), 'Format', 'png');

%% Form rose diagram for az from inflation center
d = WindRose(az2,intensity,18,'nspeeds',2,'colors',[0.6 0.6 0.6; 0.6 0.6 0.6],'titlestring','',...
    'labels',deg_labels,'ndirections',20,'freqlabelangle',45,'legendtype',0);
%set(gca,'LineWidth',1,'Color',[0 0 0])
%set(e,'LineWidth',2,'Color','k')
% x = get(e,'Xdata');
% y = get(e,'Ydata');
% p = patch(x,y,'blue');
% title('Azimuths from inflation center to events')

directory = '/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase/';
filename = 'rose_az_from_infl.png';
filename_wPath = fullfile(directory,filename);
hgexport(d, filename_wPath, hgexport('factorystyle'), 'Format', 'png');
