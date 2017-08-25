%% Loading in inflection point data
tic
close all
clear all
load('g_oridStruct.mat')
load('good_orids.mat')
erqs = fieldnames(g_oridStruct); %make sure that this is the same as the variable 'count' or problems
siteStruct = loadSiteTable('/raid/home/a/akfarrell/Uturuncu/dbplutons_alex','lat>=-23');
w_topo = 'y'; %Change to be with or without topography

%% Plotting dist vs time
p=figure; hold on;
set(p, 'Position', [1000 1000 1200 1200])
for count=1:66%numel(good_orids)
    stationz = fieldnames(g_oridStruct.(erqs{count}).refl);
    ot = g_oridStruct.(erqs{count}).time_origin(1);
    inds_phase_p = find(strcmp(g_oridStruct.(erqs{count}).phase, 'P'));
    stime = min(g_oridStruct.(erqs{count}).time_phase(inds_phase_p))-datenum(0,0,0,0,0,3);
    eq_lat = g_oridStruct.(erqs{count}).lat(1);
    eq_lon = g_oridStruct.(erqs{count}).lon(1);
    for count2 = 1:numel(stationz)
        %% caclulate the time indices, in sec*100
        ind_P = intersect(find(strcmp(g_oridStruct.(erqs{count}).sta,stationz{count2})),find(strcmp(g_oridStruct.(erqs{count}).phase,'P')));
        time_Parr = g_oridStruct.(erqs{count}).time_phase(ind_P);
        time_ot_to_P = ceil(etime(datevec(time_Parr),datevec(ot))*100);
        time_st_to_P = ceil(etime(datevec(time_Parr),datevec(stime))*100);
        P_ind = g_oridStruct.(erqs{count}).refl.(stationz{count2}).P_ind;
        
        %% Find site information, with distance from eq to site, in km
        siteInd = find(strcmp(stationz{count2},siteStruct.sta));
        siteLat = siteStruct.lat(siteInd);
        siteLon = siteStruct.lon(siteInd);
        dist = (distance(siteLat,siteLon, eq_lat, eq_lon)*111.12); %in km
        
        %% Calculate distance with depth and elevation in mind
        siteElev = siteStruct.elev(siteInd); %in km
        eqDep = g_oridStruct.(erqs{count}).depth(1); %in km
        total_elev = siteElev+eqDep;
        actual_dist = sqrt(total_elev^2+dist^2);
        if strcmp(w_topo, 'y')
            dist_backup = dist;
            dist = actual_dist;
        end
        
        %% Actually plotting stuff
        scatter(dist,time_ot_to_P/100,'+','k') %plot P wave arrival
        g_oridStruct.(erqs{count}).refl.(stationz{count2}).t_P = time_ot_to_P/100;
        if any(strcmp(fieldnames(g_oridStruct.(erqs{count}).refl.(stationz{count2})),'time'))
            for count3 = 1:numel(g_oridStruct.(erqs{count}).refl.(stationz{count2}).time)
                if strcmp(g_oridStruct.(erqs{count}).refl.(stationz{count2}).status(count3),'y')
                    t_plot = g_oridStruct.(erqs{count}).refl.(stationz{count2}).time(count3)-P_ind+time_ot_to_P;
                    g_oridStruct.(erqs{count}).refl.(stationz{count2}).t_phase(count3) = t_plot/100;
                    scatter(dist,t_plot/100,'d','f','r')
                elseif strcmp(g_oridStruct.(erqs{count}).refl.(stationz{count2}).status(count3),'m')
                    t_plot = g_oridStruct.(erqs{count}).refl.(stationz{count2}).time(count3)-P_ind+time_ot_to_P;
                    g_oridStruct.(erqs{count}).refl.(stationz{count2}).t_phase(count3) = t_plot/100;
                    scatter(dist,t_plot/100,'s','f','b')
                elseif strcmp(g_oridStruct.(erqs{count}).refl.(stationz{count2}).status(count3),'n') %if it's a 'no', do nothing
                end
                clear t_plot
            end
        end
        clear P_ind
    end
    clear ot
end
%xlim([0 30])
xlabel('Distance (km)')
ylabel('Time (s)')

%% Saving file
hold off
directory = '/home/a/akfarrell/Uturuncu/Phase/';
filename = sprintf('timeVSdist_%deqs.png',numel(erqs));
filename_wPath = fullfile(directory,filename);
hgexport(p, filename_wPath, hgexport('factorystyle'), 'Format', 'png');

savefig('velocity_vs_dist.fig')

save('g_oridStruct_wTime.mat','g_oridStruct');
toc