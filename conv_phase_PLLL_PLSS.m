%% Loading in inflection point data
tic
close all
clear all
load('conv_phase.mat');
load('g_oridStruct.mat');
erqs = conv_phase.orid; %make sure that this is the same as the variable 'count' or problems
w_topo = 'y'; %Change to be with or without topography

%% Plotting dist vs time
p=figure; hold on;
set(p, 'Position', [1000 1000 1200 1200])
for count=1:numel(conv_phase.orid)
    stationz = fieldnames(conv_phase.P);
    ot = conv_phase.OT(count);
    inds_phase_p = find(strcmp(g_oridStruct.(sprintf('eq_%d',erqs(count))).phase, 'P'));
    stime = min(g_oridStruct.(sprintf('eq_%d',erqs(count))).time_phase(inds_phase_p))-datenum(0,0,0,0,0,3);
    eq_lat = conv_phase.lat(count);
    eq_lon = conv_phase.lon(count);
    for count2 = 1:numel(stationz)
        %% caclulate the time indices, in sec*100
        P_ind = conv_phase.P.(stationz{count2})(count);
        time_Parr = datenum(0,0,0,0,0,P_ind/100)+stime;%map p arrival onto frame of reference of OT, add to start time
        time_ot_to_P = ceil(etime(datevec(time_Parr),datevec(ot))*100);
        
        %% Find site information, with distance from eq to site, in km
        siteLat = conv_phase.(sprintf('%slat',stationz{count2}));
        siteLon = conv_phase.(sprintf('%slon',stationz{count2}));
        dist = (distance(siteLat,siteLon, eq_lat, eq_lon)*111.12); %in km
        
        %% Calculate distance with depth and elevation in mind
        siteElev = conv_phase.(sprintf('%selev',stationz{count2})); %in km
        eqDep = conv_phase.depths(count); %in km
        total_elev = siteElev+eqDep;
        actual_dist = sqrt(total_elev^2+dist^2);
        if strcmp(w_topo, 'y')
            dist_backup = dist;
            dist = actual_dist;
        end
        
        %% Actually plotting stuff
        if conv_phase.P.(stationz{count2})(count) ~= 999
            scatter(dist,time_ot_to_P/100,'+','k') %plot P wave arrival
            conv_phase.vel_P.(stationz{count2})(count) = dist./(time_ot_to_P/100);
        else
            conv_phase.vel_P.(stationz{count2})(count) = NaN;
        end
        conv_phase.time_ot_to_P.(stationz{count2})(count) = time_ot_to_P/100;
        conv_phase.distance_w_topo.(stationz{count2})(count) = actual_dist;
        if conv_phase.Phase1.(stationz{count2})(count) ~= 999
            t_plot1 = conv_phase.Phase1.(stationz{count2})(count)-P_ind+time_ot_to_P;
            conv_phase.time_Phase1.(stationz{count2})(count) = t_plot1;
            scatter(dist,t_plot1/100,'d','f','r')
            conv_phase.vel_Phase1.(stationz{count2})(count) = dist./(t_plot1/100);
        else
            conv_phase.time_Phase1.(stationz{count2})(count) = 999;
            conv_phase.vel_Phase1.(stationz{count2})(count) = NaN; %nanmean
        end
        
        if conv_phase.Phase2.(stationz{count2})(count) ~= 999
            t_plot2 = conv_phase.Phase2.(stationz{count2})(count)-P_ind+time_ot_to_P;
            conv_phase.time_Phase2.(stationz{count2})(count) = t_plot2;
            scatter(dist,t_plot2/100,'s','f','b')
            conv_phase.vel_Phase2.(stationz{count2})(count) = dist./(t_plot2/100);
        else
            conv_phase.time_Phase2.(stationz{count2})(count) = 999;
            conv_phase.vel_Phase2.(stationz{count2})(count) = NaN;
        end
        
        if conv_phase.AltPhase1.(stationz{count2})(count) ~= 999
            t_plot3 = conv_phase.AltPhase1.(stationz{count2})(count)-P_ind+time_ot_to_P;
            conv_phase.time_AltPhase1.(stationz{count2})(count) = t_plot3;
            scatter(dist,t_plot3/100,'*','g')
            conv_phase.vel_AltPhase1.(stationz{count2})(count) = dist./(t_plot3/100);
        else
            conv_phase.time_AltPhase1.(stationz{count2})(count) = 999;
            conv_phase.vel_AltPhase1.(stationz{count2})(count) = NaN;
        end
        
        if conv_phase.S.(stationz{count2})(count) ~= 999
            t_plot4 = conv_phase.S.(stationz{count2})(count)-P_ind+time_ot_to_P;
            conv_phase.time_S.(stationz{count2})(count) = t_plot4;
            scatter(dist,t_plot4/100,'c','f','m')
            conv_phase.vel_S.(stationz{count2})(count) = dist./(t_plot4/100);
        else
            conv_phase.time_S.(stationz{count2})(count) = 999;
            conv_phase.vel_S.(stationz{count2})(count) = NaN;
        end
        
        clear t_plot?
    end
    clear P_ind
    clear ot
end
%xlim([0 30])
xlabel('Distance (km)')
ylabel('Time (s)')

%% Saving file
hold off
directory = '/home/a/akfarrell/Uturuncu/Phase/examples';
filename = 'timeVSdist_PLLL_PLSS.png';
filename_wPath = fullfile(directory,filename);
hgexport(p, filename_wPath, hgexport('factorystyle'), 'Format', 'png');

savefig('examples/velocity_vs_dist_PLLL_PLSS.fig')

save('conv_phase.mat','conv_phase')
toc