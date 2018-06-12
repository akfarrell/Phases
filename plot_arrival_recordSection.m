eq_id = 'eq_2166';
w = prep_for_rs(w_clean(1:3:end),eq_id,siteStruct,oridStruct);
load('phaseStruct_backup/pS1.mat');
stnz = get(w, 'station');
sta = fieldnames(phaseStruct.(eq_id));


%%
sze = 0;
EVDP = -get(w(1),'EVDP'); %convention - below sea level is negative!!!!
EVLA = get(w(1),'EVLA');
EVLO = get(w(1),'EVLO');
for count = 1:numel(sta)
    w_ind = find(strcmp(stnz,sta(count)));
    op_ind = intersect(find(strcmp(oridStruct.(eq_id).sta, sta(count))), find(strcmp(oridStruct.(eq_id).phase,'P')));
    os_ind = intersect(find(strcmp(oridStruct.(eq_id).sta, sta(count))), find(strcmp(oridStruct.(eq_id).phase,'S'))); 
    STLA = get(w(w_ind),'STLA');
    STLO = get(w(w_ind),'STLO');
    STEL = get(w(w_ind),'STEL');
    d = distance(EVLA,EVLO,STLA,STLO)*111.12; %epicentral distance
    height = STEL-EVDP; %diff in elevation btw event and station
    h = sqrt(d^2+height^2);%distance taking elevation into account
    hypee(count) = h; % all with -ee are for putting into phaseStruct
    p_arree(count) = oridStruct.(eq_id).time_phase(op_ind);
    
    hyp(sze+1:sze+numel(phaseStruct.(eq_id).(sta{count}))) = h;
    dist(sze+1:sze+numel(phaseStruct.(eq_id).(sta{count}))) = d;
    time(sze+1:sze+numel(phaseStruct.(eq_id).(sta{count})))...
        = phaseStruct.(eq_id).(sta{count});
    p_arr(sze+1:sze+numel(phaseStruct.(eq_id).(sta{count}))) ...
        = oridStruct.(eq_id).time_phase(op_ind);
    stationsies(sze+1:sze+numel(phaseStruct.(eq_id).(sta{count}))) = sta(count);
    if ~isempty(os_ind)
        s_arr(sze+1:sze+numel(phaseStruct.(eq_id).(sta{count}))) ...
            = oridStruct.(eq_id).time_phase(os_ind);
        s_arree(count) = oridStruct.(eq_id).time_phase(os_ind);
    else
        s_arr(sze+1:sze+numel(phaseStruct.(eq_id).(sta{count}))) = 1;  
        s_arree(count) = 1;
    end
    sze = sze+numel(phaseStruct.(eq_id).(sta{count}));
    clear w_ind STLA STLO STEL height h d
end

time = time*24*3600;
tp = p_arr*24*3600;
ts = s_arr*24*3600;
mint = min(tp)-1;
t_plot = time-mint;
tp_plot = tp-mint;
ts_plot = ts-mint;
c = min(time)-min(tp);
%%
close all
h1 = figure();
set(h1,'Position', [1000 1000 1400 1200])
scatter(hyp,t_plot,'+')
hold on

inds = find(ts_plot>1);
ts_plot = ts_plot(inds);
hyp_s = hyp(inds);

scatter(hyp,tp_plot,'+','m')
scatter(hyp_s,ts_plot,'+','r')
title('Using Station Elevation and Event Depth')
xlabel('Distance (km)')
ylabel('Time (s)')
ylim([0 max(ts_plot)+1])


% PLOT THE FIRST INTERPRETATION OF VELOCITIES FOR PHASES
% line([27.45 43.93],[7.85+c 9.25+c],'Color','k')
% line([27.45 43.93],[6.6+c 8.6+c],'Color','k')
% text(35,7.8,'8.6 km/s')
% text(35,9.3,'12.6 km/s')


% PLOT THE SECOND INTERPRETATION OF VELOCITIES FOR PHASES
% line([11.9 27.7],[2.9 8.4],'Color','k')
% line([16.1 28.2],[4.0 7.2],'Color','k')
% line([15.9 44.6],[3.6 9.0],'Color','k')
% text(27.5,7.9,'2.9 km/s')
% text(27.5,6.6,'3.8 km/s')
% text(40,8.7,'5.3 km/s')

% KEEP THESE LINES AND TEXTS FOR EITHER PLOT
% line([8.06 47.87],[1.53+c 8.57+c],'Color','k')
% text(6.5,1.7,'6.3 km/s')

% hline1 = gline(h1);
% hline2 = gline(h1);
% hline3 = gline(h1);
legend('Phase','P','S','Location','Best')
view([0 -90])
set(gca,'ydir','reverse')

coeffs_p = polyfit(hyp, tp_plot, 1);
% Get fitted values
fittedXp = linspace(min(hyp), max(hyp), 200);
fittedYp = polyval(coeffs_p, fittedXp);
% Plot the fitted line
plot(fittedXp, fittedYp, 'm-', 'LineWidth', 0.5);


coeffs_s = polyfit(hyp_s, ts_plot, 1);
% Get fitted values
fittedXs = linspace(min(hyp_s), max(hyp_s), 200);
fittedYs = polyval(coeffs_s, fittedXs);
% Plot the fitted line
plot(fittedXs, fittedYs, 'r-', 'LineWidth', 0.5);

text(29.9,10,sprintf('%2.2f km/s',1/coeffs_s(1)),'Color','r')
text(30,5,sprintf('%2.2f km/s',1/coeffs_p(1)),'Color','m')

%% Using values derived from the epicentral distance
% h2 = figure();
% set(h2,'Position', [1000 1000 1400 1200])
% scatter(dist,t_plot,'+')
% hold on
% scatter(dist,tp_plot,'+','m')
% scatter(dist,ts_plot,'+','r')
% title('Using Epicentral Distance')
% xlabel('Distance (km)')
% ylabel('Time (s)')
% ylim([0 max(ts_plot)+1])
% line([26.09 42.91],[7.91+c 9.25+c],'Color','k')
% line([25.98 43.15],[6.63+c 8.63+c],'Color','k')
% line([1.21 47.06],[1.33+c 8.578+c],'Color','k')
% % hline4 = gline(h2);
% % hline5 = gline(h2);
% % hline6 = gline(h2);
% legend('Phase','P','S','Location','Best')
% view([0 -90])
% set(gca,'ydir','reverse')