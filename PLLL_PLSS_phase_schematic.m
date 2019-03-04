%% Load topographic data and close previous figure
clear all;
addpath('examples/lineintersect')
load('topo_data.mat');
%load('inf_point_for2166.mat')
[lon, lat] = meshgrid(topo_data.lon, topo_data.lat);
topo_data.z = double(topo_data.z)./1000;
load('conv_phase.mat');
stationz = {'PLLL','PLSS'};
figz = 'yes';
close all;
depths = conv_phase.depths.*-1;

PLLL_velmodel = load('velmodels2/velmodel_PLLL.txt'); %can use for both PLLL and PLSS
PLLL_velmodel(1,1) = conv_phase.PLLLelev-PLLL_velmodel(1,1);
PLLL_velmodel(2,1) = PLLL_velmodel(1,1)-PLLL_velmodel(2,1);

%%
for count = 1:numel(conv_phase.EQNum)
   
     %% APMB Boundary
     for count2 = 1:numel(stationz)
         %% Interp some stuff
        lat_inds = linspace(conv_phase.lat(count),conv_phase.(sprintf('%slat',stationz{count2})),20); %20 evenly spaced samples between
        lon_inds = linspace(conv_phase.lon(count),conv_phase.(sprintf('%slon',stationz{count2})),20); %20 evenly spaced samples between

        elev_vals = zeros(1,numel(lon_inds)); %define sizes to keep box green :)
        distance_vals = elev_vals;
        azimuth_vals = elev_vals;

        if numel(lat_inds) == numel(lon_inds)
            for i = 1:numel(lat_inds)
                if strcmp(figz,'yes')
                    elev_vals(i) = interp2(lon,lat,topo_data.z,lon_inds(i),lat_inds(i));
                end
                [distance_vals(i),azimuth_vals(i)] = distance(conv_phase.lat(count),conv_phase.lon(count),lat_inds(i),lon_inds(i));
            end
            distance_vals = distance_vals.*111.12;
        end

        if strcmp(figz,'yes')
            h = figure();
            set(h, 'Position', [1000 1000 800 800])
            title({conv_phase.orid(count)})
            %% Plot eq and station
            hold on
            scatter(0,depths(count),200,'*','r')
            text(-0.7,depths(count)+0.75,sprintf('%d\n%2.2f',conv_phase.orid(count),depths(count)))
            scatter(distance_vals(end),conv_phase.(sprintf('%selev',stationz{count2})),200,'*','k')
            text(distance_vals(end)-0.7,conv_phase.(sprintf('%selev',stationz{count2}))+0.75,sprintf('%s\n%1.2f',stationz{count2},conv_phase.(sprintf('%selev',stationz{count2}))))
            line([0 distance_vals(end)],[depths(count) conv_phase.(sprintf('%selev',stationz{count2}))])
            
            %% find intersection of lines
            l1 = [0 depths(count) distance_vals(end) conv_phase.(sprintf('%selev',stationz{count2}))];
            l2 = [-2 PLLL_velmodel(1,1) distance_vals(end)+2 PLLL_velmodel(1,1)];
            [x_intersect, y_intersect] = lineintersect(l1,l2);
            hypotenus = sqrt((conv_phase.depths(count)-0.5)^2 + x_intersect^2);
            Udist = conv_phase.distance_w_topo.(stationz{count2})(count)-hypotenus;
            text(-1,-1.5,sprintf('X: %2.2f\nY: %2.2f\nHyp: %2.2f\nUdist: %2.2f',x_intersect, y_intersect,hypotenus,Udist))

            %% Plot interped stuff
            plot(distance_vals,elev_vals,'k')
            xlim([-2 distance_vals(end)+2]) 
            ylim([-5.25 6])
            
            %% Determine timing and distance relations
            time_shift = conv_phase.P.(stationz{count2})(count)-conv_phase.time_ot_to_P.(stationz{count2})(count)*100;
            PUtime = Udist/4.2;
            SUtime = Udist/PLLL_velmodel(1,3);
            PLtime = hypotenus/4.2;
            SLtime = hypotenus/PLLL_velmodel(2,3);
            total_s_time = SLtime+SUtime;
            total_p_time = PLtime+PUtime;
            VP_over_VS= total_s_time/total_p_time;
            diff_P = conv_phase.time_ot_to_P.(stationz{count2})(count)-total_p_time;
            diff_S = conv_phase.Stime_OT.(stationz{count2})(count)-total_s_time;
            text(distance_vals(end)-1,2,sprintf('PUtime: %2.2f\nSUtime: %2.2f\n PLtime: %2.2f\nSLtime: %2.2f\nDiff P: %2.2f\nDiff S: %2.2f\nSvel: %2.2f\nPvel: %2.2f\nVp/Vs %2.2f',...
                PUtime,SUtime, PLtime,SLtime,diff_P,diff_S,conv_phase.distance_w_topo.(stationz{count2})(count)/total_s_time,conv_phase.distance_w_topo.(stationz{count2})(count)/total_p_time,VP_over_VS))
            
            %% Use these timing relations to find thickness/velocity of layer
            thickness = 1; %thickness of layer, in km
            l1 = [0 depths(count) distance_vals(end) conv_phase.(sprintf('%selev',stationz{count2}))];
            l2_2 = [-2 conv_phase.(sprintf('%selev',stationz{count2}))-thickness distance_vals(end)+2 conv_phase.(sprintf('%selev',stationz{count2}))-thickness];
            [x2_intersect, y2_intersect] = lineintersect(l1,l2_2);
            hypotenus2 = sqrt(thickness^2 + (distance_vals(end)-x2_intersect)^2); %distance in layer
            hypotenusmid = Udist-hypotenus2; %distance traveling in between -0.5 km and 1 km below station
            P2time = hypotenusmid/PLLL_velmodel(1,2);
            S2time = hypotenusmid/PLLL_velmodel(1,3);
            
            %% Now looking at different paths, adding combos of P2time, S2 time, PLtime, and SLtime with velocities through 1 km thick layer
            PP = PLtime+P2time;
            SS = SLtime+S2time;
            p1_layer_time = conv_phase.P1time_OT.(stationz{count2})(count) - PP;
            p2_layer_time = conv_phase.P2time_OT.(stationz{count2})(count) - SS;
            S_Layer_Vel = hypotenus2/p1_layer_time;
            P_Layer_Vel = hypotenus2/p2_layer_time;
            
            %% Plot layer velocities
            text(distance_vals(end)-4,PLLL_velmodel(1,1)-0.5,sprintf('Pvel in layer is %1.2f and Svel is %1.2f',P_Layer_Vel,S_Layer_Vel))
            conv_phase.thickness.thickness = thickness;
            conv_phase.thickness.SVel.(stationz{count2})(count) = S_Layer_Vel;
            conv_phase.thickness.PVel.(stationz{count2})(count) = P_Layer_Vel;
            
            %% Plot velmodel stuff
            line([-2 distance_vals(end)+2], [PLLL_velmodel(1,1) PLLL_velmodel(1,1)])
            line([-2 distance_vals(end)+2], [PLLL_velmodel(2,1) PLLL_velmodel(2,1)])
            text(-1.5,PLLL_velmodel(1,1)+0.5,sprintf('%1.1f Pvel is %1.1f and Svel is %1.2f',PLLL_velmodel(1,1),PLLL_velmodel(1,2),PLLL_velmodel(1,3)))
            text(-1.5,PLLL_velmodel(2,1)+0.25,sprintf('%1.1f Pvel is %1.1f and Svel is %1.2f',PLLL_velmodel(2,1),PLLL_velmodel(2,2),PLLL_velmodel(2,3)))
            %text(-1.5,PLLL_velmodel(2,1)-0.5,sprintf('Pvel is %1.1f and Svel is %1.2f',PLLL_velmodel(3,2),PLLL_velmodel(3,3)))
            
            %% Plot phase velocities
            VP_over_VS_actual = conv_phase.vel_P.(stationz{count2})(count)/conv_phase.vel_S.(stationz{count2})(count);
            conv_phase.Vp_Vs_measured.(stationz{count2})(count) = VP_over_VS_actual;
            text(-1,3, sprintf('Vels:\nP: %1.2f\nAP1: %1.2f\nP1: %1.2f\nP2: %1.2f\nS: %1.2f\nVp/Vs %1.2f\n\nTime: %2.2f\nDist: %2.2f',conv_phase.vel_P.(stationz{count2})(count),...
                conv_phase.vel_AltPhase1.(stationz{count2})(count),conv_phase.vel_Phase1.(stationz{count2})(count),conv_phase.vel_Phase2.(stationz{count2})(count),...
                conv_phase.vel_S.(stationz{count2})(count),VP_over_VS_actual,conv_phase.time_ot_to_P.(stationz{count2})(count),conv_phase.distance_w_topo.(stationz{count2})(count)))
        end  
        %% Save figure
        hold off
        directory = sprintf('/home/a/akfarrell/Uturuncu/Phase/examples/xsections');
        if ~exist(directory,'dir')
            mkdir(directory)
        end
        filename = sprintf('%s_%d.png',stationz{count2},conv_phase.orid(count));
        filename_wPath = fullfile(directory,filename);
        hgexport(h, filename_wPath, hgexport('factorystyle'), 'Format', 'png');
        clearvars -except PLLL_velmodel lon lat topo_data conv_phase stationz figz depths count
     end
end
save('conv_phase.mat','conv_phase')
save('examples/conv_phase.mat','conv_phase')