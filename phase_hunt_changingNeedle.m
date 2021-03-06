% function high_corrs = phase_hunt_changingNeedle(allorids,oridStruct,P,fil,cutoff_val)
%orid2166 is origin 203
tic
%stylie = 'abz'; %%%%% CHANGE!!!---------0----
ch = {'HHE','HHN','HHZ'};
cutoff_val = 0.8;
overall_cutoff = 1.8;
P_pad = 12;
S_pad = 25;
SECSPERDAY = 60 * 60 * 24;
num_samps = 31;
close all
for count=203%1:numel(allorids)
    directory = '/home/a/akfarrell/Uturuncu/Phase/wf_objs';
    filename = sprintf('wf_%d.mat',allorids(count));
    filename_wPath = fullfile(directory,filename);
    if exist(filename_wPath,'file')
        load(filename_wPath)
    else
        %create and clean waveform object
        [w_raw,OrS,stations_inEq] = get_wf(allorids(count),oridStruct);
        w_clean = waveform_clean(w_raw, filterobject('b', fil, 2));
        save(filename_wPath,'w_clean', 'OrS', 'stations_inEq');
    end
    stationz = get(w_clean,'station');
    
    %% plot 
    p=figure; hold on;
    set(p, 'Position', [1000 1000 1000 1000])
    latlim = [-22.75 -21.75]; %[southern_limit northern_limit] 
    lonlim = [-67.75 -66.5]; %[western_limit eastern_limit]
    worldmap(latlim, lonlim);
    borders = shaperead('../BOL_adm0.shp', 'UseGeoCoords', true);
    colormap(flipud(colormap))
    arg_borders = shaperead('../ARG_adm0.shp', 'UseGeoCoords', true);
    try
    geoshow(borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
    end
    try
    geoshow(arg_borders, 'DefaultEdgeColor', 'black', 'DefaultFaceColor', 'white');
    end
    scatterm(oridStruct.eq_2166.lat,oridStruct.eq_2166.lon,'*','r')
    scatterm(-22.27, -67.18, 100,'^','filled','r')
    %% some more plotting and such
    for count2 = 1:3:numel(w_clean) %43-45 is PLMN 43:3:43
        if find(intersect(find(strcmp(oridStruct.(OrS).sta,stationz{count2})),find(strcmp(oridStruct.(OrS).phase,'P'))))==1
            siteInd = find(strcmp(siteSub.sta, stationz{count2}));
            site_subset_phase = {'PL03','PLLL','PLMD','PLSP'};
            site_subset_2phase = {'PLDK','PLJR','PLMK','PLMN'};
            site_subset_maybe = {'PLSM'};
            time_vals = {'', '0.51s or 1.23s', '0.47s or 0.79s', '0.45s or 0.72s'};
            if find(strcmp(site_subset_phase,stationz{count2}))
                scatterm(siteSub.lat(siteInd),siteSub.lon(siteInd),'o','r')
            elseif find(strcmp(site_subset_2phase,stationz{count2}))
                scatterm(siteSub.lat(siteInd),siteSub.lon(siteInd),'o','filled','r')
                dink = find(strcmp(site_subset_2phase,stationz{count2}));
                textm(siteSub.lat(siteInd)+0.015,siteSub.lon(siteInd)+0.01,time_vals{dink},'color','r')
            elseif find(strcmp(site_subset_maybe,stationz{count2}))
                scatterm(siteSub.lat(siteInd),siteSub.lon(siteInd),'o','filled','b')
            else
                scatterm(siteSub.lat(siteInd),siteSub.lon(siteInd),'o','filled','k')
            end
            distElev.(stationz{count2}) = sqrt((siteSub.distKM(siteInd))^2 + (siteSub.elev(siteInd)+oridStruct.eq_2166.depth(1))^2);
            distVal = sprintf('%s\n%2.2f\n%2.2f',siteSub.sta{siteInd},siteSub.distKM(siteInd),distElev.(stationz{count2}))
            if strcmp(stationz{count2},'PLCM')
                textm(siteSub.lat(siteInd)+0.045,siteSub.lon(siteInd)-0.025,distVal)
            else
                textm(siteSub.lat(siteInd)-0.04,siteSub.lon(siteInd)-0.025,distVal)
            end
                
%             %HHE = count, HHN = count+1, HHZ = count+2
%             %find P arrival to rule it out
%             ind_P = intersect(find(strcmp(oridStruct.(OrS).sta,stationz{count2})),find(strcmp(oridStruct.(OrS).phase,'P')));
%             time_Parr = oridStruct.(OrS).time_phase(ind_P);
%             try
%                 ind_S = intersect(find(strcmp(oridStruct.(OrS).sta,stationz{count2})),find(strcmp(oridStruct.(OrS).phase,'S')));
%                 time_Sarr = oridStruct.(OrS).time_phase(ind_S);
%             catch
%                 ind_S = numel(get(w_clean(count2), 'data'));
%             end
%             dnum = zeros(1,numel(get(w_clean(count2),'data')));
%             dnum(1) = datenum(get(w_clean(count2),'start'));
%             freq = get(w_clean(count2), 'freq');
%             for l = 2:numel(get(w_clean(count2),'data'))
%                 dnum(l) = datenum((l/freq)/SECSPERDAY+dnum(1));
%             end
%             try
%                 P_ind = find(dnum>=time_Parr, 1,'first');
%             catch
%                 P_ind = 0;
%             end
%             try
%                 S_ind = find(dnum>=time_Sarr, 1,'first');
%             catch
%                 S_ind = numel(get(w_clean(count2), 'data'))-num_samps;
%             end
%             for count3 = 1:numel(ch)
%                 Haystack_data.(ch{count3}) = get(w_clean(count2+count3-1),'data');
%                 needle.(ch{count3}) = Haystack_data.(ch{count3})(P_ind:P_ind+num_samps-1);
%                 %fix indexing issues
%                 lngX = length(Haystack_data.(ch{count3}));
%                 lngY = length(needle.(ch{count3}));
%                 assert(lngX >= lngY);
%                 lags = 0:(lngX-lngY);
%                 for valz = lags
%                     c.(ch{count3})(valz+1) = xcorr(Haystack_data.(ch{count3})(valz+1:valz+lngY) -...
%                         mean(Haystack_data.(ch{count3})(valz+1:valz+lngY)), needle.(ch{count3}) - mean(needle.(ch{count3})),0,'coeff');
%                 end
%                 c_backup.(ch{count3}) = c.(ch{count3}); %USE THIS VALUE TO ANNOTATE THE 3-COMPONENT PLOTS!!!!!!!!!!!-------------
%                 try
%                     if strcmp(stylie,'abz')
%                         c.(ch{count3}) = abs(c.(ch{count3}));
%                     end
%                 catch
%                 end
%                 [m,i.(ch{count3})]=max(c.(ch{count3})(P_ind+P_pad:S_ind-S_pad)); %Pad P_ind and S_ind in range for min
%                 if strcmp(stationz{count2},'PLJR') && count3 == 3
%                     [heck,darn] = min(c.HHZ(700:800));
%                 end
%                 i.(ch{count3}) = i.(ch{count3})+P_ind+P_pad;
%                 %-------------------------UNCOMMENT PLOT IF WANT TO PLOT - 2 lines
%                 %plot_xcorrs(lags, ch{count3}, Haystack_data.(ch{count3}), needle.(ch{count3}), c.(ch{count3}), stationz{count2}, ...
%                     %cutoff_val, i.(ch{count3}), allorids(count), 1,P_ind+P_pad,S_ind-S_pad) %Make sure ind padding is same as line 64!!!
%                 if m>= cutoff_val %%%only half setup for multiple returns that are greater than cutoff_val - need to finish if looking at this
%                     values = find(c.(ch{count3})>=cutoff_val);
%                     %---------------------UNCOMMENT PLOT IF WANT TO PLOT - 2 lines
%                     %plot_xcorrs(lags, ch{count3}, Haystack_data.(ch{count3}), needle.(ch{count3}), c.(ch{count3}), stationz{count2}, ...
%                         %cutoff_val, i.(ch{count3}), allorids(count), 2, P_ind,S_ind)
%                     stations_w_highcorrs.(sprintf('eq_%d',allorids(count))).(stationz{count2}).(ch{count3}).val = m;
%                     stations_w_highcorrs.(sprintf('eq_%d',allorids(count))).(stationz{count2}).(ch{count3}).index = i.(ch{count3});
%     %                 for indie = 1:numel(c.(ch{count3})>=cutoff_val)
%     %                     stations_w_highcorrs.(sprintf('eq_%d',allorids(count))).(stationz{count2}).(ch{count3}).val(indie) = m;
%     %                     stations_w_highcorrs.(sprintf('eq_%d',allorids(count))).(stationz{count2}).(ch{count3}).index(indie) = i.(ch{count3});
%     %                 end
%                 end
%                 clear m
%             end
%            % plot_xcorrs(lags, ch, Haystack_data, needle, c, stationz{count2}, cutoff_val, i, allorids(count), 4,P_ind+P_pad,S_ind-S_pad)
%             clear i
%             clear lags
%             %%
%             c_total = c.(ch{1})+c.(ch{2})+c.(ch{3});
%             inds = find(c_total) >= overall_cutoff;
%             inds_backup = inds;
%             %remove inds that correspond to P-wave pick
%             P_array = P_ind-P_pad:P_ind+P_pad;
%             for indexy = 1:numel(P_array)
%                 if intersect(inds,P_array(indexy))
%                     inds = inds(inds~=intersect(inds,P_array(indexy)));
%                 end
%             end
%             count4 = 1;
%             %remove inds that are too close to each other
%     %         while count4 <= numel(inds)
%     %             count4;
%     %             if numel(intersect(inds,inds(count4)-P_pad:inds(count4)+P_pad)) > 1
%     %                 int_vals = intersect(inds,inds(count4)-P_pad:inds(count4)+P_pad);
%     %                 [dumb,throw_val] = min(c_total(intersect(inds,int_vals)));
%     %                 inds = inds(inds~=int_vals(throw_val));
%     %                 count4 = count4+1;
%     %              end
%     %             count4 = count4+1;
%     %         end
% 
%             for counter = 1:numel(inds)
%                 data_sig(counter*numel(needle.HHZ)-numel(needle.HHZ)+1:counter*numel(needle.HHZ)) = inds(counter):inds(counter)+numel(needle.HHZ)-1;
%             end
%     %%
%     %         try
%         %         plot_xcorrs(lags, ch, Haystack_data, needle, c_total, stationz{count2}, overall_cutoff, data_sig, allorids(count), 3, P_ind,S_ind)
%         %         plot_xcorrs(lags, ch, Haystack_data, needle, c_total, stationz{count2}, overall_cutoff, data_sig, allorids(count), 3, P_ind,S_ind,inds,'best')
%         %         plot_xcorrs(lags, ch, Haystack_data, needle, c_total, stationz{count2}, overall_cutoff, data_sig, allorids(count), 3, P_ind,S_ind,inds_backup,'all')
%     % 
%     % 
%     %         catch
%     %         end
        end
        clear data_sig
        clear inds
        %clear c
        clear ind_P
        clear ind_S
        clear time_Parr
        clear time_S
        clear P_ind
        clear S_ind
        clear needle
        clear Haystack_data
    end
end
clear stylie
toc

%end