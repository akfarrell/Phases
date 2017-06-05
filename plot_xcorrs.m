function plot_xcorrs(lags, ch, Haystack, needle, c, stationz, cutoff_val, indz, orid, style,P_ind,S_ind,setting,specific)
%Style 1 vs Style 2 changes where it's saved
if nargin < 13
    setting = 'nullz';
    specific = 0;
end
h = figure();
set(h, 'Position', [1000 1000 1400 1200])
if style ~=4
subplot(4,1,1)
plot(lags,c)
ylim([min(c), max(c)])
hold on
line([lags(1) lags(end)],[cutoff_val cutoff_val],'Color','r')
line([P_ind, P_ind],[min(c), max(c)],'Color','m','LineStyle',':','LineWidth',1)
line([S_ind, S_ind],[min(c), max(c)],'Color','m','LineStyle',':','LineWidth',1)

if style ==1 || style == 2
    xlim([lags(1) numel(Haystack)])
    title(sprintf('xcorr summary %s %d %s samples %1.1f',ch,numel(needle),stationz,cutoff_val))
    
    subplot(4,1,2)
    plot(1:numel(Haystack),Haystack)
    title(sprintf('Haystack %s',ch))
    xlim([0 numel(Haystack)])
    hold on
    plot(indz:indz+numel(needle)-1, Haystack(indz:indz+numel(needle)-1),'r')
    line([P_ind, P_ind],[min(Haystack), max(Haystack)],'Color','m','LineStyle',':','LineWidth',1)
    line([S_ind, S_ind],[min(Haystack), max(Haystack)],'Color','m','LineStyle',':','LineWidth',1)
    ylim([min(Haystack), max(Haystack)])

    subplot(4,1,3)
    plot(1:numel(needle),needle)
    title(sprintf('Needle %s',ch))
    xlim([1 numel(needle)])
    ylim([min(needle), max(needle)])

    subplot(4,1,4)
    plot(1:numel(needle),needle-Haystack(indz:indz+numel(needle)-1))
    xlim([1 numel(needle)])
    title('Diff Needle and Haystack Maximum Correlation')
    filename = sprintf('xcorr_%d_%s_%dsamples_%1.1f_%s.png',orid,ch,numel(needle),cutoff_val,stationz);
    if style == 1
        directory = '/home/a/akfarrell/Uturuncu/Phase/xcorr_figs';
        filename_wPath = fullfile(directory,filename);
        hgexport(h, filename_wPath, hgexport('factorystyle'), 'Format', 'png');
        if ~exist(directory,'dir')
            mkdir(directory)
        end
    elseif style ==2
        cut_str = num2str(cutoff_val);
        cut_str = strrep(cut_str,'.','pt');
        directory = sprintf('/home/a/akfarrell/Uturuncu/Phase/xcorr_figs/%d/%s/%s',orid,stationz,cut_str);
        if ~exist(directory,'dir')
            mkdir(directory)
        end
        filename_wPath = fullfile(directory,filename);
        hgexport(h, filename_wPath, hgexport('factorystyle'), 'Format', 'png');
    end
    
    
elseif style == 3
    if ~ischar(setting)
        if ~ischar(setting)
            indstr = strread(num2str(setting),'%s');
        end
        for count2 = 1:numel(setting)
            if count2 < numel(setting)
                if setting(count2)+5 >=setting(count2+1)
                    text(setting(count2)+10,c(setting(count2))+0.25,sprintf('%s\n%1.2f',indstr{count2},c(setting(count2))))
                else
                    text(setting(count2)+10,c(setting(count2))-0.25,sprintf('%s\n%1.2f',indstr{count2},c(setting(count2))))
                end
            elseif count2 == numel(setting)
                text(setting(count2)+10,c(setting(count2))-0.25,sprintf('%s\n%1.2f',indstr{count2},c(setting(count2))))
            end
            hold on
        end
    end
    xlim([lags(1) numel(lags)+numel(needle.HHZ)-1])
    title(sprintf('xcorr summary %d %s samples %1.1f',numel(needle.HHZ),stationz,cutoff_val))
    
    for count = 1:numel(ch)
        subplot(4,1,count+1)
        plot(1:numel(Haystack.(ch{count})),Haystack.(ch{count}))
        title(sprintf('Haystack %s',ch{count}))
        xlim([0 numel(Haystack.(ch{count}))])
        hold on
        for counter2 = 1:numel(indz)/numel(needle.HHZ)
            plot(indz(counter2*numel(needle.HHZ)-numel(needle.HHZ)+1:counter2*numel(needle.HHZ)),...
                Haystack.(ch{count})(indz(counter2*numel(needle.HHZ)-numel(needle.HHZ)+1:counter2*numel(needle.HHZ))),'r')
            hold on
        end
        ylim([min(Haystack.(ch{count})), max(Haystack.(ch{count}))])
        hold off
    end

    directory = sprintf('/home/a/akfarrell/Uturuncu/Phase/xcorr_figs/total/%d',orid);
    if ~exist(directory,'dir')
        mkdir(directory)
    end
    if ~ischar(setting)
        filename = sprintf('annotated_%s_xcorr_%d_%dsamples_%1.1f_%s.png',specific,orid,numel(needle.HHZ),cutoff_val,stationz);
    else
        filename = sprintf('xcorr_%d_%dsamples_%1.1f_%s.png',orid,numel(needle.HHZ),cutoff_val,stationz);
    end
    filename_wPath = fullfile(directory,filename);
    hgexport(h, filename_wPath, hgexport('factorystyle'), 'Format', 'png');
end
elseif style == 4
   
    for count = 1:numel(ch)
        subplot(6,1,count*2-1)
        plot(lags,c.(ch{count}),'k')
        if count == 1
            title(sprintf('3-comp xcorr %d %s samples %1.1f',numel(needle.HHZ),stationz,cutoff_val))
        end
        ylim([min(c.(ch{count})), max(c.(ch{count}))])
        hold on
        line([lags(1) lags(end)],[cutoff_val cutoff_val],'Color','r')
        line([P_ind, P_ind],[min(c.(ch{count})), max(c.(ch{count}))],'Color','m','LineStyle',':','LineWidth',1)
        line([S_ind, S_ind],[min(c.(ch{count})), max(c.(ch{count}))],'Color','m','LineStyle',':','LineWidth',1)
        %xlim([lags(1) numel(lags)+numel(needle.HHZ)-1])
        if S_ind + 300 > numel(Haystack.(ch{count}))
            xlim([P_ind-300, S_ind])
        else
            xlim([P_ind-300, S_ind+300])
        end
        
        subplot(6,1,count*2)
        plot(1:numel(Haystack.(ch{count})),Haystack.(ch{count}))
        title(sprintf('Haystack %s %d %d',ch{count},orid,indz.(ch{count})))
        xlim([0 numel(Haystack.(ch{count}))])
        hold on
        plot(indz.(ch{count}):indz.(ch{count})+numel(needle.HHZ)-1, Haystack.(ch{count})(indz.(ch{count}):indz.(ch{count})+numel(needle.HHZ)-1),'r')
        line([P_ind, P_ind],[min(Haystack.(ch{count})), max(Haystack.(ch{count}))],'Color','m','LineStyle',':','LineWidth',1)
        line([S_ind, S_ind],[min(Haystack.(ch{count})), max(Haystack.(ch{count}))],'Color','m','LineStyle',':','LineWidth',1)
        ylim([min(Haystack.(ch{count})), max(Haystack.(ch{count}))])
        if S_ind + 300 > numel(Haystack.(ch{count}))
            xlim([P_ind-300, S_ind])
        else
            xlim([P_ind-300, S_ind+300])
        end
    end
    hold off
    directory = sprintf('/home/a/akfarrell/Uturuncu/Phase/xcorr_figs/total/%d',orid);
    if ~exist(directory,'dir')
        mkdir(directory)
    end
    %filename = sprintf('3comp_xcorr_%d_%dsamples_%1.1f_%s.png',orid,numel(needle.HHZ),cutoff_val,stationz);
    %filename_wPath = fullfile(directory,filename);
    %hgexport(h, filename_wPath, hgexport('factorystyle'), 'Format', 'png');
end
end