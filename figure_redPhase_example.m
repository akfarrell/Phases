function figure_redPhase_example(lags, ch, Haystack, needle, c, stationz, cutoff_val, indz, orid, style,P_ind,S_ind,setting,specific)
if style == 4
h = figure();
set(h, 'Position', [1000 1000 1400 1200])  
    for count = 1:numel(ch)
        
        subplot(3,1,count)
        plot(1:numel(Haystack.(ch{count})),Haystack.(ch{count}))
        title(sprintf('PLMK %s Event %d',ch{count},orid))
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
    directory = '/home/a/akfarrell/Uturuncu/Phase/examples';
    if ~exist(directory,'dir')
        mkdir(directory)
    end
    filename = 'example_phase_PLMK_Orid2001';
    filename_wPath = fullfile(directory,filename);
    hgexport(h, filename_wPath, hgexport('factorystyle'), 'Format', 'png');
end