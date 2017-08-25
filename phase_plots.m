function phase_plots(Haystack_data,P_ind, S_ind, azim,incd,ellip,orid,stasz)
    load('conv_phase.mat')
    pad = 150;
    k = figure();
    set(k, 'Position', [1000 1000 1400 1200])
    subplot(6,1,1)
        plot(Haystack_data.HHZ)
        hold on
        title('HHZ')
        if S_ind + pad > numel(Haystack_data.HHZ)
            xlim([P_ind-pad, S_ind])
            xlimz = [P_ind-pad, S_ind];
        else
            xlim([P_ind-pad, S_ind+pad])
            xlimz = [P_ind-pad, S_ind+pad];
        end
        ylim([min(Haystack_data.HHZ(xlimz(1):xlimz(2))),max(Haystack_data.HHZ(xlimz(1):xlimz(2)))])
    subplot(6,1,2)
        plot(azim,'k')
        ylabel('azimuth (deg clockwise from north)')
        title('Azimuth')
        if S_ind + pad > numel(azim)
            xlim([P_ind-pad, S_ind])
            xlimz = [P_ind-pad, S_ind];
        else
            xlim([P_ind-pad, S_ind+pad])
            xlimz = [P_ind-pad, S_ind+pad];
        end
        ylim([min(azim(xlimz(1):xlimz(2))),max(azim(xlimz(1):xlimz(2)))])
    subplot(6,1,3)
        plot(Haystack_data.HHE)
        if S_ind + pad > numel(Haystack_data.HHE)
            xlim([P_ind-pad, S_ind])
            xlimz = [P_ind-pad, S_ind];
        else
            xlim([P_ind-pad, S_ind+pad])
            xlimz = [P_ind-pad, S_ind+pad];
        end
        ylim([min(Haystack_data.HHE(xlimz(1):xlimz(2))),max(Haystack_data.HHE(xlimz(1):xlimz(2)))])
        title('HHE')
    subplot(6,1,4)
        plot(incd,'k')
        ylabel('Incidence angle (deg from vertical)')
        title('Incidence')
        if S_ind + pad > numel(incd)
            xlim([P_ind-pad, S_ind])
            xlimz = [P_ind-pad, S_ind];
        else
            xlim([P_ind-pad, S_ind+pad])
            xlimz = [P_ind-pad, S_ind+pad];
        end
        ylim([min(incd(xlimz(1):xlimz(2))),max(incd(xlimz(1):xlimz(2)))])
    subplot(6,1,5)
        plot(Haystack_data.HHN)
        if S_ind + pad > numel(Haystack_data.HHN)
            xlim([P_ind-pad, S_ind])
            xlimz = [P_ind-pad, S_ind];
        else
            xlim([P_ind-pad, S_ind+pad])
            xlimz = [P_ind-pad, S_ind+pad];
        end
        ylim([min(Haystack_data.HHN(xlimz(1):xlimz(2))) max(Haystack_data.HHN(xlimz(1):xlimz(2)))])
        title('HHN')
    subplot(6,1,6)
        plot(ellip,'k')
        ylabel('Ellipticity (intermediate/major')
        title('Ellipticity')
        if S_ind + pad > numel(ellip)
            xlim([P_ind-pad, S_ind])
            xlimz = [P_ind-pad, S_ind];
        else
            xlim([P_ind-pad, S_ind+pad])
            xlimz = [P_ind-pad, S_ind+pad];
        end
        ylim([min(ellip(xlimz(1):xlimz(2))),max(ellip(xlimz(1):xlimz(2)))])
    if any(find(orid==conv_phase.orid)) && any(strcmp(stasz,fieldnames(conv_phase.P)))
        indz = find(orid==conv_phase.orid);
        phases = {'P','S','Phase1','Phase2','AltPhase1'};
        colorz = {'m','m','b','r','g'};
        for count = 1:numel(phases)
            if conv_phase.(phases{count}).(stasz)(indz) == 999
            else
                subplot(6,1,1)
                    line([conv_phase.(phases{count}).(stasz)(indz), conv_phase.(phases{count}).(stasz)(indz)], ...
                        [min(Haystack_data.HHZ(xlimz(1):xlimz(2))), max(Haystack_data.HHZ(xlimz(1):xlimz(2)))], 'Color',colorz{count},'LineStyle',':','LineWidth',2)                
                subplot(6,1,2)
                    line([conv_phase.(phases{count}).(stasz)(indz), conv_phase.(phases{count}).(stasz)(indz)], ...
                        [min(azim(xlimz(1):xlimz(2))), max(azim(xlimz(1):xlimz(2)))], 'Color',colorz{count},'LineStyle',':','LineWidth',2)                
                subplot(6,1,3)
                    line([conv_phase.(phases{count}).(stasz)(indz), conv_phase.(phases{count}).(stasz)(indz)], ...
                        [min(Haystack_data.HHE(xlimz(1):xlimz(2))), max(Haystack_data.HHE(xlimz(1):xlimz(2)))], 'Color',colorz{count},'LineStyle',':','LineWidth',2)                
                subplot(6,1,4)
                    line([conv_phase.(phases{count}).(stasz)(indz), conv_phase.(phases{count}).(stasz)(indz)], ...
                        [min(incd(xlimz(1):xlimz(2))), max(incd(xlimz(1):xlimz(2)))], 'Color',colorz{count},'LineStyle',':','LineWidth',2)                
                subplot(6,1,5)
                    line([conv_phase.(phases{count}).(stasz)(indz), conv_phase.(phases{count}).(stasz)(indz)], ...
                        [min(Haystack_data.HHN(xlimz(1):xlimz(2))), max(Haystack_data.HHN(xlimz(1):xlimz(2)))], 'Color',colorz{count},'LineStyle',':','LineWidth',2)                
                subplot(6,1,6)
                    line([conv_phase.(phases{count}).(stasz)(indz), conv_phase.(phases{count}).(stasz)(indz)], ...
                        [min(ellip(xlimz(1):xlimz(2))), max(ellip(xlimz(1):xlimz(2)))], 'Color',colorz{count},'LineStyle',':','LineWidth',2)                
            end
        end
    else
        subplot(6,1,1)
            line([P_ind P_ind], [min(Haystack_data.HHZ(xlimz(1):xlimz(2))), max(Haystack_data.HHZ(xlimz(1):xlimz(2)))],...
                'Color','m','LineStyle',':','LineWidth',2) 
            line([S_ind S_ind], [min(Haystack_data.HHZ(xlimz(1):xlimz(2))), max(Haystack_data.HHZ(xlimz(1):xlimz(2)))],...
                'Color','m','LineStyle',':','LineWidth',2)
        subplot(6,1,2)
            line([P_ind P_ind], [min(azim(xlimz(1):xlimz(2))), max(azim(xlimz(1):xlimz(2)))],...
                'Color','m','LineStyle',':','LineWidth',2) 
            line([S_ind S_ind], [min(azim(xlimz(1):xlimz(2))), max(azim(xlimz(1):xlimz(2)))],...
                'Color','m','LineStyle',':','LineWidth',2) 
        subplot(6,1,3)
            line([P_ind P_ind], [min(Haystack_data.HHE(xlimz(1):xlimz(2))), max(Haystack_data.HHE(xlimz(1):xlimz(2)))],...
                'Color','m','LineStyle',':','LineWidth',2) 
            line([S_ind S_ind], [min(Haystack_data.HHE(xlimz(1):xlimz(2))), max(Haystack_data.HHE(xlimz(1):xlimz(2)))],...
                'Color','m','LineStyle',':','LineWidth',2)           
        subplot(6,1,4)
            line([P_ind P_ind], [min(incd(xlimz(1):xlimz(2))), max(incd(xlimz(1):xlimz(2)))],...
                'Color','m','LineStyle',':','LineWidth',2) 
            line([S_ind S_ind], [min(incd(xlimz(1):xlimz(2))), max(incd(xlimz(1):xlimz(2)))],...
                'Color','m','LineStyle',':','LineWidth',2)              
        subplot(6,1,5)
            line([P_ind P_ind], [min(Haystack_data.HHN(xlimz(1):xlimz(2))), max(Haystack_data.HHN(xlimz(1):xlimz(2)))],...
                'Color','m','LineStyle',':','LineWidth',2) 
            line([S_ind S_ind], [min(Haystack_data.HHN(xlimz(1):xlimz(2))), max(Haystack_data.HHN(xlimz(1):xlimz(2)))],...
                'Color','m','LineStyle',':','LineWidth',2) 
        subplot(6,1,6)
            line([P_ind P_ind], [min(ellip(xlimz(1):xlimz(2))), max(ellip(xlimz(1):xlimz(2)))],...
                'Color','m','LineStyle',':','LineWidth',2) 
            line([S_ind S_ind], [min(ellip(xlimz(1):xlimz(2))), max(ellip(xlimz(1):xlimz(2)))],...
                'Color','m','LineStyle',':','LineWidth',2)  
    end
    hold off
end