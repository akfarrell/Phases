function plot_haney_directionality(xmax,val_try, Haystack_data,P_ind,S_ind,...
    azim,incd, ellip, rpol, ppol,az, inc,needle,tvecC_z)
ch = fieldnames(Haystack_data);
pad = 100;

if az<180
    az2 = az;
else
    az2 = az-180;
end
if inc < 90
    inc2 = inc;
else
    inc2 = 180-inc;
end


figure()
subplot(8,1,1)
plot(Haystack_data.HHZ)

hold on
try
    for ct = 1:numel(needle)
        plot(needle(ct):needle(ct)+31,Haystack_data.HHZ(needle(ct):needle(ct)+31),'r')
    end
catch
end
line([P_ind P_ind], [min(Haystack_data.HHZ), max(Haystack_data.HHZ)],...
    'Color','m','LineStyle',':','LineWidth',2)
try
    line([S_ind S_ind], [min(Haystack_data.HHZ), max(Haystack_data.HHZ)],...
        'Color','m','LineStyle',':','LineWidth',2)
catch
end
line([val_try val_try], [min(Haystack_data.HHZ), max(Haystack_data.HHZ)],...
    'Color','b','LineStyle',':','LineWidth',2)
ylim([min(Haystack_data.HHZ), max(Haystack_data.HHZ)])
title('HHZ')
% try
%     xlim([P_ind-pad S_ind+pad])
% catch
%     xlim([P_ind-pad numel(Haystack_data.HHZ)])
% end
if S_ind > xmax
    xlim([P_ind-pad, S_ind+pad])
else
    xlim([P_ind-pad, xmax])
%     xlim([max(P_ind-pad,1), min(S_ind+pad,xmax)])
end

subplot(8,1,2)
xs = linspace(max(P_ind-100,1),S_ind+100,numel(azim))./100;
h = plot(tvecC_z,azim);
hold on
try
    for ct = 1:numel(needle)
        plot(needle(ct)/100:0.01:needle(ct)/100+0.31,azim(needle(ct):needle(ct)+31),'r')
    end
catch
end
if numel(azim)>100
    set(h,'YData',azim)
    hold on
    line([P_ind/100 P_ind/100], [min(azim), max(azim)],...
        'Color','m','LineStyle',':','LineWidth',2)
    line([val_try/100 val_try/100], [min(azim), max(azim)],...
        'Color','b','LineStyle',':','LineWidth',2)
    ylim([min(azim), max(azim)])
    try
        line([S_ind/100 S_ind/100], [min(azim), max(azim)],...
            'Color','m','LineStyle',':','LineWidth',2)
    catch
    end
else
    hold on
    line([P_ind/100 P_ind/100], [min(azim), max(azim)],...
        'Color','m','LineStyle',':','LineWidth',2)
    line([val_try/100 val_try/100], [min(azim), max(azim)],...
        'Color','b','LineStyle',':','LineWidth',2)
    try
        line([S_ind/100 S_ind/100], [min(azim), max(azim)],...
            'Color','m','LineStyle',':','LineWidth',2)
    catch
    end
end
title('Azimuth')
line([0 numel(Haystack_data.HHZ)],[az az],'Color','g','LineStyle',':','LineWidth',1)
line([0 numel(Haystack_data.HHZ)],[az2 az2],'Color','g','LineStyle',':','LineWidth',1)
hold off
% try
%     xlim([(P_ind-pad)/100 (S_ind+pad)/100])
% catch
%     xlim([(P_ind-pad)/100 numel(Haystack_data.HHZ)])
% end
if S_ind > xmax
    xlim([(P_ind-pad)/100, (S_ind+pad)/100])
else
    xlim([(P_ind-pad)/100, xmax/100])
%     xlim([max(P_ind-pad,1), min(S_ind+pad,xmax)])
end

subplot(8,1,3)
plot(Haystack_data.(ch{2}))
hold on
try
    for ct = 1:numel(needle)
        plot(needle(ct):needle(ct)+31,Haystack_data.(ch{2})(needle(ct):needle(ct)+31),'r')
    end
catch
end
line([P_ind P_ind], [min(Haystack_data.(ch{2})), max(Haystack_data.(ch{2}))],...
    'Color','m','LineStyle',':','LineWidth',2)
try
    line([S_ind S_ind], [min(Haystack_data.(ch{2})), max(Haystack_data.(ch{2}))],...
        'Color','m','LineStyle',':','LineWidth',2)
catch
end
line([val_try val_try], [min(Haystack_data.(ch{2})), max(Haystack_data.(ch{2}))],...
    'Color','b','LineStyle',':','LineWidth',2)
ylim([min(Haystack_data.(ch{2})), max(Haystack_data.(ch{2}))])
title(sprintf('%s',ch{2}))
if S_ind > xmax
    xlim([P_ind-pad, S_ind+pad])
else
    xlim([P_ind-pad, xmax])
%     xlim([max(P_ind-pad,1), min(S_ind+pad,xmax)])
end

subplot(8,1,4)          
plot(tvecC_z,incd)
hold on
try
    for ct = 1:numel(needle)
        plot(needle(ct)/100:0.01:needle(ct)/100+0.31,incd(needle(ct):needle(ct)+31),'r')
    end
catch
end
line([P_ind/100 P_ind/100], [min(incd), max(incd)],...
    'Color','m','LineStyle',':','LineWidth',2)
line([val_try/100 val_try/100], [min(incd), max(incd)],...
    'Color','b','LineStyle',':','LineWidth',2)
line([0 numel(Haystack_data.HHZ)],[inc inc],'Color','g','LineStyle',':','LineWidth',1)
line([0 numel(Haystack_data.HHZ)],[inc2 inc2],'Color','g','LineStyle',':','LineWidth',1)
try
    line([S_ind/100 S_ind/100], [min(incd), max(incd)],...
        'Color','m','LineStyle',':','LineWidth',2)
catch
end
hold off
if S_ind > xmax
    xlim([(P_ind-pad)/100, (S_ind+pad)/100])
else
    xlim([(P_ind-pad)/100, xmax/100])
%     xlim([max(P_ind-pad,1), min(S_ind+pad,xmax)])
end
title('Incidence angle')

subplot(8,1,5)
plot(Haystack_data.(ch{3}))
hold on
try
    for ct = 1:numel(needle)
        plot(needle(ct):needle(ct)+31,Haystack_data.(ch{3})(needle(ct):needle(ct)+31),'r')
    end
catch
end
line([P_ind P_ind], [min(Haystack_data.(ch{3})), max(Haystack_data.(ch{3}))],...
    'Color','m','LineStyle',':','LineWidth',2)
try
    line([S_ind S_ind], [min(Haystack_data.(ch{3})), max(Haystack_data.(ch{3}))],...
        'Color','m','LineStyle',':','LineWidth',2)
catch
end
line([val_try val_try], [min(Haystack_data.(ch{3})), max(Haystack_data.(ch{3}))],...
    'Color','b','LineStyle',':','LineWidth',2)
ylim([min(Haystack_data.(ch{3})), max(Haystack_data.(ch{3}))])
title(sprintf('%s',ch{3}))
if S_ind > xmax
    xlim([P_ind-pad, S_ind+pad])
else
    xlim([P_ind-pad, xmax])
%     xlim([max(P_ind-pad,1), min(S_ind+pad,xmax)])
end
            
            
subplot(8,1,6)          
plot(tvecC_z,ellip)
hold on
try
    for ct = 1:numel(needle)
        plot(needle(ct)/100:0.01:needle(ct)/100+0.31,ellip(needle(ct):needle(ct)+31),'r')
    end
catch
end
line([P_ind/100 P_ind/100], [min(ellip), max(ellip)],...
    'Color','m','LineStyle',':','LineWidth',2)
line([val_try/100 val_try/100], [min(ellip), max(ellip)],...
    'Color','b','LineStyle',':','LineWidth',2)
line([0 numel(Haystack_data)],[1 1],'Color','g','LineStyle',':','LineWidth',1)
try
    line([S_ind/100 S_ind/100], [min(ellip), max(ellip)],...
        'Color','m','LineStyle',':','LineWidth',2)
catch
end
hold off
if S_ind > xmax
    xlim([(P_ind-pad)/100, (S_ind+pad)/100])
else
    xlim([(P_ind-pad)/100, xmax/100])
%     xlim([max(P_ind-pad,1), min(S_ind+pad,xmax)])
end
title('Ellipticity')           
           
           
subplot(8,1,7)          
plot(tvecC_z,rpol)
hold on
try
    for ct = 1:numel(needle)
        plot(needle(ct)/100:0.01:needle(ct)/100+0.31,rpol(needle(ct):needle(ct)+31),'r')
    end
catch
end
line([P_ind/100 P_ind/100], [min(rpol), max(rpol)],...
    'Color','m','LineStyle',':','LineWidth',2)
line([val_try/100 val_try/100], [min(rpol), max(rpol)],...
    'Color','b','LineStyle',':','LineWidth',2)
try
    line([S_ind/100 S_ind/100], [min(rpol), max(rpol)],...
        'Color','m','LineStyle',':','LineWidth',2)
catch
end
hold off
if S_ind > xmax
    xlim([(P_ind-pad)/100, (S_ind+pad)/100])
else
    xlim([(P_ind-pad)/100, xmax/100])
%     xlim([max(P_ind-pad,1), min(S_ind+pad,xmax)])
end
title('Rectilinearity')

subplot(8,1,8)    
plot(tvecC_z,ppol)
hold on
title('Planarity')
try
    for ct = 1:numel(needle)
        plot(needle(ct)/100:0.01:needle(ct)/100+0.31,ppol(needle(ct):needle(ct)+31),'r')
    end
catch
end
line([P_ind/100 P_ind/100], [min(ppol), max(ppol)],...
    'Color','m','LineStyle',':','LineWidth',2)
line([val_try/100 val_try/100], [min(ppol), max(ppol)],...
    'Color','b','LineStyle',':','LineWidth',2)
try
    line([S_ind/100 S_ind/100], [min(ppol), max(ppol)],...
        'Color','m','LineStyle',':','LineWidth',2)
catch
end
hold off
if S_ind > xmax
    xlim([(P_ind-pad)/100, (S_ind+pad)/100])
else
    xlim([(P_ind-pad)/100, xmax/100])
%     xlim([max(P_ind-pad,1), min(S_ind+pad,xmax)])
end
% 
% % Now create destination graph
% fn = findobj('type','figure');
% n = length(fn);
% figure(n+1)
% ax = zeros(8,1);
% for i = 1:8
%     ax(i)=subplot(8,1,i);
% end
% % Now copy contents of each figure over to destination figure
% % Modify position of each axes as it is transferred
% for i = 1:8
%     figure(n-8+i)
%     h = get(gcf,'Children');
%     %newh = copyobj(h,9);
%     newh = copyobj(h,n+1);
%     for j = 1:length(newh)
%         posnewh = get(newh(j),'Position');
%         possub  = get(ax(i),'Position');
%         set(newh(j),'Position',...
%         [posnewh(1) possub(2) posnewh(3) possub(4)])
%     end
%     delete(ax(i));
% end
% figure(n+1)
% 
% close(n-8+1:n)
% fn2 = findobj('type','figure');
% n2 = length(fn2);
% figure(n2)
% copyobj(allchild(n+1),n2)
% close(n+1)

