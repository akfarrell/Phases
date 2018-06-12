function plot_threecomp_directionality(xmax,val_try, Haystack_data,P_ind,S_ind,...
    tcrp_az,tcrp_inc, tcrp_rec, tcrp_plan, tcrp_en,az, inc,varargin)
ch = fieldnames(Haystack_data);
pad = 100;
figure
plot(Haystack_data.HHZ)

hold on
try
    s = varargin{2};
    for ct = 1:numel(s)
        plot(s(ct):s(ct)+31,Haystack_data.HHZ(s(ct):s(ct)+31),'r')
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
                xlim([P_ind-pad, xmax])
            else
                xlim([P_ind-pad, S_ind+pad])
            end


h = plot(tcrp_az);
hold on
try
    for ct = 1:numel(s)
        tcrp_az_data = varargin{1};
        plot(s(ct)/100:0.01:s(ct)/100+0.31,tcrp_az_data(s(ct):s(ct)+31),'r')
    end
catch
end
if numel(varargin{1})>100
    set(h,'YData',varargin{1})
    hold on
    line([P_ind/100 P_ind/100], [min(varargin{1}), max(varargin{1})],...
        'Color','m','LineStyle',':','LineWidth',2)
    line([val_try/100 val_try/100], [min(varargin{1}), max(varargin{1})],...
        'Color','b','LineStyle',':','LineWidth',2)
    ylim([min(varargin{1}), max(varargin{1})])
    try
        line([S_ind/100 S_ind/100], [min(varargin{1}), max(varargin{1})],...
            'Color','m','LineStyle',':','LineWidth',2)
    catch
    end
else
    hold on
    line([P_ind/100 P_ind/100], [min(get(tcrp_az,'data')), max(get(tcrp_az,'data'))],...
        'Color','m','LineStyle',':','LineWidth',2)
    line([val_try/100 val_try/100], [min(get(tcrp_az,'data')), max(get(tcrp_az,'data'))],...
        'Color','b','LineStyle',':','LineWidth',2)
    try
        line([S_ind/100 S_ind/100], [min(get(tcrp_az,'data')), max(get(tcrp_az,'data'))],...
            'Color','m','LineStyle',':','LineWidth',2)
    catch
    end
end
line([0 numel(Haystack_data.HHZ)],[az az],'Color','g','LineStyle',':','LineWidth',1)
hold off
% try
%     xlim([(P_ind-pad)/100 (S_ind+pad)/100])
% catch
%     xlim([(P_ind-pad)/100 numel(Haystack_data.HHZ)])
% end
           if S_ind > xmax
                xlim([(P_ind-pad)/100, xmax/100])
            else
                xlim([(P_ind-pad)/100, (S_ind+pad)/100])
            end

figure()
plot(Haystack_data.(ch{2}))
hold on
try
    for ct = 1:numel(s)
        plot(s(ct):s(ct)+31,Haystack_data.(ch{2})(s(ct):s(ct)+31),'r')
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
                xlim([P_ind-pad, xmax])
            else
                xlim([P_ind-pad, S_ind+pad])
            end
plot(tcrp_inc)
hold on
try
    for ct = 1:numel(s)
        tcrp_inc_data = get(tcrp_inc,'data');
        plot(s(ct)/100:0.01:s(ct)/100+0.31,tcrp_inc_data(s(ct):s(ct)+31),'r')
    end
catch
end
line([P_ind/100 P_ind/100], [min(get(tcrp_inc,'data')), max(get(tcrp_inc,'data'))],...
    'Color','m','LineStyle',':','LineWidth',2)
line([val_try/100 val_try/100], [min(get(tcrp_inc,'data')), max(get(tcrp_inc,'data'))],...
    'Color','b','LineStyle',':','LineWidth',2)
line([0 numel(Haystack_data.HHZ)],[inc inc],'Color','g','LineStyle',':','LineWidth',1)
try
    line([S_ind/100 S_ind/100], [min(get(tcrp_inc,'data')), max(get(tcrp_inc,'data'))],...
        'Color','m','LineStyle',':','LineWidth',2)
catch
end
hold off
           if S_ind > xmax
                xlim([(P_ind-pad)/100, xmax/100])
            else
                xlim([(P_ind-pad)/100, (S_ind+pad)/100])
            end


figure()
plot(Haystack_data.(ch{3}))
hold on
try
    for ct = 1:numel(s)
        plot(s(ct):s(ct)+31,Haystack_data.(ch{3})(s(ct):s(ct)+31),'r')
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
                xlim([P_ind-pad, xmax])
            else
                xlim([P_ind-pad, S_ind+pad])
            end

plot(tcrp_rec)
hold on
try
    for ct = 1:numel(s)
        tcrp_rec_data = get(tcrp_rec,'data');
        plot(s(ct)/100:0.01:s(ct)/100+0.31,tcrp_rec_data(s(ct):s(ct)+31),'r')
    end
catch
end
line([P_ind/100 P_ind/100], [min(get(tcrp_rec,'data')), max(get(tcrp_rec,'data'))],...
    'Color','m','LineStyle',':','LineWidth',2)
line([val_try/100 val_try/100], [min(get(tcrp_rec,'data')), max(get(tcrp_rec,'data'))],...
    'Color','b','LineStyle',':','LineWidth',2)
line([0 numel(Haystack_data)],[1 1],'Color','g','LineStyle',':','LineWidth',1)
try
    line([S_ind/100 S_ind/100], [min(get(tcrp_rec,'data')), max(get(tcrp_rec,'data'))],...
        'Color','m','LineStyle',':','LineWidth',2)
catch
end
hold off
           if S_ind > xmax
                xlim([(P_ind-pad)/100, xmax/100])
            else
                xlim([(P_ind-pad)/100, (S_ind+pad)/100])
            end


plot(tcrp_plan)
hold on
try
    for ct = 1:numel(s)
        tcrp_plan_data = get(tcrp_plan,'data');
        plot(s(ct)/100:0.01:s(ct)/100+0.31,tcrp_plan_data(s(ct):s(ct)+31),'r')
    end
catch
end
line([P_ind/100 P_ind/100], [min(get(tcrp_plan,'data')), max(get(tcrp_plan,'data'))],...
    'Color','m','LineStyle',':','LineWidth',2)
line([val_try/100 val_try/100], [min(get(tcrp_plan,'data')), max(get(tcrp_plan,'data'))],...
    'Color','b','LineStyle',':','LineWidth',2)
try
    line([S_ind/100 S_ind/100], [min(get(tcrp_plan,'data')), max(get(tcrp_plan,'data'))],...
        'Color','m','LineStyle',':','LineWidth',2)
catch
end
hold off
           if S_ind > xmax
                xlim([(P_ind-pad)/100, xmax/100])
            else
                xlim([(P_ind-pad)/100, (S_ind+pad)/100])
            end


plot(tcrp_en)
hold on
try
    for ct = 1:numel(s)
        tcrp_en_data = get(tcrp_en,'data');
        plot(s(ct)/100:0.01:s(ct)/100+0.31,tcrp_en_data(s(ct):s(ct)+31),'r')
    end
catch
end
line([P_ind/100 P_ind/100], [min(get(tcrp_en,'data')), max(get(tcrp_en,'data'))],...
    'Color','m','LineStyle',':','LineWidth',2)
line([val_try/100 val_try/100], [min(get(tcrp_en,'data')), max(get(tcrp_en,'data'))],...
    'Color','b','LineStyle',':','LineWidth',2)
try
    line([S_ind/100 S_ind/100], [min(get(tcrp_en,'data')), max(get(tcrp_en,'data'))],...
        'Color','m','LineStyle',':','LineWidth',2)
catch
end
hold off
           if S_ind > xmax
                xlim([(P_ind-pad)/100, xmax/100])
            else
                xlim([(P_ind-pad)/100, (S_ind+pad)/100])
            end

% Now create destination graph
fn = findobj('type','figure');
n = length(fn);
figure(n+1)
ax = zeros(8,1);
for i = 1:8
    ax(i)=subplot(8,1,i);
end
% Now copy contents of each figure over to destination figure
% Modify position of each axes as it is transferred
for i = 1:8
    figure(n-8+i)
    h = get(gcf,'Children');
    %newh = copyobj(h,9);
    newh = copyobj(h,n+1);
    for j = 1:length(newh)
        posnewh = get(newh(j),'Position');
        possub  = get(ax(i),'Position');
        set(newh(j),'Position',...
        [posnewh(1) possub(2) posnewh(3) possub(4)])
    end
    delete(ax(i));
end
figure(n+1)

close(n-8+1:n)
fn2 = findobj('type','figure');
n2 = length(fn2);
figure(n2)
copyobj(allchild(n+1),n2)
close(n+1)

end