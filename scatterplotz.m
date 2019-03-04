clear all;
load('conv_phase.mat');
load('g_oridStruct.mat');
close all
figure()
scatter(conv_phase.time_ot_to_P.PLLL,conv_phase.Stime_OT.PLLL-conv_phase.time_ot_to_P.PLLL)
hold on;
scatter(conv_phase.time_ot_to_P.PLSS,conv_phase.Stime_OT.PLSS-conv_phase.time_ot_to_P.PLSS,'r')

%%
orids = conv_phase.orid;
for count = 1:numel(orids)
    eqs = sprintf('eq_%d',orids(count));
    time_phase = g_oridStruct.(eqs).time_phase;
    phase = g_oridStruct.(eqs).phase;
    sta = g_oridStruct.(eqs).sta;
    unique_sta = unique(sta);
    otime = g_oridStruct.(eqs).time_origin(1);
    staz = {0};
    for count2 = 1:numel(unique_sta)
        if numel(find(strcmp(unique_sta(count2),sta))) == 2
            staz{numel(staz)+1} = sta{find(strcmp(unique_sta(count2),sta))};
        end
    end
    staz=staz(2:end);
    for count3 = 1:numel(staz)
        inds = find(strcmp(staz{count3},sta));
        P_val(count3) = time_phase(inds(find(strcmp('P',phase(inds)))));
        S_val(count3) = time_phase(inds(find(strcmp('S',phase(inds)))));;
        S_minus_P(count3) = etime(datevec(S_val(count3)),datevec(P_val(count3)));
        OT_to_P(count3) = etime(datevec(P_val(count3)),datevec(otime));
    end
    figure()
    scatter(OT_to_P,S_minus_P)
    hold on;
    for count4 = 1:numel(staz)
        text(OT_to_P(count4)-0.1,S_minus_P(count4),staz{count4})
    end
    myfit = polyfit(OT_to_P,S_minus_P,1);
    x=linspace(min(OT_to_P),max(OT_to_P),10);
    y=myfit(1)*x+myfit(2);
    plot(x,y)
    text(mean(x),mean(y),sprintf('%1.2f',myfit(1)+1))

    Vp_Vs(count) = myfit(1)+1;
    
    title(sprintf('Earthquake %d',orids(count)))
    xlabel('P Arrival (s)')
    ylabel('S - P Time (s)')
    hold off;
    clear staz;clear otime;clear time_phase;clear phase;clear sta;clear inds;clear P_val;clear S_val;clear S_minus_P;clear OT_to_P;clear myfit;clear x;clear y;
end
save('Vp_Vs.mat','Vp_Vs')