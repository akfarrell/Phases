%find which orid has the most phases picked
y=fields(oridStruct);
for i=1:numel(y)
   p(i)=numel(oridStruct.(y{i}).phase);
end
clear i
[val,ind]=max(p)
%returns val=50 and ind=223, eq 2166
%explore that motherfucker
Or=2166;
[w_raw,OrS,stations_inEq] = get_wf(Or,oridStruct)
%% - same as big_stat_shit.eqs.eq_2166
% clear stat_shit
% statz={'timeres','smajax','sminax','sdepth','stime'};
% for s=1:numel(statz)
%     if s==1
%         stat_shit.(sprintf('max_%s',statz{s}))=max(oridStruct.(OrS).(statz{s}));
%         stat_shit.(sprintf('min_%s',statz{s}))=min(oridStruct.(OrS).(statz{s}));
%         stat_shit.(sprintf('std_%s',statz{s}))=std(oridStruct.(OrS).(statz{s}));
%         stat_shit.(sprintf('mean_%s',statz{s}))=mean(oridStruct.(OrS).(statz{s}));
%         stat_shit.(sprintf('range_%s',statz{s}))=range(oridStruct.(OrS).(statz{s}));
%     else
%         stat_shit.(statz{s})=oridStruct.(OrS).(statz{s})(1);
%     end
% end
% clear s
%%
%compare to all of the shit for all orids
clear big_stat_shit
for numz=1:numel(y)
    for s=1:numel(statz)
        if strcmp(statz{s},'timeres')
            big_stat_shit.eqs.(y{numz}).(sprintf('max_%s',statz{s}))=max(oridStruct.(y{numz}).(statz{s}));
            big_stat_shit.eqs.(y{numz}).(sprintf('min_%s',statz{s}))=min(oridStruct.(y{numz}).(statz{s}));
            big_stat_shit.eqs.(y{numz}).(sprintf('std_%s',statz{s}))=std(oridStruct.(y{numz}).(statz{s}));
            big_stat_shit.eqs.(y{numz}).(sprintf('mean_%s',statz{s}))=mean(oridStruct.(y{numz}).(statz{s}));
            big_stat_shit.eqs.(y{numz}).(sprintf('range_%s',statz{s}))=range(oridStruct.(y{numz}).(statz{s}));
        else
            big_stat_shit.eqs.(y{numz}).(statz{s})=oridStruct.(y{numz}).(statz{s})(1);
            big_stat_shit.eqs.(y{numz}).(statz{s})=oridStruct.(y{numz}).(statz{s})(1);
            big_stat_shit.eqs.(y{numz}).(statz{s})=oridStruct.(y{numz}).(statz{s})(1);
            big_stat_shit.eqs.(y{numz}).(statz{s})=oridStruct.(y{numz}).(statz{s})(1);
            big_stat_shit.eqs.(y{numz}).(statz{s})=oridStruct.(y{numz}).(statz{s})(1);
        end
    end
end
clear numz
clear s
for s=1:numel(statz)
    for numz=1:numel(y)
        h.smajax(numz)=oridStruct.(y{numz}).smajax(1);
        h.sminax(numz)=oridStruct.(y{numz}).sminax(1);
        h.sdepth(numz)=oridStruct.(y{numz}).sdepth(1);
        h.stime(numz)=oridStruct.(y{numz}).stime(1);
        h.timeres(numz)=mean(oridStruct.(y{numz}).timeres);
    end
    big_stat_shit.stats.(statz{s}).(sprintf('max_%s',statz{s}))=max(h.(statz{s}));
    big_stat_shit.stats.(statz{s}).(sprintf('min_%s',statz{s}))=min(h.(statz{s}));
    big_stat_shit.stats.(statz{s}).(sprintf('std_%s',statz{s}))=std(h.(statz{s}));
    big_stat_shit.stats.(statz{s}).(sprintf('mean_%s',statz{s}))=mean(h.(statz{s}));
    big_stat_shit.stats.(statz{s}).(sprintf('range_%s',statz{s}))=range(h.(statz{s}));
end
clear numz s h