clear all; close all; clc;

[NUM,TXT,RAW] = xlsread('Exotic_Phases_code.xlsx');
TXT_cut = TXT(2:end, :);
[orid_sort, Index] = sort(NUM(:,1));
NUM_s = NUM(Index,:); TXT_s = TXT_cut(Index,:); RAW_s = RAW(Index,:);

eq_orid = NUM_s(:,1);
station = TXT_s(:,2);

orid = unique(eq_orid);
orid_str = {1:length(orid)};

for count = 1:length(orid)
    orid_str{count} = ['o', num2str(orid(count))];
end

stas = {};
for count = 1:length(orid)
    for count2 = 1:length(eq_orid)
        if orid(count) == eq_orid(count2)
            stas{count2} = station{count2};
        else
            continue
        end
    end
    events.(orid_str{count}).station = stas(find(~cellfun(@isempty,stas)));
    clear stas
end