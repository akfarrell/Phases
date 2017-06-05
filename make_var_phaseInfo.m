clear all
addpath('/raid/home/a/akfarrell/')
addpath('/raid/home/a/akfarrell/Uturuncu')
addpath('/home/a/akfarrell/Uturuncu/Phase/ind_eq_code')
[oridStruct, allorids] = get_eq_info();

load('good_orids.mat')
for count = 1:14%numel(good_orids)
    directory = '/home/a/akfarrell/Uturuncu/Phase/wf_objs';
    filename = sprintf('wf_%d.mat',good_orids(count));
    filename_wPath = fullfile(directory,filename);
    load(filename_wPath)
    stationz = get(w_clean,'station');
    fname = sprintf('eq_%d',good_orids(count));
    eq_info = str2func(fname);
    P_S = eq_info();
    clear eq_info
    or_ind = find(allorids,good_orids(count));
    %% Put all metrics into g_oridStruct, it's the final structure!
    g_oridStruct.(sprintf('eq_%d',good_orids(count))) = oridStruct.(sprintf('eq_%d',good_orids(count)));
    for count2 = 1:numel(stationz)
        %% load data from xcorr code
        fname2 = sprintf('corr_%d_%s_min.mat',good_orids(count),stationz{count2});
        directory2 = sprintf('/home/a/akfarrell/Uturuncu/Phase/corrs/%d',good_orids(count));
        filename_wPath2 = fullfile(directory2,fname2);
        if ~exist(filename_wPath2)
        else
            load(filename_wPath2)

            %% Find where in g_oridStruct the station is
    %         indz = find(strcmp(g_oridStruct.(sprintf('eq_%d',good_orids(count))).sta,stationz{count2}));
    %         if numel(indz)>1 %if there's more than one instance of the station (i.e. P and S picks)
    %             varz = find(strcmp(g_oridStruct.(sprintf('eq_%d',good_orids(count))).phase(indz),'P'));
    %             indz = indz(varz);
    %             clear varz
    %         end
            if find(strcmp(stationz{count2},fieldnames(P_S)))
                g_oridStruct.(sprintf('eq_%d',good_orids(count))).refl.(stationz{count2}) = P_S.(stationz{count2});
                g_oridStruct.(sprintf('eq_%d',good_orids(count))).refl.(stationz{count2}).P_ind = P_ind;
                g_oridStruct.(sprintf('eq_%d',good_orids(count))).refl.(stationz{count2}).S_ind = S_ind;
                g_oridStruct.(sprintf('eq_%d',good_orids(count))).refl.(stationz{count2}).lags = lags;
                g_oridStruct.(sprintf('eq_%d',good_orids(count))).refl.(stationz{count2}).P_pad = P_pad;
                g_oridStruct.(sprintf('eq_%d',good_orids(count))).refl.(stationz{count2}).S_pad = S_pad;
                g_oridStruct.(sprintf('eq_%d',good_orids(count))).refl.(stationz{count2}).c = c;
    %         clear indz
            end
        end
    end
end
save('g_oridStruct.mat','g_oridStruct')