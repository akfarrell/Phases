%% This is scavenged code, not functional

tic
addpath('/raid/home/a/akfarrell/')
addpath('/raid/home/a/akfarrell/Uturuncu')
clear all;
clc;
load('conv_phase.mat')
fil=[2 25];
[oridStruct, allorids] = get_eq_info(); 

ch = {'HHE','HHN','HHZ'};
num_samps = 15;
%close all
orids = conv_phase.orid;
close all
stasz = {'PLLL','PLSS'};
phases = {'P','AltPhase1','Phase1','Phase2','S'};

%% loop
for count = 1:numel(orids)
    %try
    directory = '/home/a/akfarrell/Uturuncu/Phase/wf_objs';
    filename = sprintf('wf_%d.mat',orids(count));
    filename_wPath = fullfile(directory,filename);
    if exist(filename_wPath,'file')
        load(filename_wPath)
    else
        create and clean waveform object
        [w_raw,OrS,stations_inEq] = get_wf(orids(count),oridStruct);
        w_clean = waveform_clean(w_raw, filterobject('b', fil, 2));
        save(filename_wPath,'w_clean', 'OrS', 'stations_inEq');
    end
    stnz = get(w_clean,'station');
    for count2 = 1:numel(stasz)
        stationz = stasz{count2};
        indz = find(strcmp(stnz,stationz));
        dataHHE_raw = get(w_clean(indz(1)),'data');
        dataHHN_raw = get(w_clean(indz(2)),'data');
        dataHHZ_raw = get(w_clean(indz(3)),'data');
        
        for count3 = 1:numel(phases)
            indexval = conv_phase.(phases{count3}).(stationz);
            if indexval == 999
                clear indexval
            else
                dataHHE_struct.(phases{count3}) = dataHHE_raw(indexval:indexval+num_samps-1);
                dataHHN_struct.(phases{count3}) = dataHHN_raw(indexval:indexval+num_samps-1);
                dataHHZ_struct.(phases{count3}) = dataHHZ_raw(indexval:indexval+num_samps-1);
                clear indexval
            end
        end
        clear count3
        pnamez = fieldnames(dataHHE_struct);
        %%
        for count3 = 1:numel(pnamez)
            dataHHE = dataHHE_struct.(pnamez{count3});
            dataHHN = dataHHN_struct.(pnamez{count3});
            dataHHZ = dataHHZ_struct.(pnamez{count3});
%% From old directionality code
            N = numel(dataHHE);
            ZZ = (1/N)*sum(dataHHZ.*dataHHZ);
            ZN = (1/N)*sum(dataHHZ.*dataHHN);
            ZE = (1/N)*sum(dataHHZ.*dataHHE);

            NN = (1/N)*sum(dataHHN.*dataHHN);
            NE = (1/N)*sum(dataHHN.*dataHHE);

            EE = (1/N)*sum(dataHHE.*dataHHE);
            
            correlation_matrix = [ ZZ ZN ZE ;...
                                                 ZN NN NE;...
                                                 ZE NE EE];

            [seig_vec, seig_mat] = eig(correlation_matrix);
            su1 = seig_vec(:,3);
            su2 = seig_vec(:,2);
            su3 = seig_vec(:,1);
            p_inc.(sprintf('eq_%d',orids(count))).(stasz{count2}).(pnamez{count3}) = acosd(abs(su1(1)));
            p_az.(sprintf('eq_%d',orids(count))).(stasz{count2}).(pnamez{count3}) = atand(su1(2)*sign(su1(1)))/(su1(3)*sign(su1(1)));
            %[m, I] = nanmax(data);
            %data_subsetE = dataE()
            ZZ2.(pnamez{count3})(count) = ZZ;
            ZN2.(pnamez{count3})(count) = ZN;
            ZE2.(pnamez{count3})(count) = ZE;

            NN2.(pnamez{count3})(count) = NN;
            NE2.(pnamez{count3})(count) = NE;

            EE2.(pnamez{count3})(count) = EE;
            clear Z? N? EE dataHH? correlation_matrix seig_vec seig_mat su?
        end
    end
end
clear count3
%% Sum it all up
for count3 = 1:numel(pnamez)
    num_sta = numel(ZZ2.(pnamez{count3}));
    ZZ_sum.(pnamez{count3}) = sum(ZZ2.(pnamez{count3}))/num_sta;
    ZN_sum.(pnamez{count3}) = sum(ZN2.(pnamez{count3}))/num_sta;
    ZE_sum.(pnamez{count3}) = sum(ZE2.(pnamez{count3}))/num_sta;

    NN_sum.(pnamez{count3}) = sum(NN2.(pnamez{count3}))/num_sta;
    NE_sum.(pnamez{count3}) = sum(NE2.(pnamez{count3}))/num_sta;

    EE_sum.(pnamez{count3}) = sum(EE2.(pnamez{count3}))/num_sta;

    sum_corr_matrix.(pnamez{count3}) = [ ZZ_sum.(pnamez{count3}) ZN_sum.(pnamez{count3}) ZE_sum.(pnamez{count3});...
                        ZN_sum.(pnamez{count3}) NN_sum.(pnamez{count3}) NE_sum.(pnamez{count3});...
                        ZE_sum.(pnamez{count3}) NE_sum.(pnamez{count3}) EE_sum.(pnamez{count3})];
    [eig_vec, eig_mat] = eig(sum_corr_matrix.(pnamez{count3}));
    eig_vec.(pnamez{count3}) = eig_vec;
    eig_mat.(pnamez{count3}) = eig_mat;

%%
    % l1 = eig_mat(3,3);
    % l2 = eig_mat(2,2);
    % l3 = eig_mat(1,1);

    u1.(pnamez{count3}) = eig_vec.(pnamez{count3})(:,3);
    u2.(pnamez{count3}) = eig_vec.(pnamez{count3})(:,2);
    u3.(pnamez{count3}) = eig_vec.(pnamez{count3})(:,1);

    phase_az.(pnamez{count3}) = atand((u1.(pnamez{count3})(2)*sign(u1.(pnamez{count3})(1)))/(u1.(pnamez{count3})(3)*sign(u1.(pnamez{count3})(1))));
    phase_inc.(pnamez{count3}) = acosd(abs(u1.(pnamez{count3})(1)));
    clear eig_vec eig_mat
end
%end
