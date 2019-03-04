
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
        %%
        for count3 = 1:numel(fieldnames(dataHHE_struct))
            pnamez = fieldnames(dataHHE_struct);
            dataHHE = dataHHE_struct.(pnamez{count3});
            dataHHN = dataHHN_struct.(pnamez{count3});
            dataHHZ = dataHHZ_struct.(pnamez{count3});
            
            h = figure();
            set(h, 'Position', [1000 1000 1250 1250])
            subplot(2,2,1)
            plot3(dataHHE,dataHHN,dataHHZ)
            hold on
            scatter3(dataHHE,dataHHN,dataHHZ,'o')
            c = 1:numel(dataHHE);
            d = num2cell(c);
            %e = cellstr(d);
            dx = 10; dy = 10;dz=10;
            text(dataHHE+dx,dataHHN+dy,dataHHZ+dz,d);
            grid on
            box on
            xlabel '- W, + E'
            ylabel '- S, + N'
            zlabel '- down, + up'
            title(sprintf('Earthquake %d %s %s',orids(count),stationz,pnamez{count3}))

            subplot(2,2,2)
            plot(dataHHE, dataHHN) %Change so later in time is different color
            hold on
            scatter(dataHHE, dataHHN, 'o')
            text(dataHHE+dx,dataHHN+dy,d);
            grid on
            box on
            xlabel '- W, + E'
            ylabel '- S, + N'

            subplot(2,2,3)
            plot(dataHHN, dataHHZ) %Change so later in time is different color
            hold on
            scatter(dataHHN, dataHHZ, 'o')
            text(dataHHN+dy,dataHHZ+dz,d);
            grid on
            box on
            xlabel '- S, + N'
            ylabel '- down, + up'


            subplot(2,2,4)
            plot(dataHHE, dataHHZ) %Change so later in time is different color
            hold on
            scatter(dataHHE, dataHHZ, 'o')
            text(dataHHE+dx,dataHHZ+dz,d);
            grid on
            box on
            xlabel '- W, + E'
            ylabel '- down, + up'
            hold off
            clear dataHHZ dataHHN dataHHE dx dy
            
            %% Saving file
            hold off
            directory = '/home/a/akfarrell/Uturuncu/Phase/examples/directions';
            filename = sprintf('direction_%d_%s_%s.png',orids(count),stationz,pnamez{count3});
            filename_wPath = fullfile(directory,filename);
            hgexport(h, filename_wPath, hgexport('factorystyle'), 'Format', 'png');
            
            %% Checking direction of waveform
            
        end
    end
    %catch
    %end
    clearvars -except stasz count orids oridStruct allorids fil num_samps ch phases conv_phase
end
toc