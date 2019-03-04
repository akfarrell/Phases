function timez = automated_arrivals(siteSub,eq_depth)
tic
% durr=load('velmodel_toEdit.txt');
%deps = [46,56,66,76];
deps = [3];
%deps=deps-1
phases={'P','S'};
clear timez
for s=1:numel(siteSub.sta)
    sta_int=siteSub.sta{s};
    stInd = find(strcmp(siteSub.sta, sta_int));
%     fid=fopen('velmodel_toEdit2.txt','w');
%     durr(1)=siteSub.elev(stInd)+0.1;
%     fprintf(fid,'%2.4f %2.4f %2.4f\n',[durr(:,1)';durr(:,2)';durr(:,3)']);
%     fclose(fid);  %these steps done in makin_files.m
    for d=1:numel(deps)
        for p=1:numel(phases)
            result = perl('tt_reflected.pl', fullfile('velmodels2',sprintf('velmodel_%s.txt',sta_int)), num2str(siteSub.distKM(stInd)), num2str(siteSub.elev(stInd)+eq_depth), phases{p},num2str(deps(d)),'0');
            timez.(siteSub.sta{s}).refl.(phases{p})(d) = str2double(result((strfind(result,'TravelTime[s]   ')+16):(strfind(result,' -- Takeoff')-1)));
            if isnan(timez.(siteSub.sta{s}).refl.(phases{p})(d))
                timez.(siteSub.sta{s}).refl.(phases{p})(d) = str2double(result((strfind(result,'TravelTime[s]  ')+15):(strfind(result,' -- Takeoff')-1)));
            end
            if d==numel(deps)
                result2 = perl('tt_direct.pl', fullfile('velmodels2',sprintf('velmodel_%s.txt',sta_int)), num2str(siteSub.distKM(stInd)), num2str(siteSub.elev(stInd)+eq_depth), phases{p},'0');
                timez.(siteSub.sta{s}).d.(phases{p}) = str2double(result2((strfind(result2,'TravelTime[s]   ')+16):(strfind(result2,' -- Takeoff')-1)));
                if isnan(timez.(siteSub.sta{s}).d.(phases{p}))
                    timez.(siteSub.sta{s}).d.(phases{p}) = str2double(result2((strfind(result2,'TravelTime[s]  ')+15):(strfind(result2,' -- Takeoff')-1)));
                end
            end
        end
    end
end
timez.order=deps;
toc
timez2=timez;
if isequal(timez2,timez) == 0
    clear timez
    error('There are NaN values in timez - clearing timez')
end
clear timez2
end