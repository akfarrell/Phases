durr=load('velmodels2/utu_vels.txt');
deps = [3];
%deps=deps-1
phases={'P','S'};
clear times
for s=1:numel(siteStruct.sta)
    sta_int=siteStruct.sta{s};
    stInd = find(strcmp(siteStruct.sta, sta_int));
    filez = fullfile('velmodels2',sprintf('velmodel_%s.txt',sta_int));
    fid=fopen(filez,'w');
    durr(1)=siteStruct.elev(stInd)+0.5;
    fprintf(fid,'%2.4f %2.4f %2.4f\n',[durr(:,1)';durr(:,2)';durr(:,3)']);
    fclose(fid);
end