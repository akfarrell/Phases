load('3D_1.75_Utur_0.1km.mat');
vals=zeros(numel(modlats),numel(modlats));
inds=zeros(numel(modlats),numel(modlats));
figure()
for i=1:numel(modlats)
    for k=1:numel(modlats)
        [val,ind]=max(pmodel(i,k,40:80));
        ind=ind+40;
        vals(i,k)=val;
        inds(i,k)=ind;
        plot(squeeze(pmodel(i,k,40:80)),-d.zmod(40:80))
        hold on
    end
end
vals
inds