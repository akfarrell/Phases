load('3D_1.75_Utur_0.1km.mat');
vals=zeros(numel(modlats),numel(modlats));
inds=zeros(numel(modlats),numel(modlats));
for i=1:numel(modlats)
    for k=1:numel(modlats)
        [val,ind]=min(pmodel(i,k,100:300));
        ind=ind+100;
        vals(i,k)=val;
        inds(i,k)=ind;
    end
end
vals
inds