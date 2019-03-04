load('good_evids.mat')
addpath(genpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu/Phase'))
addpath('/Users/alexandrafarrell/Desktop/akfarrell/')
addpath('/Users/alexandrafarrell/Desktop/akfarrell/Uturuncu')
[evidStruct, allorids] = get_eq_info();
ml = zeros(1,numel(good_evids));
for count = 1:numel(good_evids)
    eq = sprintf('eq_%d',good_evids(count));
    ml(count) = evidStruct.(eq).ml(1);
end
%%
inds = find(ml~=-999);
min(ml(inds))
max(ml(inds))
numel(find(ml(inds)<=0))