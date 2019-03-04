% rng( 0 ,'twister')  %Set seed to 0
% %Set your starting hyperparameters
% Modality = 3;
% Sigmas = [ 2 1 3]; %Standard Deviations
% Means =  [-5 1 9];
% MixCoefs = [.4 .25 .35]; 
% 
% Length = 1000;
% %Get matrix with individual distributions (non-effective way just for illustr.)
% X = randn(Length,Modality).*repmat(Sigmas,Length,1)+repmat(Means,Length,1);
% %Mix the distributions to your MixCoefs to get the final mixture.
% Y = [ X(1: (Length* MixCoefs(1)),1); 
%       X(1: (Length* MixCoefs(2)),2); 
%       X(1: (Length* MixCoefs(3)),3) ]; 
% 
% %Uncomment to visually inspect your empirical pdf and check multimodality
% %ksdensity(Y)
% 
% %Fit a trimodal to your data
% fitted_obj3 = gmdistribution.fit(Y, 3);
% %Fit a bimodal to your data
% fitted_obj2 = gmdistribution.fit(Y, 2);
% 
% %Check which if the bimodal fit is better based on AIC
% fitted_obj2.AIC < fitted_obj3.AIC
% %Check which if the bimodal fit is better based on BIC
% fitted_obj2.BIC < fitted_obj3.BIC
% %Unsurprisingly the trimodal is better is both cases.
% 
% best_fitted_obj = fitted_obj3; %Watch it the ordering probably it is not 
% % the one you started with
% 
% Est_Sigmas = sqrt( best_fitted_obj.Sigma);  %[ 1.03  3.00  2.00]
% Est_Means =  best_fitted_obj.mu;            %[ 1.09  9.15 -5.11]
% 
% %You can get the  mixing proportions directly by using:
% best_fitted_obj.PComponents
% 
% %Or get an estimate for each reading about the possibility of being part in
% %a specific distribution by using the .posterior() functionality :
% 
% %Posterior "ownership" probabilities of a "-11" reading 
% best_fitted_obj.posterior(-11)
% %ans =
% %    0.0000    0.0000    1.0000 %So pretty certainly it is on the third
% %    component (unsurprisingly) so from the estimated N( -5.11 , 2.00^2)
% 
% %Posterior "ownership" probabilities of a "3" reading 
% best_fitted_obj.posterior(3)
% %ans =
% %    0.7557    0.2434    0.0009 %Probably on the first component ie. from
% %    the N( 1.09, 1.03^2) but actually the second distro is not too
% %    unlikely either

%%
figure()
scatter(1:numel(vel_diamond),sort(vel_diamond))

vel_diamond2 = vel_diamond';
d_3_vel_dist = gmdistribution.fit(vel_diamond2,3);
d_2_vel_dist = gmdistribution.fit(vel_diamond2,2);
d_4_vel_dist = gmdistribution.fit(vel_diamond2,4);
d_2_vel_dist.AIC < d_3_vel_dist.AIC
d_2_vel_dist.BIC < d_3_vel_dist.BIC