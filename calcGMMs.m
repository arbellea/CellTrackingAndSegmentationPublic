function [gmmFG, gmmBG] = calcGMMs(Data, Labels, Kfg, Kbg)

gmmOpts = statset('MaxIter',500);
gmmBG = fitgmdist(Data(Labels(:)==0,:),Kbg,'Options',gmmOpts, 'RegularizationValue', 0.01);
%%
gmmFG = fitgmdist(Data(Labels(:)>0,:),Kfg,'Options',gmmOpts,  'RegularizationValue', 0.01);
%%
gmmBG = gmdistribution(gmmBG.mu,gmmBG.Sigma,ones(1,Kbg)./Kbg);
gmmFG = gmdistribution(gmmFG.mu,gmmFG.Sigma,ones(1,Kfg)./Kfg);

