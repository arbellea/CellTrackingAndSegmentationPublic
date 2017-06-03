function P_FG = calcGMMProb(I,gmmFG, gmmBG, alpha)

Data = calcFeatures({I});
dens_BG = pdf(gmmBG,Data);
dens_FG = pdf(gmmFG,Data);
S = (1-alpha)*dens_BG+ alpha*dens_FG;
P_FG = alpha*dens_FG./S;
P_FG = reshape(P_FG,size(I));
;

