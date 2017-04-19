function Phi = SDF(BW)
if all(BW(:))
    Phi = zeros(size(BW));
    return;
end
Cont = abs(BW - imerode(BW,ones(3)));
Dist = bwdist(Cont);
Phi = Dist.*BW -Dist.*(~BW);
%Phi = bwdist(~BW)-bwdist(BW)+(~BW);