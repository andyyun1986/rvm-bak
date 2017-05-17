function [BASIS, Scales] = SB2_PreProcessBasis(BASIS)
[~,M]	= size(BASIS);
Scales	= sqrt(sum(BASIS.^2));
Scales(Scales==0)	= 1;
for m=1:M
  BASIS(:,m)	= BASIS(:,m) / Scales(m);
end
