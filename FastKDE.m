function dens =  FastKDE(data,x,varargin)
%%
if isempty(varargin)
sig = 1.06*(numel(data))^(1/5);
else
    sig = varargin{1};
end

h = hist(data,x);
f = -ceil(4*sig):ceil(4*sig); f = 1./(sig*sqrt(2*pi))*exp(-0.5*(f/sig).^2); f= f./sum(f);
dens = conv(f,h); dens= dens(ceil(4*sig)+1:end-ceil(4*sig));
dens = dens./sum(dens);

