function varargout = CalcBGLighting(I,varargin)
if length(varargin)>=1
    M = varargin{1};
    if isempty(M)
         M = true(size(I));
    end
else
    M = true(size(I));
end

h = fspecial('gaussian',500,100);
mu = median(I(M(:)));
I2= I;
I2(~M) = mu;
B = imfilter(I2,h,'symmetric');
ICorrected = I-B;


%ICorrected(M) = I(M)-B(M);
%ICorrected(M) = I(M)-B(M);

varargout{1} = ICorrected;
varargout{2} = B;


