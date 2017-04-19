function [Cropped, ind] = CropImage(I,c,h,w,varargin)
y1 = max(round(c(1)-h/2),1);
y2 = min(round(c(1)+h/2),size(I,1));
x1 = max(round(c(2)-w/2),1);
x2 = min(round(c(2)+w/2),size(I,2));

Cropped = I(y1:y2,x1:x2);
if nargout>1
    [X,Y] = meshgrid(x1:x2,y1:y2);
    ind = sub2ind(size(I),Y(:),X(:));
end
end