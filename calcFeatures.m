function Data_norm = calcFeatures(Images)

h = fspecial('gaussian', [5,5], 0.7);
sobh = fspecial('sobel');
lOgh = fspecial('log');
laph = fspecial('laplacian');
Data = [];
for i = 1:numel(Images)
    I = Images{i};
    std1 = stdfilt(I,ones(3));
    std2 = stdfilt(I,ones(11));
    lap = abs(imfilter(I,laph));
    lOg = imfilter(I,lOgh);
    sob1 = imfilter(I,sobh);
    sob2 = imfilter(I,sobh');
    sob = sqrt(sob1.^2+sob2.^2);    
    g = imfilter(I,h);
    Data_i = cat(2,I(:),g(:),std1(:),std2(:),lap(:),lOg(:),sob(:));
    Data = cat(1,Data,Data_i);
end
Data_std = std(Data,0,1);
Data_norm = bsxfun(@rdivide, Data, Data_std);
