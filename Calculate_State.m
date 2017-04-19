function states = Calculate_State(I,L,varargin)

[Height,Width]=size(I);
hcrop = Height/8;
wcrop = Width/8;
[X,Y]=meshgrid(1:Width,1:Height);
if length(varargin)==1
    Kalmans= varargin{1};
else
    Kalmans=[];
end


if isfield(Kalmans,'U')
   % profile on;
    enabled = logical([Kalmans.enabled]);
    enKalmans = Kalmans(enabled);
    U = arrayfun(@(x) fullSingle(x.U),enKalmans,'UniformOutput',false);
    %Us = arrayfun(@(x) sparseSingle(x.U),enKalmans,'UniformOutput',false);
    %Cell_Size = cellfun(@(u) sum(u(:).*I(:)),U,'UniformOutput',false);
    %MUXY = cellfun(@(u,cell_size)  u(:)'*[X(:),Y(:)]/sum(u(:)),U,Cell_Size,'UniformOutput',false);
    Cell_Size = cellfun(@(u) sum(u(:)),U,'UniformOutput',false);
    weightedSize =  cellfun(@(u,i) (I(:)'*u(:)),U,'UniformOutput',false);
    MUXY = cellfun(@(u,cell_size)  u(:)'*[X(:),Y(:)]/cell_size,U,Cell_Size,'UniformOutput',false);
    
    MU = cellfun(@(mu)  [mu,I(round(mu(2)),round(mu(1)))],MUXY,'UniformOutput',false);
    SIGMA = cellfun(@(u,mu,cell_size)  sqrt(u(:)'*(I(:)-mu(3)).^2/cell_size),U,MU,Cell_Size,'UniformOutput',false);
    BWs = arrayfun(@(x) sparseSingle(L==x.ID),enKalmans,'UniformOutput',false);
    CONTOURs = cellfun(@(bw) sparseSingle(fullSingle(bw) - imerode(fullSingle(bw),ones(3))) ,BWs,'UniformOutput',false);
    prev_CONTOUR = arrayfun(@(x) x.Contour,enKalmans,'UniformOutput',false);
    prev_STATE = arrayfun(@(x) x.prev_state,enKalmans,'UniformOutput',false);
    HDs = cellfun(@Calc_HD,CONTOURs,prev_CONTOUR,MU,prev_STATE,'UniformOutput',false);
    STATE = cellfun(@(mu,sigma,cell_size) [mu,sigma,cell_size,0,0,0,0,0],MU,SIGMA,Cell_Size,'UniformOutput',false);
    [states(1:length(STATE)).kalman_state] = deal(STATE{:});
    [states(:).HD] = deal(HDs{:});
    [states(:).BW] = deal(BWs{:});
    [states(:).Contour] = deal(CONTOURs{:});
    [states(:).size] = deal(Cell_Size{:});
    [states(:).weightedSize] = deal(weightedSize{:});
    
else
   uniqueL = unique(L(L>0));
   states = struct('ID',[],'kalman_state',[],'BW',[],'Contour',[]);
for i = 1:length(uniqueL);
    l = uniqueL(i);
    if l==0
        continue
    end
    states(i).ID = l;
    T = zeros(size(L));
    idx = L(:)==l;
    if isempty(idx)
      states(i).kalman_state =[];
      states(i).SDF =[];
      states(i).means=[];
      states(i).BW=L==l;
      continue;
    end
      
    T(idx)=I(idx);
    [y,x,g] = find(T);
    muxy = mean([x,y],1);
    mu = [muxy,I(round(muxy(2)),round(muxy(1)))];
    
    
    gl_sigma = std(g);
    BW = (L==l);
    cell_size = sum(BW(:));
    weightedSize =  I(:)'*BW(:);
    
    states(i).size = cell_size;
    states(i).BW = sparseSingle(BW);
    states(i).Contour = sparseSingle(BW-imerode(BW,ones(3)));
    states(i).kalman_state = [mu,gl_sigma,cell_size,0,0,0,0,0];
    states(i).weightedSize = weightedSize;
    
end
end

end

function hd = Calc_HD(Contour,prev_Contour,mu,prev_state)
        prevContour = fullSingle(prev_Contour);
        [p1y,p1x] = find(fullSingle(prev_Contour)); [p2y,p2x]=find(fullSingle(Contour));
        p2x = round(p2x - mu(1)+prev_state(1)) ;
        p2y = round(p2y -mu(2)+prev_state(2)) ;
        valid = (p2x>0&p2x<size(prevContour,2))&p2y>0&p2y<size(prevContour,1);
        p2x = p2x(valid);
        p2y = p2y(valid);
        p1 = [p1y,p1x];p2=[p2y,p2x];
        
        boxy = min([p1y;p2y]):max([p1y;p2y]);
       
        boxx = min([p1x;p2x]):max([p1x;p2x]);
        
        prevContour = prevContour(boxy,boxx);
        dist = bwdist(prevContour);
        ind = sub2ind(size(prevContour),(p2y-boxy(1)+1),(p2x-boxx(1)+1));
        
        hd = mean(dist(ind));
        
        %hd = HausdorffDist(p1,p2);
        if hd==0 
            hd = eps;
        end
        if isempty(hd)
            hd = 1e2;
        end
end

function Cropped = CropImage_prev(I,c,c_n,h,w,fillval)
yn1 = max(round(c_n(1)-h/2),1);
yn2 = min(round(c_n(1)+h/2),size(I,1));
xn1 = max(round(c_n(2)-w/2),1);
xn2 = min(round(c_n(2)+w/2),size(I,2));

Cropped = fillval.*ones(yn2-yn1+1,xn2-xn1+1);

y1 = max(round(c(1)-h/2),1);
y2 = min(round(c(1)+h/2),size(I,1));
x1 = max(round(c(2)-w/2),1);
x2 = min(round(c(2)+w/2),size(I,2));

Y1 = 1; Y2 = size(Cropped,1);
X1 = 1; X2 = size(Cropped,2);
if (y2-y1)<(yn2-yn1)
    if y1==1
        Y1 = (yn2-yn1)-(y2-y1)+1;
    else
        Y2 =(y2-y1)+1;
    end
elseif (y2-y1)>(yn2-yn1)
    if yn1==1
        y1 = y2-Y2+1;
    else
        y2 = y1 + Y2-1;
    end
end

if (x2-x1)<(xn2-xn1)
    if x1==1
        X1 = (xn2-xn1)-(x2-x1)+1;
    else
        X2 = (x2-x1)+1;
    end
elseif (x2-x1)>(xn2-xn1)
    if xn1==1
        x1 = x2-X2+1;
    else
        x2 = x1 + X2-1;
    end
end
Cropped(Y1:Y2,X1:X2) = I(y1:y2,x1:x2);
end
function PadI = PadImage(I,c,h,w,H,W,padv)
PadI = padv*ones(H,W);

y1 = max(round(c(1)-h/2),1);
y2 = min(round(c(1)+h/2),H);
x1 = max(round(c(2)-w/2),1);
x2 = min(round(c(2)+w/2),W);

if (y2-y1+1)==size(I,1)&&(x2-x1+1)==size(I,2)
    PadI(y1:y2,x1:x2) = I;
end

end
function Cropped = CropImage(I,c,h,w)
y1 = max(round(c(1)-h/2),1);
y2 = min(round(c(1)+h/2),size(I,1));
x1 = max(round(c(2)-w/2),1);
x2 = min(round(c(2)+w/2),size(I,2));
Cropped = I(y1:y2,x1:x2);
end

