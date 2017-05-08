function [L,L_New_Cell,Kalmans,z_pred,z_pred_orig,cog_diff,DEBUG] = Fuzzy_Segmentation(Tracking,Kalmans,I,I_prev,params,save_debug,t) %#ok<INUSD>
[Height,Width]=size(I);
[X,Y] = meshgrid(1:Width,1:Height);
conserveMemory = true;
global counter
XY = [Y(:),X(:)]';

if isfield(params,'patchSize')
    patchSize = params.patchSize;
else
    patchSize = 100;
end
if isfield(params,'minErr')
    min_err = params.minErr;
else
    min_err = 15;
end
if isfield(params,'maxItr')
    max_itr = params.maxItr;
else
    max_itr =5;
end

if isfield(params,'minCellSize')
    min_cell_size = params.minCellSize;
else
    min_cell_size = 100;
end

if isfield(params,'solidityThr')
    Solidity_Thr = params.solidityThr;
else
    Solidity_Thr = 0.8;
end

if isfield(params,'cellPrior')
    alpha = params.cellPrior;
else
    alpha = 0.2;
end
if isfield(params,'regImages')
    regImages = params.regImages;
else
    regImages = false;
end
if isfield(params,'edgeQ')
    q = params.edgeQ;
else
    q = 2;
end
if isfield(params,'edgeK')
    k = params.edgeK;
else
    k = 20;
end

if isfield(params,'geoOnly')
    geoOnly = params.geoOnly;
else
    geoOnly =false;
end

enabled = [Kalmans.enabled]'&([Kalmans.size]>=ceil(min_cell_size/10))';
enKalmans = Kalmans(enabled);
if save_debug
    DEBUG = cell(max_itr,1);
else
    DEBUG = [];
end
try
    if params.HDGMM
        nBG = calcGMMProb(I,Tracking.gmmFG, Tracking.gmmBG, alpha);
        
    else
        pGrayGlobal = alpha*Tracking.dens_cells./(alpha*(Tracking.dens_cells)+(1-alpha)*Tracking.dens_BG);
        %PBG = Tracking.dens_BG(round(I)+1);    
        %P_Cell = Tracking.dens_cells(round(I)+1);
        nBG = pGrayGlobal(round(I)+1);

        %PBG_prev = Tracking.dens_BG(round(I_prev)+1);
        %P_Cell_prev = Tracking.dens_cells(round(I_prev)+1);


        %nBG_prev = pGrayGlobal(round(I_prev)+1);
    end
    if regImages
    [optimizer, metric]  = imregconfig('monomodal');
    tform = imregtform(I, I_prev, 'translation', optimizer, metric);
    cog_diff = -tform.T(3,1:2);
    else
        cog_diff = [0 0];
    end
    
    [z_pred,~, P_pred] = arrayfun(@(x) predict(x.kalman),Kalmans,'UniformOutput',false);
    z_pred_orig = z_pred;
    z_pred_mat = cell2mat(shiftdim(z_pred(enabled),1));
    
    new_enabled = ~(z_pred_mat(:,1)<3|z_pred_mat(:,2)<3|z_pred_mat(:,1)>Width-3|z_pred_mat(:,2)>Height-3);
    enabled(enabled) = new_enabled;
    z_pred = z_pred(enabled);
    
    z_pred = cellfun(@(z) z+[cog_diff,0,0,0,0,0,0,0,0],z_pred,'UniformOutput',false);
    mu = cellfun(@(z) flipud(z(1:2)'), z_pred,'UniformOutput',false);
    g = cellfun(@(z) (z(3)), z_pred,'UniformOutput',false);
    est_mu = mu;
    mu_prev = mu;
    P_pred = P_pred(new_enabled);
    enKalmans = enKalmans(new_enabled);
    [enKalmans.z_pred] = deal(z_pred{:});
    cell_area = arrayfun(@(x) x.size,enKalmans,'UniformOutput',false);
    cellSize = arrayfun(@(x) x.weightedSize,enKalmans,'UniformOutput',false);
    BWs = arrayfun(@(x) x.BW,enKalmans,'UniformOutput',false);
    HDs = arrayfun(@(x) max(x.HD*5,0.01),enKalmans,'UniformOutput',false);
    last_mu = arrayfun(@(x) (flipud(x.prev_state(1:2))),enKalmans,'UniformOutput',false);
    hg = fspecial('gauss');
    gradI = calc_Grad_I(I);
    invG = max(1./(1+(gradI./(k*std(gradI(:)))).^q),0.01);
    localGray = false;
    
    if isfield(enKalmans,'dens_cells')
        pGrayLocal = arrayfun(@(x) alpha*x.dens_cells./((1-alpha)*x.dens_BG+alpha*x.dens_cells),enKalmans,'UniformOutput',false);
        emptyGrayLocal = cellfun(@isempty,pGrayLocal);
        [pGrayLocal{emptyGrayLocal}] = deal(pGrayGlobal);
        localGray = true;
    end
    
    if conserveMemory
        clear gradI;
    end
    sigG = P_pred{1}(3,3);

    %%
    
    %%
    
    
    itr = 0;
    err = min_err+1;
    
    
    %%
    NKalmans = sum(enabled);
    U = cell(1,NKalmans);

    SDF_moved = cell(1,NKalmans);

    Phi_norm =cell(1,NKalmans);
    Phi_cropped =cell(1,NKalmans);
    gradI_cropped = cell(1,NKalmans);
    invG_cropped = cell(1,NKalmans);
    nBG_cropped = cell(1,NKalmans);
    mu_cropped = cell(1,NKalmans);
    Speed_cropped = cell(1,NKalmans);
    D_cropped =cell(1,NKalmans);
    P_cropped = cell(1,NKalmans);
    P = cell(1,NKalmans);
    U_cropped =cell(1,NKalmans);
    
    while err>=min_err&&itr<max_itr
        
        
        itr=itr+1;
        if save_debug
            DEBUG{itr}.itr = itr;
            DEBUG{itr}.est_mu = est_mu;
            DEBUG{itr}.P_pred = P_pred;
            DEBUG{itr}.last_mu = last_mu;
            DEBUG{itr}.I = I;
            DEBUG{itr}.nBG = nBG;
            DEBUG{itr}.invG = invG;
        end
        if itr==1
            %Gauss = cellfun(@(m,p) reshape(mvnpdf([Y(:),X(:)],m',rot90(p(1:2,1:2)*10,2)),size(X)),est_mu,P_pred,'uniformoutput',false);
       
            Dist = cellfun(@(m,p) calcDistMap(X,Y,m,p(1)*10,patchSize) ,est_mu,P_pred,'uniformoutput',false);
            sumDist =reshape(sum(cat(2,Dist{:}),2)+eps,size(X));
            Dist_stat = cellfun(@(m,p) calcDistMap(X,Y,m,p(1)*10,patchSize) ,last_mu,P_pred,'uniformoutput',false);
            sumDist_stat =reshape(sum(cat(2,Dist_stat{:}),2)+eps,size(X));
            last_mu_orig = last_mu;
            %mu = cellfun(@(m,p) FindMostLiklymuGray(nBG,invG,m,p),mu,P_pred,'uniformoutput',0);
            mu = cellfun(@(d) FindMostLiklymuGray2(nBG,invG,reshape(d.^2,size(X)),sumDist,ones(size(X))),Dist,'uniformoutput',0);
            last_mu = cellfun(@(d) FindMostLiklymuGray2(nBG,invG,reshape(d.^2,size(X)),sumDist_stat,ones(size(X))),Dist_stat,'uniformoutput',0);
            
            if conserveMemory
                clear Dist;
            end
        end
        
        counter = 0;
        [Phi_moved,~,~,~,~,valid] = cellfun(@(bw,em,mm,l,sigma,hd) movePhi(fullSingle(bw),em,mm,l,hd,sigma,patchSize),BWs,est_mu,mu,last_mu_orig,P_pred,HDs,'uniformoutput',0);
        [Phi_stat,~,~,~,~,valid_stat] = cellfun(@(bw,em,mm,l,sigma,hd) movePhi(fullSingle(bw),em,mm,l,hd,sigma,patchSize),BWs,last_mu,last_mu,last_mu,P_pred,HDs,'uniformoutput',0);
        valid = [valid{:}];
        valid_stat = [valid_stat{:}];
        valid = valid_stat&valid;
        
        if any(~valid)
            
            HDs = HDs(valid);
            Phi_moved = Phi_moved(valid);
            Phi_stat = Phi_stat(valid);
            %BWs_moved = BWs_moved(valid);
            BWs = BWs(valid);
            %BWs_stat = BWs_stat(valid);
            %BW_cropped = BW_cropped(valid);
            %SDF_stat_cropped = SDF_stat_cropped(valid);
            
            %BW_cropped_stat = BW_cropped_stat(valid);
            %SDF_moved_cropped = SDF_moved_cropped(valid);
            %Phi_moved_cropped = Phi_moved_cropped(valid);
            %Phi_stat_cropped = Phi_stat_cropped(valid);
            mu = mu(valid);
            %mu_prev = mu_prev(valid);
            est_mu = est_mu(valid);
            last_mu = last_mu(valid);
            last_mu_orig = last_mu_orig(valid);
            cell_area = cell_area(valid);
            enKalmans = enKalmans(valid);
            P_pred = P_pred(valid);
            if localGray
            pGrayLocal = pGrayLocal(valid);
            end
            %Pi = Pi(valid);
            g = g(valid);
            cellSize = cellSize(valid);
            enabled(enabled) = logical(valid);
            if exist('U','var')
                U = U(valid);
            end
            
            
            
        end
        
        %Phimat = reshape(full(cell2mat(Phi_moved)),Height,Width,[]);
        %Sphi = sum(Phimat,3) + (1-max(Phimat,[],3));
        %Sphi(Sphi(:)==0)=1;
        if save_debug
            continue;   
        end
        if itr==1
            %Phi_norm = cellfun(@(phi) phi./Sphi,Phi_moved,'uniformoutput',0);
            %Phi_cropped = cellfun(@(Im,m) CropImage(Im,m,patchSize,patchSize),Phi_norm,mu,'uniformoutput',0);
            
            [Phi_cropped, cropInd] = cellfun(@(Im,m) CropImage(Im,m,patchSize,patchSize),Phi_moved,mu,'uniformoutput',0);
            Phi_cropped_stat = cellfun(@(Im,m) CropImage(Im,m,patchSize,patchSize),Phi_stat,mu,'uniformoutput',0);
            if save_debug
                DEBUG{itr}.Phi_cropped = Phi_cropped;
            end
            invG_cropped = cellfun(@(m) CropImage(invG,m,patchSize,patchSize),mu,'uniformoutput',0);
            I_cropped = cellfun(@(m) CropImage(I,m,patchSize,patchSize),mu,'uniformoutput',0);
            if localGray
                nBG_cropped = cellfun(@(p,Ic) p(round(Ic)+1),pGrayLocal,I_cropped,'uniformoutput',0);
            else
                nBG_cropped = cellfun(@(m) CropImage(nBG,m,patchSize,patchSize),mu,'uniformoutput',0);
            end
            %Pi_cropped = cellfun(@(m,p) CropImage(p,m,patchSize,patchSize),mu,Pi,'uniformoutput',0);
            
            
            PGray_cropped = cellfun(@(gg,i,csize,area)  reshape(sqrt(2*pi)*csize./area*normpdf(i(:),gg,csize./area),size(i)), g,I_cropped,cellSize,cell_area,'uniformoutput',false);
    
            %nEdge_cropped = cellfun(@(m) CropImage(nEdge,m,patchSize,patchSize),mu,'uniformoutput',0);
            if save_debug
                DEBUG{itr}.nBG_cropped = nBG_cropped;
            end
            %mu_cropped = cellfun(@FindMostLiklymu,nBG_cropped,Phi_cropped,'uniformoutput',0);
            mu_cropped =  cellfun(@(m)CropMu(m,patchSize,patchSize),mu,'uniformoutput',false);
            mu_cropped_stat =  cellfun(@(mc,m,lm) (round(mc-m+lm)) ,mu_cropped,mu,last_mu,'uniformoutput',false);
            counter = 0;
            mu_cropped = cellfun(@FindMostLiklymuGray,nBG_cropped,invG_cropped,mu_cropped,P_pred,PGray_cropped,'uniformoutput',0);
            if conserveMemory
                %clear PGray_cropped
            end
            tmp = cellfun(@isempty,mu_cropped);
            if any(tmp)
                save('TMP.mat');
                error('Empty mu_corpped');
            end
            if save_debug
                DEBUG{itr}.mu_cropped = mu_cropped;
            end
            Speed_cropped = cellfun(@(nbg,invg,phi) max((nbg).*(invg),1e-8),nBG_cropped,invG_cropped,Phi_cropped,'uniformoutput',0);
           % Speed_cropped = cellfun(@fixSpeedIm,Speed_cropped,mu_cropped,P_pred,'uniformoutput',0);
            if save_debug
                DEBUG{itr}.Speed_cropped= Speed_cropped;
            end
             if conserveMemory
                clear invG_cropped;
            end
            
            if ~geoOnly
            [D_cropped,Dg,De] = cellfun(@CreateDiffDist,Speed_cropped,mu_cropped,'uniformoutput',0);
            D_cropped_stat = cellfun(@CreateDiffDist,Speed_cropped,mu_cropped_stat,'uniformoutput',0);
            else
            [~,D_cropped] = cellfun(@CreateDiffDist,Speed_cropped,mu_cropped,'uniformoutput',0);
            [~,D_cropped_stat] = cellfun(@CreateDiffDist,Speed_cropped,mu_cropped_stat,'uniformoutput',0);
                
            end
            if save_debug
                DEBUG{itr}.D_cropped = D_cropped;
            end
            P_cropped = cellfun(@(phi,d) 1./(d+1),Phi_cropped,D_cropped,'uniformoutput',0);
            P_cropped_tot = cellfun(@(phi,d,p) (p+1./(d+1))*0.5,Phi_cropped_stat,D_cropped_stat,P_cropped,'uniformoutput',0);
            if conserveMemory
                clear D_cropped_stat;
                clear D_cropped;
                clear P_cropped
                clear Phi_cropped
                clear Phi_cropped_stat
            end
            
            %P = cellfun(@(Im,m) PadImage(Im,m,patchSize,patchSize,Height,Width,0),P_cropped_tot,mu,'uniformoutput',0);
            %P_Weighted = cellfun(@(p,nbg,i,cellsize) RegularizeIntensity(p,ones(size(p)),nbg,i,cellsize,1.5),P_cropped_tot,nBG_cropped,I_cropped,cellSize,'uniformoutput',false);
            
            %P = cellfun(@(Im,m) PadImage(Im,m,patchSize,patchSize,Height,Width,0),P_cropped_tot,mu,'uniformoutput',0);
            sparseP = cellfun(@(p,ind) sparse(ind,ones(size(ind)),p(:),Height*Width,1),P_cropped_tot,cropInd,'uniformoutput',false); 
            Pmat = cat(2,sparseP{:});
            
            [BGy,BGx] = find(nBG<0.5);
        rng(1);
        %r = randperm(numel(BGx),ceil(numel(BGx)./10));
        perim = zeros(size(nBG));
        perim(1:2,:) = 1;
        perim(end-1:end,:) = 1;
        perim(:,1:2) = 1;
        perim(:,end-1:end) = 1;
        [perimy,perimx] = find(perim); 
        %SE = (1-nBG).*invG;
        %DBGg = max(graydist(1./SE,BGx(r),BGy(r),'quasi-euclidean'),eps);
        %DBGe = max(graydist(ones(size(SE)),BGx(r),BGy(r),'quasi-euclidean'),eps);
        SE = nBG.*nBG>0.5;
        DBGg = max(graydist(1./SE,perimx,perimy,'quasi-euclidean'),eps);
        DBGe = max(double(bwdist(perim,'quasi-euclidean')),eps);
        
        if ~geoOnly
             DBGDiff= DBGg-DBGe;
        else
             DBGDiff= DBGg;
        end
        pBG = 1./(1+DBGDiff);
        else
            %Phi_norm = cellfun(@(phi) phi./Sphi,Phi_moved,'uniformoutput',0);
            %Phi_cropped = cellfun(@(Im,m) CropImage(Im,m,patchSize,patchSize),Phi_norm,mu,'uniformoutput',0);
            
            [Phi_cropped cropInd] = cellfun(@(Im,m) CropImage(Im,m,patchSize,patchSize),Phi_moved,mu,'uniformoutput',0);
             Phi_cropped_stat = cellfun(@(Im,m) CropImage(Im,m,patchSize,patchSize),Phi_stat,mu,'uniformoutput',0);
            I_cropped = cellfun(@(m) CropImage(I,m,patchSize,patchSize),mu,'uniformoutput',0);
            if save_debug
                DEBUG{itr}.Phi_cropped = Phi_cropped;
            end
            invG_cropped = cellfun(@(m) CropImage(invG,m,patchSize,patchSize),mu,'uniformoutput',0);
            if localGray
                nBG_cropped = cellfun(@(p,Ic) p(round(Ic)+1),pGrayLocal,I_cropped,'uniformoutput',0);
            else
                nBG_cropped = cellfun(@(m) CropImage(nBG,m,patchSize,patchSize),mu,'uniformoutput',0);
            end
            if save_debug
                DEBUG{itr}.nBG_cropped = nBG_cropped;
            end
            U_cropped = cellfun(@(Im,m) CropImage(reshape(Im,Height,Width),m,patchSize,patchSize),U,mu,'uniformoutput',0);
            %mu_cropped = cellfun(@(u,phi)FindMuCOM(u.*phi),U_cropped,Phi_cropped,'uniformoutput',0);
            mu_cropped =  cellfun(@(m)CropMu(m,patchSize,patchSize),mu,'uniformoutput',false);
            if itr==2
                PGray_cropped = cellfun(@(gg,i,csize,area)  reshape(sqrt(2*pi)*csize./area*normpdf(i(:),gg,csize./area),size(i)), g,I_cropped,cellSize,cell_area,'uniformoutput',false);
                mu_cropped = cellfun(@FindMostLiklymuGray,nBG_cropped,invG_cropped,mu_cropped,P_pred,PGray_cropped,'uniformoutput',0);
            end
                
            mu_cropped_stat =  cellfun(@(mc,m,lm) (round(mc+-m+lm)) ,mu_cropped,mu,last_mu,'uniformoutput',false);
            if save_debug
                DEBUG{itr}.mu_cropped = mu_cropped;
            end
            Speed_cropped = cellfun(@(nbg,invg,phi) max((nbg).*(invg),1e-8),nBG_cropped,invG_cropped,Phi_cropped,'uniformoutput',0);
            Speed_cropped_stat = cellfun(@(nbg,invg,phi) max((nbg).*(invg),1e-8),nBG_cropped,invG_cropped,Phi_cropped_stat,'uniformoutput',0);
             if conserveMemory
                 clear invG_cropped
             end
            %Speed_cropped = cellfun(@fixSpeedIm,Speed_cropped,mu_cropped,P_pred,'uniformoutput',0);
            if save_debug
                DEBUG{itr}.Speed_cropped= Speed_cropped;
            end
            if ~geoOnly
            [D_cropped, ~,De] = cellfun(@CreateDiffDist,Speed_cropped,mu_cropped,'uniformoutput',0);
            [D_cropped_stat, ~, De_stat] = cellfun(@CreateDiffDist,Speed_cropped_stat,mu_cropped_stat,'uniformoutput',0);
            else
                [~,D_cropped,De] = cellfun(@CreateDiffDist,Speed_cropped,mu_cropped,'uniformoutput',0);
                [~,D_cropped_stat] = cellfun(@CreateDiffDist,Speed_cropped,mu_cropped_stat,'uniformoutput',0);
            
            end
             if conserveMemory
                 clear Speed_cropped
                 clear Speed_cropped_stat
             end
            if save_debug
                DEBUG{itr}.D_cropped = D_cropped;
            end
            if itr >2
            P_cropped = cellfun(@(phi,d) phi./(d+1),Phi_cropped,D_cropped,'uniformoutput',0);
            P_cropped_stat = cellfun(@(phi,d) phi./(d+1),Phi_cropped_stat,D_cropped_stat,'uniformoutput',0);
            else
            P_cropped = cellfun(@(phi,d) 1./(d+1),Phi_cropped,D_cropped,'uniformoutput',0);
            P_cropped_stat = cellfun(@(phi,d) 1./(d+1),Phi_cropped_stat,D_cropped_stat,'uniformoutput',0);   
            end
            P_cropped_tot = cellfun(@createPTot ,P_cropped,P_cropped_stat,U_cropped,mu_cropped,mu_cropped_stat,'uniformoutput',false);
            
            if conserveMemory
                
            end
            
            %P_cropped_tot = cellfun(@(phi,d,p,u,m,ms) (p.*u(m(1),m(2))+u(ms(1),ms(2))./(d+1))./(u(m(1),m(2))+u(ms(1),ms(2))),Phi_cropped_stat,D_cropped_stat,P_cropped,U_cropped,mu_cropped,mu_cropped_stat,'uniformoutput',0);
            %P = cellfun(@(Im,m) PadImage(Im,m,patchSize,patchSize,Height,Width,0),P_cropped_tot,mu,'uniformoutput',0);
            
            counter = 0;
            alpha = min([cellSize{~isinf([cellSize{:}])}]);
            %alpha = 1.2;
            P_Weighted = cellfun(@(p,de,u, nbg, i, cellsize) RegularizeIntensity(p,u,de,nbg,i,cellsize,alpha),P_cropped,De,U_cropped,nBG_cropped,I_cropped,cellSize,'uniformoutput',false);
            P_Weighted_stat = cellfun(@(p,de,u, nbg, i, cellsize) RegularizeIntensity(p,u,de,nbg,i,cellsize,alpha),P_cropped_stat,De,U_cropped,nBG_cropped,I_cropped,cellSize,'uniformoutput',false);
            P_cropped_tot = cellfun(@createPTot ,P_Weighted,P_Weighted_stat,U_cropped,mu_cropped,mu_cropped_stat,'uniformoutput',false);
            
            %P = cellfun(@(Im,m) PadImage(Im,m,patchSize,patchSize,Height,Width,0),P_Weighted,mu,'uniformoutput',0);
            sparseP = cellfun(@(p,ind) sparse(ind,ones(size(ind)),p(:),Height*Width,1),P_cropped_tot,cropInd,'uniformoutput',false); 
            Pmat = cat(2,sparseP{:});
            if conserveMemory
                clear P_Weighted;
                clear I_cropped;
                clear D_cropped_stat;
                clear D_cropped;
                clear P_cropped;
                clear P_cropped_stat;
                clear Phi_cropped;
                clear Phi_cropped_stat;
            end
            
        end
        
        %Pmat = reshape(full(cell2mat(P)),Height,Width,[]);
        
        
        %BGSkel = bwmorph(nBG<0.5,'skel',inf);
        S = max(sum(Pmat,2) + 1e-4*pBG(:),eps);

        

        %% This part is new June 19th 
        [U,U_cropped] = cellfun(@(p,m) NormalizeAndCropU(p,S,nBG,m,patchSize),sparseP,mu,'uniformoutput',0);
        
        %U_cropped = cellfun(@corU, U_cropped,Phi_cropped,'uniformoutput',false);
        %U  = cellfun(@(Im,m) PadImage(Im,m,patchSize,patchSize,Height,Width,0),U_cropped,mu,'uniformoutput',0);
        %% END
        UBG = 1-nBG;
        if save_debug
            DEBUG{itr}.UBG = UBG;
            DEBUG{itr}.P_cropped = P_cropped;
            DEBUG{itr}.Phi_cropped = Phi_cropped;
            DEBUG{itr}.U_cropped = U_cropped;
            DEBUG{itr}.mu = mu;
            DEBUG{itr}.Speed_cropped= Speed_cropped;
            DEBUG{itr}.D_cropped = D_cropped;
            DEBUG{itr}.pBG = pBG;
            DEBUG{itr}.enKalmans = enKalmans;
        else
            DEBUG = [];
        end
        mu_prev = mu;
        
        if itr>1
        mu = cellfun(@(u,phi,phis) XY(:,(nBG(:)>0.5))*(u(nBG(:)>0.5))./max(sum(u(nBG(:)>0.5)),eps),U,Phi_moved,Phi_stat,'UniformOutput',false);
        end
        prev_size = cell_area;
        cellSize = cellfun(@(u) round(I(:)'*u(:)),U,'uniformoutput',false);
        cell_area = cellfun(@(u) sum(u(:)),U,'uniformoutput',false);
        [err,err_idx] = max(sqrt(sum((cell2mat(cell_area)-cell2mat(prev_size)).^2,1)));
        if any([cell_area{:}]==0)
            remove_cell = full([cell_area{:}]>0);
            BWs = BWs(remove_cell);
            %BWs_stat = BWs_stat(remove_cell);
            %BW_cropped_stat = BW_cropped_stat(remove_cell);
            P_pred = P_pred(remove_cell);
            HDs = HDs(remove_cell);
            %BWs_moved = BWs_moved(remove_cell);
            mu = mu(remove_cell);
            mu_prev = mu_prev(remove_cell);
            est_mu = est_mu(remove_cell);
            SDF_moved = SDF_moved(remove_cell);
            Phi_moved = Phi_moved(remove_cell);
            Phi_stat = Phi_stat(remove_cell);
            Phi_norm = Phi_norm(remove_cell);
            if localGray
            pGrayLocal = pGrayLocal(remove_cell);
            end

            U = U(remove_cell);
            last_mu = last_mu(remove_cell);
            last_mu_orig = last_mu_orig(remove_cell);
            
            [enKalmans(~remove_cell).enabled] = deal(false);
            enabled(enabled) = logical(remove_cell);
            enKalmans= enKalmans(remove_cell);
            cell_area = cell_area(remove_cell);
            prev_size = prev_size(remove_cell);
            g = g(remove_cell);
            cellSize = cellSize(remove_cell);
        end
        
    end

    Us = cellfun(@(u,phi,phis) sparseSingle(u),U,Phi_moved,Phi_stat,'UniformOutput',false);
    [enKalmans(:).U] = deal(Us{:});
    
    Kalmans(enabled) = enKalmans;
    enabledcell = num2cell(logical(enabled));
    [Kalmans(:).enabled] = deal(enabledcell{:});

    
    [~,L] = max(cat(2,1-nBG(:),U{:}),[],2);
    L = reshape(L,Height,Width);

    L = L-1;
    labels = [enKalmans.ID];
    L(L>0) = labels(L(L>0));

    L_cropped = cellfun(@(m) CropImage(L,m,patchSize,patchSize),mu,'uniformoutput',0);
    NotBG = false(size(I));
    nBGinvG = nBG.*invG;
    NotBG(L==0) = logical(nBG(L==0)>0.5);

    NotBGFull = NotBG;
    NotBG_Candidates = [];
    links = {};
    strelOpen = strel('disk',4);
    for l = unique(L(L>0))'
        BW = L==l;
        BW = imopen(BW,strelOpen);
        L(L==l&~BW) = 0;
        if ~any(BW(:))
            continue;
            Kalmans([Kalmans.ID]==l).enabled=false;
        end
        NotBG_Candidates_l = regionprops(BW,'Solidity','Area','PixelIdxList','Centroid');
        cs = fliplr(cat(1,NotBG_Candidates_l(:).Centroid));
        k = find([enKalmans.ID]==l);
        ml = mu_prev{k}';
        r = sum((bsxfun(@minus,cs,ml)).^2,2)>25*det(P_pred{k}(1:2,1:2));
        BW_cropped = CropImage(BW,mu_prev{k},patchSize,patchSize);
        candidates_cropped = regionprops(BW_cropped,'Solidity','Area','PixelIdxList','Centroid');
        s = size(BW_cropped);
        %Solidity_Thr = min(1,Solidity_Thr*1.2);
        removeCandidate = NotBG_Candidates_l(([NotBG_Candidates_l.Solidity]<=Solidity_Thr|[NotBG_Candidates_l.Area]<=min_cell_size)|r');
        L(cat(1,removeCandidate(:).PixelIdxList)) = 0;
        NotBG(cat(1,removeCandidate(:).PixelIdxList)) = 1;
        NotBGFull(cat(1,removeCandidate(:).PixelIdxList)) = 1;
        candidates_cropped = candidates_cropped([candidates_cropped.Solidity]>Solidity_Thr&[candidates_cropped.Area]>min_cell_size);
        NotBG_Candidates_l =  NotBG_Candidates_l([NotBG_Candidates_l.Solidity]>Solidity_Thr&[NotBG_Candidates_l.Area]>min_cell_size);
        if length(NotBG_Candidates_l)==0
            Kalmans([Kalmans.ID]==l).enabled=false;
        end
        if length(NotBG_Candidates_l)<=1
            continue
        end
        flags = false(length(NotBG_Candidates_l),1);
        if length(NotBG_Candidates_l)>1
        for i = 1:length(NotBG_Candidates_l)
            [idxy,idxx]= ind2sub([Height,Width],NotBG_Candidates_l(i).PixelIdxList);
            [idxyC,idxxC]= ind2sub(s,candidates_cropped(i).PixelIdxList);
            if any(idxx==1)||any(idxy==1)||any(idxx==Width)||any(idxy==Height)||any(idxxC==1)||any(idxyC==1)||any(idxxC==s(2))||any(idxyC==s(1))
                flags(i) = true;
                NotBG_Candidates = cat(1,NotBG_Candidates,NotBG_Candidates_l(i));
                NotBGFull(NotBG_Candidates_l(i).PixelIdxList) = 1;
                L(NotBG_Candidates_l(i).PixelIdxList) = 0;
            end
        end
        end

        if sum(~flags)>1
            NotBG_Candidates = cat(1,NotBG_Candidates,NotBG_Candidates_l(~flags));
            NotBGFull(L==l) = 1;
            L(L==l) = 0;
            Kalmans([Kalmans(:).ID]==l).enabled = false;
            NotBG_Candidates_l = NotBG_Candidates_l(~flags);
            Kalmans([Kalmans(:).ID]==l).Children = NotBG_Candidates_l;
            
           
        end
    end
    NotBG = imerode(NotBG,ones(3));
    NotBG_Candidates_bg = regionprops(NotBG,'Solidity','Area','PixelIdxList','Centroid');
    NotBG_Candidates = cat(1,NotBG_Candidates,NotBG_Candidates_bg);
    New_Cells = NotBG_Candidates([NotBG_Candidates.Solidity]>Solidity_Thr&[NotBG_Candidates.Area]>min_cell_size);
    if ~isempty(New_Cells)
        BW_New_Cell = arrayfun(@(x) Create_New_Cell_L(NotBGFull,x.PixelIdxList),New_Cells,'UniformOutput',false);
        BWMAT = cat(3,BW_New_Cell{:});
        indBW = find(BWMAT);
        [indy,indx,indz] = ind2sub(size(BWMAT),indBW);
        L_New_Cell = zeros(size(L));
        L_New_Cell(sub2ind(size(L),indy,indx))=indz;
    else
        L_New_Cell=[];
    end
    killed = setdiff(unique(L(L>0)),[enKalmans.ID]');
    
catch err
    save('Debug.mat','DEBUG','-v7.3');
    rethrow(err);
end

end

function BW = Create_New_Cell_L(BW_Full,Idx)
BW = zeros(size(BW_Full));
BW(Idx(:)) = BW_Full(Idx(:));
%BW = imdilate(BW,ones(3));
end

function Cropped_mu = CropMu(c,h,w) %#ok<*DEFNU>
y1 = max(round(c(1)-h/2),1);
x1 = max(round(c(2)-w/2),1);
Cropped_mu = round(c-[y1-1;x1-1]);
end

function mu = unCropMu(c,mu_cropped,h,w)

y1 = max(round(c(1)-h/2),1)-1;
%y2 = min(round(c(1)+h/2),H);
x1 = max(round(c(2)-w/2),1)-1;
%x2 = min(round(c(2)+w/2),W);
mu = mu_cropped+[y1;x1];
end


function PadI = PadImage(I,c,h,w,H,W,padv)
if padv==0
    PadI = sparse(H,W);
else
    
    PadI = padv*ones(H,W);
end
y1 = max(round(c(1)-h/2),1);
y2 = min(round(c(1)+h/2),H);
x1 = max(round(c(2)-w/2),1);
x2 = min(round(c(2)+w/2),W);

if (y2-y1+1)==size(I,1)&&(x2-x1+1)==size(I,2)
    PadI(y1:y2,x1:x2) = I;
end

end

function gradI = calc_Grad_I(I)

hy = fspecial('sobel'); %[-1;0;1];
hx = hy';
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');
gradI = (Ix.^2 + Iy.^2);
end

function mu = FindMostLiklymu(PG,Phi)
%[X,Y] = meshgrid(1:size(PG,2),1:size(PG,1));
if size(PG)==size(Phi)
    [muy,mux]=find((PG.*Phi)==max(PG(:).*Phi(:)),1);
    %mu = round([Y(:),X(:)]'*(PG(:).*Phi(:))./sum(PG(:).*Phi(:)));
    mu = [muy;mux];
else
    [muy,mux]=find(PG==max(PG(:)),1);
    mu = [muy;mux];
end

end

function mu = FindMostLiklymuGray(PG,invG,mu,sigma,PGray)
global counter
counter = counter+1;
[X,Y] = meshgrid(1:size(PG,2),1:size(PG,1));
G = reshape(mvnpdf([Y(:),X(:)],mu',rot90(10*sigma(1:2,1:2),2)),size(X));
[~,ind] = max(PGray(:).*PG(:).*invG(:).*G(:));
[muy mux]=ind2sub(size(PG),ind);
%[muy,mux]=find((PG.*invG.*G.*pi)==max(pi(:).*PG(:).*invG(:).*G(:)),1);
%mu = round([Y(:),X(:)]'*(PG(:).*Phi(:))./sum(PG(:).*Phi(:)));
mu = [muy;mux];

end
function mu = FindMostLiklymuGray2(PG,invG,G,Gsum,Pi)
mapLikely = Pi.*PG.*invG.*G./Gsum;
[~,ind] = max(mapLikely(:));
[muy,mux] = ind2sub(size(G),ind);
%[muy,mux]=find(()==max(Pi(:).*PG(:).*invG(:).*G(:)./Gsum(:)),1);
%mu = round([Y(:),X(:)]'*(PG(:).*Phi(:))./sum(PG(:).*Phi(:)));
mu = [muy;mux];

end

function mu = FindMuCOM(P)
[X,Y] = meshgrid(1:size(P,2),1:size(P,1));
mu = round([Y(:),X(:)]'*P(:)./sum(P(:)));
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

function [DiffD,Dg,De] = CreateDiffDist(SE,mu)
%[X,Y] = meshgrid(1:size(SE,2),1:size(SE,1));
%De = max(sqrt((Y-mu(1)).^2+(X-mu(2)).^2),eps);
if any(isnan(SE(:)))||any(isnan(mu))||any(mu<1)||mu(1)>size(SE,1)||mu(2)>size(SE,2)
    DiffD = inf(size(SE));
    Dg = inf(size(SE));
    De = inf(size(SE));
    return
    %fprintf('Fast Marching will crash!')
end
Dg = max(graydist(1./SE,mu(2),mu(1),'quasi-euclidean'),eps);
BW = zeros(size(SE));
BW(mu(1),mu(2)) = 1;
De = max(double(bwdist(BW,'quasi-euclidean')),eps);
%Dg = max(msfm2d(SE,mu,true,true),eps);
DiffD = max(Dg-De,0);
end

function [Speed_G] = fixSpeedIm(Speed,mu,C)
[X,Y] = meshgrid(1:size(Speed,2),1:size(Speed,1));
D = (X(:)-mu(2)).^2./C(1,1)+(Y(:)-mu(1)).^2./C(2,2);
P = reshape(exp(-0.5*D),size(Speed));
Speed_G = max(Speed,P);
end

function [Phi_smoothed] = smoothPhi(Phi,C)
w = max(C(1,1),C(1,2))*5;
[X,Y] = meshgrid(1:w,1:w);
mu = round(size(X)./2);
D = (X(:)-mu(2)).^2./C(1,1)+(Y(:)-mu(1)).^2./C(2,2);
P = reshape(exp(-0.5*D),size(X));
P = P./sum(P(:));
Phi_smoothed = imfilter(Phi,P,'same','conv');
end

function [bw] = moveBW(bw,m,lastm)

[y,x] = find(bw);
y = round(y-lastm(1)+m(1));
x = round(x-lastm(2)+m(2));
invalid = (y<1|y>size(bw,1))|(x<1|x>size(bw,2));
y(invalid) = [];
x(invalid) = [];
idx = sub2ind(size(bw),y,x);
bw(:) = 0;

bw(idx) = 1;
end

function [Phi_moved,BWs_moved,BW_cropped,SDF_moved_cropped,Phi_moved_cropped,valid] = movePhi(BWs,est_mu,meas_mu,last_mu,HD,sigma,patchSize)
global counter
counter = counter+1;
BWs_moved=  moveBW(BWs,meas_mu,last_mu);
[H,W] = size(BWs_moved);
BW_cropped= CropImage(BWs_moved,meas_mu,patchSize,patchSize);
mu_cropped =  CropMu(meas_mu,patchSize,patchSize);
[X,Y] = meshgrid(1:size(BW_cropped,2),1:size(BW_cropped,1));
%G = reshape(mvnpdf([X(:),Y(:)],mu_cropped',sigma(1:2,1:2)),size(X));
SDF_moved_cropped = SDF(BW_cropped);
Phi_moved_cropped = 1./(1+exp(-SDF_moved_cropped./(2.*pi).*sqrt(3)));
%PhiG = imfilter(G,double(Phi_moved_cropped),'same');
PhiG = Phi_moved_cropped;%.*mvnpdf(meas_mu',est_mu',sigma(1:2,1:2))./mvnpdf(est_mu',est_mu',sigma(1:2,1:2));
Phi_moved = PadImage(PhiG,meas_mu,patchSize,patchSize,H,W,0);
%Phi_moved = PadImage(Phi_moved_cropped,meas_mu,patchSize,patchSize,H,W,0);
valid = any(BW_cropped(:));
end

function [U,U_cropped] = NormalizeAndCropU(P,S,nBG,mu,patchSize)
U = P./S.*nBG(:);
U_cropped= CropImage(reshape(U,size(nBG)),mu,patchSize,patchSize);
end

function ptot = createPTot(p,ps,u,m,ms)
m = round(m);
ms = round(ms);
if all(m'>1 & m'<size(u))
    f1 = u(m(1),m(2));
else
    f1 = 0;
end
if all(ms'>1 & ms'<size(u))
    f2 = u(ms(1),ms(2));
else
    f2 = 0;
end
if f1<eps&&f2<eps, f1=0.5;f2=0.5; end
ptot = (p*f1+ps*f2)./(f1+f2);
end


function mu = calcMu(mu_old,u,phi,phis,h,w)
[X,Y] = meshgrid(1:size(u,2),1:size(u,1));
YX = [Y(:),X(:)];
s = u(:).*(phis(:)+phi(:));
mu_cropped = YX'*s./max(sum(s),eps);
mu = unCropMu(mu_old,mu_cropped,h,w);

end


function uout = corU(u,phi)
cor = filter2(phi./sum(phi(:)),u);
[peaky, peakx] = find(cor==max(cor(:))); 
peak = zeros(size(cor));
peak(peaky,peakx)=1;
uout = filter2(flipud(fliplr(phi)),peak);
end

function [PWeighted,PixelWeights] = RegularizeIntensity(P_cropped,U_cropped,De,nBG_cropped,I_cropped,cellSize,alpha)
global counter
counter = counter+1;
%[~,sid] = sort(U_cropped(:)./(1+De(:)),'descend');
[~,sid] = sort(P_cropped(:),'descend');

sortedSum = zeros(size(P_cropped));
sortedSum(sid) = cumsum(I_cropped(sid).*(nBG_cropped(sid)));
sortedSum = reshape(sortedSum,size(P_cropped));
%PixelWeights = max(1e-4,min(1,(sortedSum./((1-alpha).*cellSize)-(alpha./(1-alpha)))));
if isinf(cellSize)
    cellSize = sortedSum;
end
PixelWeights = max(0,(sortedSum - cellSize)./alpha+1);

%PixelWeights = max(eps,min(1,(-(sortedSum-cellSize)./(alpha./10))+1));
%PWeighted = PixelWeights.*P_cropped;
PWeighted = (0.5*P_cropped).^PixelWeights;
end
            
function distMap = calcDistMap(X,Y,mu,sig,patchSize)
Y = (Y(:)- mu(1));
X = (X(:) - mu(2));
infLoc = X.^2>patchSize.^2|Y.^2>patchSize.^2;
G = mvnpdf([X(~infLoc),Y(~infLoc)],[0,0],diag([sig,sig]));
%X(infLoc) = inf;
%Y(infLoc) = inf;
distMap = sparse(find(~infLoc),ones(sum(~infLoc),1),G,size(X,1),size(X,2));
end