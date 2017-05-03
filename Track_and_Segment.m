function [save_dir_name,Kalmans] = Track_and_Segment(data,Tracking,Params,varargin)
rng(1);
if ~isempty(varargin)&&isstruct(varargin{1})
    taggedData = varargin{1};
end

Link = struct();
Save_images = Params.Flags.WriteVideo;
SaveCheckPoints = Params.Flags.SaveCheckPoints;
LoadCheckPoints = Params.Flags.LoadCheckPoints;
segParams = Params.parameters;
save_debug = Params.Flags.SaveDebug;
deleteIfErr = Params.Flags.deleteIfErr;
ISBI = Params.Flags.ISBI;
t = datestr(now);
t(t==' ')='_';
t(t==':')='-';
if SaveCheckPoints||LoadCheckPoints
    t = 'CheckPoints';
end
if Save_images
    save_dir_name =fullfile(getenv('HOME'),'Outputs',sprintf('Results_%s_%s',Params.General.Name,t));
    save_dir_vis = fullfile(save_dir_name,'Visualize');
    save_dir_res = fullfile(save_dir_name,'Results');
    mkdir(save_dir_vis);
    mkdir(save_dir_res);
    mkdir(save_dir_name);
    if isfield(Params.Flags,'saveCode')&&Params.Flags.saveCode
        ticSaveCode = tic;
    zipPath  = fullfile(save_dir_name,'Code.zip');
    mainpath = fileparts(which('main.m'));
    filelist = dir(mainpath);
    hidden = arrayfun(@(f) (f.name(1)~='.')&&all(~strcmpi(f.name,{'mitodix.log','license'})),filelist);
    fnames = arrayfun(@(f) fullfile(mainpath,f.name),filelist(hidden),'uniformoutput',false);
    zip(zipPath,fnames);
    tocSaveCode = toc(ticSaveCode);
    fprintf('Done Saving Code in %0.3f seconds...\n',tocSaveCode)
    end
end
save_dir_checkp = fullfile(save_dir_name,'CheckPoints');
if SaveCheckPoints
    
    mkdir(save_dir_checkp);
end
if any(save_debug)
    mkdir(fullfile(save_dir_name,'Debug'))
end
min_cell_size = 1;
Height = data.Height+20;
Width = data.Width+20;
Kalmans = Tracking.Kalmans;
stopFrame = Tracking.stop_frame;
try  
    if LoadCheckPoints&& exist(save_dir_checkp,'dir')
        file_list = dir(fullfile(save_dir_checkp,'*.mat'));
        fileNames = {file_list.name};
        expr = 'CheckPoint_(\d+).mat';
        tokens = cellfun(@(str) regexp(str,expr,'tokens'),fileNames,'uniformoutput',false);
        emptyTokens = cellfun(@isempty, tokens);
        tokens = tokens(~emptyTokens);
        token = cellfun(@(token) str2double(token{1}{1}),tokens);
        if ~isempty(token)
            if islogical(LoadCheckPoints)
                t = max(token);
            elseif isscalar(LoadCheckPoints)
                tokensSort = sort(token);
                t =tokensSort(find((tokensSort-LoadCheckPoints)<=0,1,'last'));
                if isempty(t)
                    t = tokensSort(1);
                end
            end
            stopFrame = Tracking.stop_frame;
            Tracking = load(fullfile(save_dir_checkp,sprintf('CheckPoint_%d.mat',t)));
            Kalmans = Tracking.Kalmans;
            if isfield(Tracking,'Link')
                Link = Tracking.Link;
            end
            Tracking.current_t = t;
        end
        
    end
    for t =Tracking.current_t:stopFrame
        if SaveCheckPoints>0&& mod(t,SaveCheckPoints)==0 
            Tracking.Link = Link;
            save(fullfile(save_dir_name,'CheckPoints',sprintf('CheckPoint_%d.mat',t)),'-struct','Tracking');
        end
        tStartFrame = tic;
        [~,mbgidx] = max(Tracking.dens_BG);
        mBG1 = Tracking.dens_x(mbgidx);
        I = double(imread(data.Frame_name{t}));
        I = I-Tracking.B;
        I = max(I,0);
        I = step(Tracking.med_filt,I);
        %figure(1); imshow(I,[]);
        I_prev = double(imread(data.Frame_name{t-1}));
        I_prev = I_prev-Tracking.B;
        I_prev = max(I_prev,0);
        I_prev = step(Tracking.med_filt,I_prev);
        
        tSeg = tic;
        fprintf('Start Segmentation of frame %d...\n',t);
        [L,L_New_Cells,Kalmans,z_pred,z_pred_orig,cog_diff,Debug] = Fuzzy_Segmentation(Tracking,Kalmans,I,I_prev,segParams,any(save_debug));
        disabeledKalmans = Kalmans(~[Kalmans.enabled]);
        z_pred_orig = z_pred_orig(~[Kalmans.enabled]);
        L_New_Cells_orig = L_New_Cells;
        timeSeg = toc(tSeg);
        fprintf('Done Segmentation of frame %d in %f seconds...\n',t,timeSeg);
        
        Obj_num = length(Kalmans);
        disabeledKalmans = Kalmans(~[Kalmans.enabled]);
        Kalmans = Kalmans([Kalmans.enabled]);
        tCalc = tic;
        states = Calculate_State(I,L,Kalmans);
        timeCalc = toc(tCalc);
        fprintf('Done Calculate State of frame %d in %0.3f seconds...\n',t,tCalc);
        tUp = tic;
        for n = 1:length(Kalmans)
            
            if isempty(states(n).kalman_state)
                Kalmans(n).state = Kalmans(n).z_pred';
                Kalmans(n).state(6:end)=0;
                Kalmans(n).HD = 10;
                Kalmans(n).num_props=endiKalmans(n).num_props+1;
                continue
            end
            %states(n).kalman_state(1:2) = states(n).kalman_state(1:2)-cog_diff;
            Kalmans(n).state =states(n).kalman_state';
            Kalmans(n).state_err(end+1,:) = Kalmans(n).z_pred-states(n).kalman_state;
            Kalmans(n).Contour = states(n).Contour;
            Kalmans(n).BW = states(n).BW;
            Kalmans(n).HD = states(n).HD;
            Kalmans(n).weightedSize = states(n).weightedSize;
            
            p=1;
            Kalmans(n).state(6:end) = p*(Kalmans(n).state(1:5)-Kalmans(n).prev_state(1:5))+(1-p)*(Kalmans(n).prev_state(6:end));
            Kalmans(n).state(6:7) = Kalmans(n).state(6:7)-cog_diff';
            if Kalmans(n).num_props
                Kalmans(n).state(6:end) = 0;
            end
            Kalmans(n).num_props = 0;
            Kalmans(n).size = states(n).size;
            Kalmans(n).cycle = Kalmans(n).cycle+1;
            Tracking.ISBI_RES(Tracking.ISBI_RES(:,1)==Kalmans(n).ID,3) = Tracking.ISBI_RES(Tracking.ISBI_RES(:,1)==Kalmans(n).ID, 3)+1;
            
        end
        timeUp = toc(tUp);
        fprintf('Done State update of frame %d in %f seconds...\n',t,timeUp);
        
        if ~isfield(Kalmans,'state')
            Kalmans(1).state =[];
        end
        
        if ~isempty(L_New_Cells)&&any(L_New_Cells(:)>0)
            tnCalc = tic;
            states =Calculate_State(I,L_New_Cells);
            timeNCalc = toc(tnCalc);
            fprintf('Done New State Calculation of frame %d in %f seconds...\n',t,timeNCalc);
            tnUp = tic;
            uniqe_L = unique(L_New_Cells(L_New_Cells>0));
            m =length(Kalmans);
            for n = 1:size(states,2)
                mu = states(n).kalman_state(1:2);
                if mu(1)<=1||mu(2)<=1||mu(1)>(Width)||mu(2)>(Height)
                    continue;
                end
                for l = numel(Link):-1:1
                    
                    if ~isfield(Link(1),'Mother')||Link(l).Time<t
                        break;
                    end
                    Link(l).Children(Link(l).Children==n) =  Tracking.maxCellID +1;
                    
                end
                m = m+1;
                Kalmans(m).new = true;
                Kalmans(m).count =m;
                Kalmans(m).ID =Tracking.maxCellID + uniqe_L(n);
                
                Kalmans(m).enabled = true;
                Kalmans(m).kalman = Create_New_Kalman(states(n).kalman_state,8^2*eye(10),2*eye(10));
                Kalmans(m).num_props = 0;
                Kalmans(m).Children =[];
                Kalmans(m).prev_state = states(n).kalman_state;
                Kalmans(m).weightedSize = inf;
                Kalmans(m).state = states(n).kalman_state;
                Kalmans(m).HD =100;
                Kalmans(m).size = states(n).size;
                Kalmans(m).z_pred = zeros(1,10);
                Kalmans(m).state_err = [];
                Kalmans(m).BW = states(n).BW;
                Kalmans(m).U = (states(n).BW);
                Kalmans(m).Contour = states(n).Contour;
                Kalmans(m).Mother = [];
                Kalmans(m).cycle = 1;
                
                L(fullSingle(Kalmans(m).BW)) = Kalmans(m).ID;
                Tracking.ISBI_RES = cat(1,Tracking.ISBI_RES,[Kalmans(m).ID,t-1,t-1,0]);
            end
            Tracking.maxCellID  = Tracking.maxCellID + uniqe_L(n);
            timeNUp = toc(tnUp);
            fprintf('Done New State Update of frame %d in %f seconds...\n',t,timeNUp);
        end
    
        if ~isempty(disabeledKalmans)
            tmitLink = tic;
            for n = 1:length(disabeledKalmans)
                
                if ~isempty(disabeledKalmans(n).Children)
                    motherID = disabeledKalmans(n).ID;
                    children = zeros(length(disabeledKalmans(n).Children),1);
                    for child = 1:length(disabeledKalmans(n).Children)

                        Kalmans([Kalmans.ID]==L(disabeledKalmans(n).Children(child).PixelIdxList(1))).Mother=disabeledKalmans(n).ID;
                        children(child) =L(disabeledKalmans(n).Children(child).PixelIdxList(1)) ;
                        Tracking.ISBI_RES(Tracking.ISBI_RES(:,1)==children(child),4)=motherID;
                    end
                    if ~isfield(Link,'Mother')
                        Link(1).Mother = motherID;
                        Link(1).Children = children;
                        Link(1).Time = t;
                    else
                        Link(end+1).Mother = motherID;
                        Link(end).Children = children;
                        Link(end).Time = t;
                    end
                    
                end
                
                Tracking.ISBI_RES(Tracking.ISBI_RES(:,1)==disabeledKalmans(n).ID,3)=t-2;
            end
            timeMitLink = toc(tmitLink);
            fprintf('Done Mitosis Link of frame %d in %f seconds...\n',t,timeMitLink);
        end
    
        tKUp = tic;
        for n = 1:length(Kalmans)
            if  Kalmans(n).enabled
                [z_corr(:,n),~,~] = correct(Kalmans(n).kalman,Kalmans(n).state);
                Kalmans(n).prev_state = z_corr(:,n);
                mes_mus(n,:)=Kalmans(n).state(1:2)';%.*[Width,Height];
                if mes_mus(n,1)<=1||mes_mus(n,2)<=1||mes_mus(n,1)>=(Width)||mes_mus(n,2)>=(Height)||Kalmans(n).size<min_cell_size % mes_mus(n,1)<=10||mes_mus(n,2)<=10||mes_mus(n,1)>=(Width-10)||mes_mus(n,2)>=(Height-10)||Kalmans(n).size<min_cell_size
                    Kalmans(n).enabled = false;
                    
                end
            end
            
            
        end
        z_mat= cell2mat(z_pred')';
        pred_mus = z_mat(1:2,:)';
        cor_mus = z_corr(1:2,:)';
        timeKUp = toc(tKUp);
        fprintf('Done Kalman Update of frame %d in %f seconds...\n',t,timeKUp);
        [Kalmans([Kalmans.new]).new] = deal(false);
        tSave = tic;
        if Save_images
            
            if ISBI
                frame_name = fullfile(save_dir_res,sprintf('mask%03d.tif',t-1));
            else
                [~,fname,ext]=fileparts(data.Frame_name{t});
                frame_name = sprintf('Seg_%s%s',fname,ext);
                frame_name = fullfile(save_dir_res,frame_name);
            end
            imwrite(uint16(L),frame_name);
            
            
        end
        timeSave = toc(tSave);
        
        fprintf('Done Save frame %d in %f seconds...\n',t,timeSave);
        Tracking.L = L;
       
        
        if ismember(t,save_debug)
            Kalmans_debug= Kalmans; %#ok<*NASGU>
            Seg_Debug = Debug;
            save(fullfile(save_dir_name,'Debug',sprintf('Kalmans_debug_frame_%d.mat',t)),'Kalmans_debug','Seg_Debug','Tracking','-v7.3');
        end
       
        tKDE = tic;
        if ~exist('DensCellPoints','var')&&~exist('DensBGPoints','var')
        DensCellPoints =[];
        DensBGPoints = [];
        DensEdgePoints = [];
        end
         LCells = L>0;
         
         LBG = ~(LCells);
            
        BG_est_refresh = Params.parameters.BG_est_refresh;
        if mod(t,BG_est_refresh)==0&&exist('DensCellPoints','var')&&exist('DensBGPoints','var')
            LCells = L>0;
            
            LBG = ~(LCells);
            
            
            DensCellPoints = cat(1,DensCellPoints,I(LCells&I<Tracking.maxgray));
            DensBGPoints = cat(1,DensBGPoints,I(LBG));
            
           
            if isfield(Params.parameters,'useGMM')&&Params.parameters.useGMM
                Kbg = Params.parameters.Kbg;
                Kfg = Params.parameters.Kfg;
                gmmOpts = statset('MaxIter',500);
                gmmBG = fitgmdist(DensBGPoints(DensBGPoints>0),Kbg,'Options',gmmOpts);
                gmmFG = fitgmdist(DensCellPoints(DensCellPoints>0),Kfg,'Options',gmmOpts);
                gmmBG = gmdistribution(gmmBG.mu,gmmBG.Sigma,ones(1,Kbg)./Kbg);
                gmmFG = gmdistribution(gmmFG.mu,gmmFG.Sigma,ones(1,Kfg)./Kfg);
                Tracking.dens_BG = pdf(gmmBG,Tracking.dens_x');
                Tracking.dens_cells = pdf(gmmFG,Tracking.dens_x');
                if isfield(Params.parameters,'useLocalGL')&&Params.parameters.useLocalGL
                    
                    
                    Kbg = Params.parameters.Kbg;
                    Kfg = Params.parameters.Kfg;
                    gmmOpts = statset('MaxIter',500);
                    for n = 1:numel(Kalmans)
                        if Kalmans(n).enabled
                            if ~isfield(Kalmans(n),'DensCellPoints')
                                Kalmans(n).DensCellPoints =[];
                                Kalmans(n).DensBGPoints =[];
                            end
                            cent = Kalmans(n).state(2:-1:1);
                            L_cropped = CropImage(L,cent,Params.parameters.patchSize,Params.parameters.patchSize);
                            I_cropped = CropImage(I,cent,Params.parameters.patchSize,Params.parameters.patchSize);
                            Kalmans(n).DensCellPoints = cat(1,Kalmans(n).DensCellPoints,I_cropped(L_cropped(:)>0));
                            Kalmans(n).DensBGPoints = cat(1,Kalmans(n).DensBGPoints,I_cropped(L_cropped(:)==0));
                            if isempty(Kalmans(n).DensBGPoints(Kalmans(n).DensBGPoints>0))||isempty(Kalmans(n).DensCellPoints(Kalmans(n).DensCellPoints>0))
                             Kalmans(n).dens_BG  = Tracking.dens_BG;
                             Kalmans(n).dens_cells = Tracking.dens_cells;
                            else
                            gmmBG = fitgmdist(Kalmans(n).DensBGPoints(Kalmans(n).DensBGPoints>0),Kbg,'Options',gmmOpts);
                            gmmFG = fitgmdist(Kalmans(n).DensCellPoints(Kalmans(n).DensCellPoints>0),Kfg,'Options',gmmOpts);
                            gmmBG = gmdistribution(gmmBG.mu,gmmBG.Sigma,ones(1,Kbg)./Kbg);
                            gmmFG = gmdistribution(gmmFG.mu,gmmFG.Sigma,ones(1,Kfg)./Kfg);
                            Kalmans(n).dens_BG = pdf(gmmBG,Tracking.dens_x');
                            Kalmans(n).dens_cells = pdf(gmmFG,Tracking.dens_x');
                            end
                            Kalmans(n).DensCellPoints =[];
                            Kalmans(n).DensBGPoints =[];
                        end
                    end
                end
                
            else
                u = (4/(3*min(numel(DensBGPoints)+numel(DensCellPoints))))^(1./5)*max(std(DensCellPoints),std(DensBGPoints));
                dens_cells = FastKDE(DensCellPoints,Tracking.dens_x,u);
                dens_BG = FastKDE(DensBGPoints,Tracking.dens_x,u);
                zz = find(dens_cells==0&dens_cells==0);
                zd = knnsearch([DensCellPoints;DensBGPoints],zz');
                dens_cells(zz(zd<=numel(DensCellPoints))) = eps;
                dens_BG(zz(zd>numel(DensCellPoints))) = eps;
                mucells = mean(DensCellPoints);
                mubg = mean(DensBGPoints);
                Tracking.dens_cells = dens_cells;
                Tracking.dens_BG = dens_BG;
                Tracking.priorBG = sum(L(:)==0)./length(L(:));
                Tracking.priorCell = 1-Tracking.priorBG;
                
            end
        else
            if isfield(Params.parameters,'useLocalGL')&&Params.parameters.useLocalGL
                for n = 1:numel(Kalmans)
                if Kalmans(n).enabled
                    if ~isfield(Kalmans(n),'DensCellPoints')
                        Kalmans(n).DensCellPoints =[];
                        Kalmans(n).DensBGPoints =[];
                    end
                    cent = Kalmans(n).state(2:-1:1);
                    L_cropped = CropImage(L,cent,Params.parameters.patchSize,Params.parameters.patchSize);
                    I_cropped = CropImage(I,cent,Params.parameters.patchSize,Params.parameters.patchSize);
                    Kalmans(n).DensCellPoints = cat(1,Kalmans(n).DensCellPoints,I_cropped(L_cropped(:)>0));
                    Kalmans(n).DensBGPoints = cat(1,Kalmans(n).DensBGPoints,I_cropped(L_cropped(:)==0));
                end
                end
            end
                DensCellPoints = cat(1,DensCellPoints,I(LCells&I<Tracking.maxgray));
                DensBGPoints = cat(1,DensBGPoints,I(LBG));
            
        end
        timeKDE = toc(tKDE);
        
        fprintf('Done KDE Update %d in %f seconds...\n',t,timeKDE);
        tEndFrame = toc(tStartFrame);
        Kalmans = Kalmans([Kalmans.enabled]);
        Tracking.Kalmans = Kalmans;
      
        whos Kalmans
        fprintf('Elapsed time for Frame %d is %2.5f\n',t,tEndFrame)
        
    end
    if ISBI
        dlmwrite(fullfile(save_dir_name,'res_track.txt'),Tracking.ISBI_RES,' ')
    end
    if Save_images
        save(fullfile(save_dir_name,'Link.mat'),'Link');
    end
catch err
    if SaveCheckPoints
            Tracking.Link = Link;
            save(fullfile(save_dir_name,'CheckPoints',sprintf('CheckPoint_%d.mat',t)),'-struct','Tracking');
    end
        
    if ismember(t,save_debug)
        mkdir([save_dir_name,'_Debug']);
        save(fullfile([save_dir_name,'_Debug'],sprintf('Kalmans_debug_frame_%d.mat',t)),'-v7.3');
    end
    if deleteIfErr&&isdir(save_dir_name)
        try
            rmdir(save_dir_name,'s');
        catch err1
            disp('Could not delete output directory')
        end
        
    end
    rethrow(err);
end

end
