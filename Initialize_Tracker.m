function [Tracking] = Initialize_Tracker(data,tagged_data,Params)

med_filt = vision.MedianFilter;

if strcmp(data.Type,'uint8')
    dens_x = 0:(2^8-1);
    Tracking.maxgray = 2^8-1;
else
    dens_x = 0:(2^16-1);
    Tracking.maxgray = 2^16-1;
end
DensCellPoints = [];
DensBGPoints = [];
Tracking.B  = 0;

for t = 1:tagged_data.Frame_Num
    %% Pre-process Image
    % Read Image
    I = double(imread(data.Frame_name{t}));
    L = double(imread(tagged_data.Frame_name{t}));
    Tracking.L = L;
    % Calculate and remove Background lighting
    [I,B] =  CalcBGLighting(I,L==0);
    
    I = max(I,0);
    LCells = L>0;
    LBG = ~LCells;
    DensCellPoints = cat(1,DensCellPoints,I(LCells&I<Tracking.maxgray));
    DensBGPoints = cat(1,DensBGPoints,I(LBG));
    %% Create Kalman and predict
    if t==1
        states = Calculate_State(I,L);
        Kalmans = struct('ID',num2cell([states.ID]),'count',num2cell(1:size(states,2)),'enabled',true,'kalman',[],'num_props',0,...
            'prev_state',[],'children',[],'Contour',[],...
            'size',[],'HD',0.1,'new',true,'z_pred',[],'state_err',[],'U',[],'Mother',[],'Children',[],'cycle',0);
        Kalmans = arrayfun(@(s,k) fill_Kalman(s,k),states,Kalmans);
        Tracking.B = Tracking.B +B;
        Tracking.maxCellID = max(L(:));
        continue;
    end
    %% Predict and Correct
    prev_states = states;
    states =Calculate_State(I,L);
    Tracking.maxCellID = max(max(L(:)),Tracking.maxCellID);
    
    for n = 1:length(states)
        m = find([prev_states.ID]==states(n).ID);
        if isempty(m)
            mcount = max([Kalmans(:).count]);
            Kalmans(end+1).count = mcount+1;
            Kalmans(end) = fill_Kalman(states(n),Kalmans(end));
            Kalmans(end).ID = states(n).ID;
            Kalmans(end).HD = 0.1;
            Kalmans(end).num_props = 0;
            Kalmans(end).new = true;
            Kalmans(end).enabled = true;
            Kalmans(end).cycle = 0;
            Kalmans(end).weightedSize = states(n).weightedSize;
            continue;
        end
        states(n).kalman_state(6:10) =states(n).kalman_state(1:5)-prev_states(m).kalman_state(1:5);
        [z_pred, x_pred, P_pred] = predict(Kalmans(m).kalman);
        [z_corr,x_corr,P_corr] = correct(Kalmans(m).kalman,states(n).kalman_state);
        
        [p1y,p1x] = find(fullSingle(Kalmans(m).Contour)); [p2y,p2x]=find(fullSingle(states(n).Contour));
        p2x = p2x - states(n).kalman_state(1)+Kalmans(m).prev_state(1) ;
        p2y = p2y - states(n).kalman_state(2)+Kalmans(m).prev_state(2) ;
        p1 = [p1y,p1x];p2=[p2y,p2x];
        hd = HausdorffDist(p1,p2);
        if isempty(hd)
            hd = 0;
        end
        
        Kalmans(m).HD = Kalmans(m).HD*(t-2)/(t-1)+hd/(t-1) ;
        Kalmans(m).prev_state = z_corr';
        Kalmans(m).Contour = states(n).Contour;
        Kalmans(m).BW = states(n).BW;
        Kalmans(m).cycle = Kalmans(m).cycle+1;
        Kalmans(m).state = states(n).kalman_state;
    end
    Tracking.B = Tracking.B +B;
    
end

if isfield(Params.parameters,'useGMM')&&Params.parameters.useGMM
   
    Kbg = Params.parameters.Kbg;
    Kfg = Params.parameters.Kfg;
    gmmOpts = statset('MaxIter',500);
    gmmBG = fitgmdist(DensBGPoints(DensBGPoints>0),Kbg,'Options',gmmOpts);
    gmmFG = fitgmdist(DensCellPoints(DensCellPoints>0),Kfg,'Options',gmmOpts);
    gmmBG = gmdistribution(gmmBG.mu,gmmBG.Sigma,ones(1,Kbg)./Kbg);
    gmmFG = gmdistribution(gmmFG.mu,gmmFG.Sigma,ones(1,Kfg)./Kfg);
    Tracking.dens_BG = pdf(gmmBG,dens_x');
    Tracking.dens_cells = pdf(gmmFG,dens_x');
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
                gmmBG = fitgmdist(Kalmans(n).DensBGPoints(Kalmans(n).DensBGPoints>0),Kbg,'Options',gmmOpts);
                gmmFG = fitgmdist(Kalmans(n).DensCellPoints(Kalmans(n).DensCellPoints>0),Kfg,'Options',gmmOpts);
                gmmBG = gmdistribution(gmmBG.mu,gmmBG.Sigma,ones(1,Kbg)./Kbg);
                gmmFG = gmdistribution(gmmFG.mu,gmmFG.Sigma,ones(1,Kfg)./Kfg);
                Kalmans(n).dens_BG = pdf(gmmBG,dens_x');
                Kalmans(n).dens_cells = pdf(gmmFG,dens_x');
            end
        end
    end
else
    u = (4/(3*min(numel(DensBGPoints)+numel(DensCellPoints))))^(1./5)*std(cat(1,DensCellPoints,DensBGPoints));
    dens_cells = FastKDE(DensCellPoints,dens_x,u);
    dens_BG = FastKDE(DensBGPoints,dens_x,u);
    Tracking.dens_cells = dens_cells;
    Tracking.dens_BG = dens_BG;
end
Tracking.priorBG = sum(L(:)==0)./length(L(:));
Tracking.priorCell = 1-Tracking.priorBG;
Tracking.dens_x = dens_x;
Tracking.B = Tracking.B./t;

Tracking.Kalmans = Kalmans;
Tracking.current_t = t+1;
Tracking.med_filt = med_filt;


