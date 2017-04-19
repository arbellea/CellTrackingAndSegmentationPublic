function CreateOutputVideo(inputInfo,saveDir)
%{
filelist = dir(fullfile(segDir,ext));
base = ['Seg_%d',ext(2:end)];
t = arrayfun(@(s) sscanf(s.name,base,1),filelist);
[~,sort_idx] = sort(t);
sort_list = filelist(sort_idx);
%}
%segInfo = Load_Data(segDir,segExt);
%%
writerObj = VideoWriter(fullfile(saveDir,'Results.avi'));
writerObj.FrameRate=2;
open(writerObj);
c = true;
%%
load(fullfile(saveDir,'Link.mat'));
longName = {};
Base = [];
Gen = {};
Num = 0;

for i = 1:length(Link)
    if ~isfield(Link(i),'Mother')||isempty(Link(i).Mother)
        continue;
    end
    
    if length(longName)<Link(i).Mother||isempty(longName{Link(i).Mother})
       longName{Link(i).Mother} = sprintf('%d',Link(i).Mother);
       Base(Link(i).Mother) = Link(i).Mother;
       Gen{Link(i).Mother} = char(64);
       Num(Link(i).Mother) = 0;
    end
    for n = 1:length(Link(i).Children)
        %longName{Link(i).Children(n)} = sprintf('%s.%d',longName{Link(i).Mother},n);
        gen = char(Gen{Link(i).Mother}+1);
        maxNum = max(Num(Base==Base(Link(i).Mother)&cellfun(@(g) strcmpi(g,gen),Gen)));
        if isempty(maxNum)
            maxNum = 0;
        end
        
        longName{Link(i).Children(n)} = sprintf('%d.%s.%d',Base(Link(i).Mother),gen,maxNum+1);
        Gen{Link(i).Children(n)} = gen;
        Num(Link(i).Children(n)) = maxNum+1;
        Base(Link(i).Children(n)) = Base(Link(i).Mother);
    end
end
%%
try
    
    idprev = [];
    for t = 1:inputInfo.Frame_Num
        disp(t);
        I = imread(inputInfo.Frame_name{t});
        I_Raw = I;
        I = double(I);
        p = prctile(I(:),[0.1,99.9]);
        I = max(I,p(1)); I = min(I,p(2));
        I =(I-min(I(:)))/(max(I(:))-min(I(:)))*255;
        I = uint8(I);
        
        
        
        
        [~,fname,fext]=fileparts(inputInfo.Frame_name{t});
        segfile = fullfile(saveDir,'Results',sprintf('Seg_%s%s',fname,fext));
        
        if ~exist(segfile,'file') && c
            continue;
        elseif ~exist(segfile,'file') && ~c
            break;
        end
        c = false;
        Seg = imread(segfile); %Seg = Seg(11:end-10,11:end-10);
        reg = regionprops(Seg,I_Raw,'Centroid','Area','MeanIntensity');
        
        cents = [reg(:).Centroid];
        cents(isnan(cents))=[];
        cents = reshape(cents,2,[]);
        id = unique(Seg(Seg>0));
        Contours = arrayfun(@(l) find(bwperim(Seg==l)),id,'uniformoutput',0);
        mval =intmax(class(I));
        
        IR = uint8(I); IG = uint8(I); IB = uint8(I);
        C = cat(1,Contours{:});
        CNew = (cat(1,Contours{~ismember(id,idprev)}));
        C2d = unique(C);
        CNew2d = unique(CNew);
        IR(C2d(:))=mval;IG(C2d(:))=0;IB(C2d(:))=0;IB(CNew2d(:)) = mval;
        
        text_loc = [cents(1,:);cents(2,:)-5]';
        cell_num=num2str(id);
        cell_num_cell=cellstr(cell_num);
        for iid = 1:length(id) 
            if id(iid)>length(longName)
                continue;
            end
            if ~isempty(longName{id(iid)})
                cell_num_cell{iid} = longName{id(iid)};
            end
        end
        frame = insertText(IB,int32(text_loc),cell_num_cell','FontSize',20,'TextColor','white','BoxOpacity',0);
        frame = cat(3,IR,IG,frame(:,:,3));
        idprev = id;
        visFile = fullfile(saveDir,'Visualize',sprintf('Vis_%s%s',fname,fext));
        imwrite(frame,visFile);
        writeVideo(writerObj,frame);
    end
    close(writerObj);
catch err
    close(writerObj);
    rethrow(err);
end
end