function [fldname] = replace_wspace(strcell,start)
fldname = strcell{start};
for i = (start+1):length(strcell)
    if ~isempty(strcell{i})
        fldname = strcat(fldname,'_',strcell{i});
    end
end

    