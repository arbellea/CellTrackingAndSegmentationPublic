function [OutStruct] = read_parameters_txt(fname)

fprintf('Reading parameters file... Please wait...\n');



OutStruct = struct;

fid = fopen(fname);  

while ~feof(fid)
    tline = fgetl(fid);
    tline = strtrim(tline);
    if isempty(tline)|| tline(1) == '%'
        continue;
    end
    splt_line = strsplit(tline);
    if splt_line{1} == '*'
        fldname = replace_wspace(splt_line,2);
    else
        splt_line = strsplit(tline,'=');
        splt_line = strtrim(splt_line);
        sub_splt = strsplit(splt_line{1});
        subfldname = replace_wspace(sub_splt,1);
        val = eval(splt_line{2});
        OutStruct.(fldname).(subfldname) = val;
    end    
end

fprintf('Done Reading parameters file\n');


        
        
        
        
