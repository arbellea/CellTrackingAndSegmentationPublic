function A = fullSingle(S)
switch S.class
    case 1 
        c = 'logical';
    case 2 
        c = 'single';
    case 3 
        c = 'double';
        
    case 4 
        c = 'int8';
        
    case 5 
        c = 'int16';
        
    case 6 
        c = 'int32';
        
    case 7
        c = 'int64';
        
    case 8
        c = 'uint8';
        
    case 9
        c = 'uint16';
        
    case 10
        c = 'uint32';
        
    case 11
        c = 'uint64';
    otherwise
        error('Cannot create create full matrix');
end
if S.class==1
    A = false(S.size);
    A(S.idx) = logical(S.val);
else
    A = zeros(S.size,c);
    A(S.idx) = cast(S.val,c);

end
end
