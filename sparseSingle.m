function S = sparseSingle(A)
switch class(A)
    case 'logical'
        c = 1;
    case 'single'
        c = 2;
    case 'double'
        c = 3;
    case 'int8'
        c = 4;
    case 'int16'
        c = 5;
    case 'int32'
        c = 6;
    case 'int64'
        c = 7;
    case 'uint8'
        c = 8;
    case 'uint16'
        c = 9;
    case 'uint32'
        c = 10;
    case 'uint64'
        c = 11;
    otherwise
        error('Cannot create create sparse matrix for type %s',class(A));
end
S.class = uint8(c);
s = size(A);
if all(s<intmax('uint8'))
    S.size = uint8(size(A));
elseif all(s<intmax('uint16'))
    S.size = uint16(size(A));
elseif all(s<intmax('uint32'))
    S.size = uint32(size(A));
else
    S.size = uint64(size(A));
end
idx = find(A>0);
midx = max(idx);
if midx<intmax('uint8')
    S.idx = uint8(idx);
elseif midx<intmax('uint16')
    S.idx = uint16(idx);
elseif midx<intmax('uint32')
    S.idx = uint32(idx);
else
    S.idx = uint64(idx);
end
S.val = single(full(A(idx)));
end
