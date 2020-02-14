function FA = range_mod(A,d)
% A: input set, given as list of indices

if isempty(A)
    FA = 0;
else
    FA = max(d-1,0) + max(A)-min(A)+1;
end

end