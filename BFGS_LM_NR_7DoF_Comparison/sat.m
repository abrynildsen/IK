function q_sat = sat(q,qlims)
% outputs a value of q that is between the lower uand upper bounds
for j = 1:length(q)
    q_sat(j,1) = min(max(q(j),qlims(j,1)),qlims(j,2));
end

