function [y] = correct_wei(x, p_equ, delta)
% compute the corrected weight for the weighted preconditioner
% delta + abs(x).^(p-2) or (delta + abs(x)).^(p-2)
if p_equ == 2
    y = ones(size(x));
elseif p_equ > 2
    y = delta + (abs(x).^(p_equ-2));
else % p_equ < 2
    y = (delta + abs(x)).^(p_equ-2);
end
end