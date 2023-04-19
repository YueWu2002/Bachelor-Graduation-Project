function [y] = correct_fullwei(x, p_equ, delta)
% compute the corrected weight for the full-weighted preconditioner
% delta + abs(x).^(p-4) or (delta + abs(x)).^(p-4)
if p_equ == 4
    y = ones(size(x));
elseif p_equ > 4
    y = delta + (abs(x).^(p_equ-4));
else % p_equ < 4
    y = (delta + abs(x)).^(p_equ-4);
end
end