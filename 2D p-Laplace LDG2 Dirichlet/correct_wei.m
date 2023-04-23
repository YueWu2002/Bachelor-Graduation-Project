function [y] = correct_wei(x, p_equ, delta)
% compute the corrected weight for the weighted preconditioner
% delta + abs(x).^(p-2) or (delta + abs(x)).^(p-2)
if p_equ == 2.0
    y = ones(size(x));
elseif p_equ > 2.0
    y = delta + (abs(x).^(p_equ-2));
elseif p_equ > 1.0 && p_equ < 2.0
    y = (delta + abs(x)).^(p_equ-2);
else
    error('correct_wei: error p_equ!');
end
end