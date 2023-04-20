function [val] = modified_sign(a)
if a >= 0.0
    val = 1;
else
    val = -1;
end
end