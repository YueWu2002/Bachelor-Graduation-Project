function [val] = modified_sign(a)
if a >= 0.0
    val = 1.0;
else
    val = -1.0;
end
end