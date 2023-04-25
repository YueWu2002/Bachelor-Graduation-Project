p_equ = 1.5;
q_equ = p_equ / (p_equ-1.0);

cd exact
% exact1 % only for p_equ == 2
% sigma = 0;
% exact2
% a = 0.3;
% exact3
exact4
cd ..

g_D = uu;
g_N = @(x,y) [fluxx(x,y), fluxy(x,y)];

% first, we test our program with Poisson's equation
% then, we generalize to p-Laplace

%{
cd("mesh");
mymesh1
raw_mesh{1} = msh;
mymesh2
raw_mesh{2} = msh;
mymesh3
raw_mesh{3} = msh;
mymesh4
raw_mesh{4} = msh;
mymesh5
raw_mesh{5} = msh;
mymesh6
raw_mesh{6} = msh;
cd ..
%}
cd("newmesh");
mymesh1
raw_mesh{1} = msh;
mymesh2
raw_mesh{2} = msh;
mymesh3
raw_mesh{3} = msh;
mymesh4
raw_mesh{4} = msh;
mymesh5
raw_mesh{5} = msh;
mymesh6
raw_mesh{6} = msh;
cd ..

for k = 1: 6
    mesh{k} = read_mesh(raw_mesh{k}, false);
end

show_spy = false;
p = 6;      % advise: 1<=p<=6
vel = [1.0, 0]; % auxiliary H(div) function
basis_type = 'Bernstein'; % 'Taylor' or 'Bernstein' (very recommended due to floating-point stability)

m_test = 6;
err = zeros(m_test, 3);
for k = 1: m_test
    sol = pLaplace2D_LDG(mesh{k}, ff, g_D, g_N, p, p_equ, vel, show_spy, basis_type, true);
    % remark: please DO NOT use parfor because it is significantly slower than vectorization.

    err(k,:) = compute_error(mesh{k}, sol, p, p_equ, uu, qqx, qqy, fluxx, fluxy, basis_type);

    print(1, '-depsc2', '-vector', ['Jh_', num2str(k)], '-r300');
    print(2, '-depsc2', '-vector', ['rho_', num2str(k)], '-r300');
    print(3, '-depsc2', '-vector', ['gradnorm_', num2str(k)], '-r300');
end

ord = nan(m_test, 3);
ord(2:end, :) = log2(err(1:end-1,:)./err(2:end,:));

%{
format short e
disp(err)

format short
disp(ord)
%}

for i = 1: size(err, 1)
    fprintf('%d ', i);
    for j = 1: size(err, 2)
        fprintf("& %.4e & %.4f ", [err(i,j), ord(i,j)]);
    end
    fprintf('\\\\\n');
end