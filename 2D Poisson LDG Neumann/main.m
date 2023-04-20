uu = @(x,y) exp(sin(pi*x))*cos(pi*(x+y));
qqx = @(x,y) exp(sin(pi*x))*pi*(cos(pi*x)*cos(pi*(x+y))-sin(pi*(x+y)));
qqy = @(x,y) -exp(sin(pi*x))*pi*sin(pi*(x+y));
ff = @(x,y) exp(sin(pi*x))*pi*pi*(cos(pi*(x+y))*(2.0-cos(pi*x)^2+sin(pi*x))+2*cos(pi*x)*sin(pi*(x+y)));
g_D = uu;
g_N = @(x,y) [qqx(x,y), qqy(x,y)];
c0 = -0.454719441676668752704457919840;

show_spy = false;
p = 1;      % advise: 1<=p<=9
vel = [1.0, 0];
basis_type = 'Bernstein'; % 'Taylor' or 'Bernstein' (very recommended due to floating-point stability)

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

for k = 1: 6
    mesh{k} = read_mesh(raw_mesh{k}, false);
end

m_test = 6;
err = zeros(m_test, 2);
for k = 1: m_test
    sol = poisson2D_MDLDG(mesh{k}, ff, g_D, g_N, c0, p, vel, show_spy, basis_type, 'direct');
    err(k,:) = compute_error(mesh{k}, sol, p, uu, qqx, qqy, basis_type);
end
% Remark: 使用Schur补对精度有一定影响(尤其密网格高阶数时)，对计算速度可以有效提高。
% 问题：几乎所有的时间都用在组装矩阵上了。(通过我重写程序中的边界积分，时间迅速下降了，例如k=6的最密网格，时间从640秒缩到7.5秒)；使用 Schur 补会提高数值精度，且可以加速计算(k=7最密网格，使用Schur补的直接解法需要4.3秒，原本解法需要28.6秒)

ord = nan(m_test, 2);
ord(2:end, :) = log2(err(1:end-1,:)./err(2:end,:));


format short e
disp(err)

format short
disp(ord)


for i = 1: size(err, 1)
    fprintf('%d ', i);
    for j = 1: size(err, 2)
        fprintf("& %.4e & %.4f ", [err(i,j), ord(i,j)]);
    end
    fprintf('\\\\\n');
end