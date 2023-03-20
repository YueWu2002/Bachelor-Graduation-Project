uu = @(x,y) exp(sin(pi*x))*cos(pi*(x+y));
qqx = @(x,y) exp(sin(pi*x))*pi*(cos(pi*x)*cos(pi*(x+y))-sin(pi*(x+y)));
qqy = @(x,y) -exp(sin(pi*x))*pi*sin(pi*(x+y));
ff = @(x,y) exp(sin(pi*x))*pi*pi*(cos(pi*(x+y))*(2.0-cos(pi*x)^2+sin(pi*x))+2*cos(pi*x)*sin(pi*(x+y)));
g_D = uu;
g_N = @(x,y) [qqx(x,y), qqy(x,y)];

show_spy = false;
p = 2;      % advise: 1<=p<=7
vel = [1.0, 0];

mymesh1
mesh{1} = read_mesh(msh, false);
mymesh2
mesh{2} = read_mesh(msh, false);
mymesh3
mesh{3} = read_mesh(msh, false);
mymesh4
mesh{4} = read_mesh(msh, false);
mymesh5
mesh{5} = read_mesh(msh, false);
mymesh6
mesh{6} = read_mesh(msh, false);

m_test = 6;
err = zeros(m_test, 2);
for k = 1: m_test
    sol = poisson2D_MDLDG(mesh{k}, ff, g_D, g_N, p, vel, show_spy, 'reduce');
    err(k,1:2) = compute_error(mesh{k}, sol, p, uu, qqx, qqy);
end
% Remark: 使用Schur补对精度都暂时观察不到影响，对计算速度可以提高。
% 问题：几乎所有的时间都用在组装矩阵上了。(通过我重写程序中的边界积分，时间迅速下降了，例如k=6的最密网格，时间从640秒缩到7.5秒)；使用 Schur 补会提高数值精度，且可以加速计算(k=7最密网格，使用Schur补的直接解法需要4.3秒，原本解法需要28.6秒)

ord = nan(m_test, 2);
ord(2:end, :) = log2(err(1:end-1,:)./err(2:end,:));


format short e
disp(err)

format short
disp(ord)
