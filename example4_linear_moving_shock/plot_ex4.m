function plot_ex4(K, n, limiter);
% para: n, P_n polynomial;

addpath('../src/');
format long;
f = @(x) 1.*(x >= 3.2 & x <= 4.2);

switch n;
case 120;
  numer1 = load(['example4_Nx120_K',num2str(K),'_PP',num2str(limiter),'.dat']);
  x1 = numer1(:,1); y1 = numer1(:,3);
  plot(x1, f(x1), '-k', x1, y1, 'ro');
  min(y1)
case 40;
  numer2 = load(['example4_Nx40_K',num2str(K),'.dat']);
  x2 = numer2(:,1); y2 = numer2(:,3);
  plot(x2, f(x2), '-k', x2, y2, 'ro');
  ex2 = f(x2);
case 80;
  numer3 = load(['example4_Nx80_K',num2str(K),'.dat']);
  x3 = numer3(:,1); y3 = numer3(:,3);
  plot(x3, f(x3), '-k', x3, y3, 'ro');
case 160;
  numer4 = load(['example4_Nx160_K',num2str(K),'.dat']);
  x4 = numer4(:,1); y4 = numer4(:,3);
  plot(x4, f(x4), '-k', x4, y4, 'ro');
case 320;
  numer5 = load(['example4_Nx320_K',num2str(K),'.dat']);
  x5 = numer5(:,1); y5 = numer5(:,3);
  plot(x5, f(x5), '-k', x5, y5, 'ro');
case 0;
  numer1 = load(['example4_Nx20_K',num2str(K),'.dat']);
  numer2 = load(['example4_Nx40_K',num2str(K),'.dat']);
  numer3 = load(['example4_Nx80_K',num2str(K),'.dat']);
  numer4 = load(['example4_Nx160_K',num2str(K),'.dat']);
  %numer5 = load(['example3_Nx320_K',num2str(K),'.dat']);
  x1 = numer1(:,1); y1 = numer1(:,3);
  x2 = numer2(:,1); y2 = numer2(:,3);
  x3 = numer3(:,1); y3 = numer3(:,3);
  x4 = numer4(:,1); y4 = numer4(:,3);
  %x5 = numer5(:,1); y5 = numer5(:,3);
  ex1 = f(x1);
  ex2 = f(x2);
  ex3 = f(x3);
  ex4 = f(x4);
  %ex5 = f(x5);
  plot(x4, f(x4), 'ok', x1, y1, 'o', x2, y2, '*', x3, y3, '--', x4, ...
  y4, '^');%, x5, y5, 'v');
  legend('exact', '20', '40', '80', '160');
  %nx = 320; h = 1/nx;
  %min_y = [min(y1);min(y2);min(y3);min(y4);min(y5)];
  min_y = [min(y1);min(y2);min(y3);min(y4);];%min(y5)];
otherwise
  disp('Wrong choice of plot');
end

