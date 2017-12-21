function plot_ex1(K, n, limiter);
% para: n, P_n polynomial;

addpath('../src/')
format long ;
f = @(x) (sin(4*x)-8*sin(2*x)+12*x)/32;
err = zeros(5,3);

switch n;
case 20;
  numer1 = load(['example1_Nx20_K',num2str(K),'.dat']);
  x1 = numer1(:,1); y1 = numer1(:,3);
  plot(exact(:,1), exact(:,2), '-k', x1, y1, 'ro');
  ex1 = f(x1);
  nx = 20; h = 1/nx;
  err(1,1) = cal_norm(ex1-y1, numer1(:,2), h, 2);
  err(1,2) = cal_norm(ex1-y1, numer1(:,2), h, 1);
  err(1,3) = cal_norm(ex1-y1, numer1(:,2), h, 'inf');
  disp(err);
case 40;
  numer2 = load(['example1_Nx40_K',num2str(K),'.dat']);
  x2 = numer2(:,1); y2 = numer2(:,3);
  plot(exact(:,1), exact(:,2), '-k', x2, y2, 'ro');
  ex2 = f(x2);
  nx = 40; h = 1/nx;
  err(2,1) = cal_norm(ex2-y2, numer2(:,2), h, 2);
  err(2,2) = cal_norm(ex2-y2, numer2(:,2), h, 1);
  err(2,3) = cal_norm(ex2-y2, numer2(:,2), h, 'inf');
  disp(err);
case 80;
  numer3 = load(['example1_Nx80_K',num2str(K),'.dat']);
  x3 = numer3(:,1); y3 = numer3(:,3);
  plot(exact(:,1), exact(:,2), '-k', x3, y3, 'ro');
  ex3 = f(x3);
  nx = 80; h = 1/nx;
  err(3,1) = cal_norm(ex3-y3, numer3(:,2), h, 2);
  err(3,2) = cal_norm(ex3-y3, numer3(:,2), h, 1);
  err(3,3) = cal_norm(ex3-y3, numer3(:,2), h, 'inf');
  disp(err);
case 160;
  numer4 = load(['example1_Nx160_K',num2str(K),'.dat']);
  x4 = numer4(:,1); y4 = numer4(:,3);
  plot(exact(:,1), exact(:,2), '-k', x4, y4, 'ro');
  ex4 = f(x4);
  nx = 160; h = 1/nx;
  err(4,1) = cal_norm(ex4-y4, numer4(:,2), h, 2);
  err(4,2) = cal_norm(ex4-y4, numer4(:,2), h, 1);
  err(4,3) = cal_norm(ex4-y4, numer4(:,2), h, 'inf');
  disp(err);
case 320;
  numer5 = load(['example1_Nx320_K',num2str(K),'.dat']);
  x5 = numer5(:,1); y5 = numer5(:,3);
  plot(exact(:,1), exact(:,2), '-k', x5, y5, 'ro');
  ex5 = f(x5);
  nx = 320; h = 1/nx;
  err(5,1) = cal_norm(ex5-y5, numer5(:,2), h, 2);
  err(5,2) = cal_norm(ex5-y5, numer5(:,2), h, 1);
  err(5,3) = cal_norm(ex5-y5, numer5(:,2), h, 'inf');
  disp(err);
case 0;
  numer1 = load(['example1_Nx20_K',num2str(K) ,'_PP',num2str(limiter),'.dat']);
  numer2 = load(['example1_Nx40_K',num2str(K) ,'_PP',num2str(limiter),'.dat']);
  numer3 = load(['example1_Nx80_K',num2str(K) ,'_PP',num2str(limiter),'.dat']);
  numer4 = load(['example1_Nx160_K',num2str(K),'_PP',num2str(limiter),'.dat']);
  numer5 = load(['example1_Nx320_K',num2str(K),'_PP',num2str(limiter),'.dat']);
  x1 = numer1(:,1); y1 = numer1(:,3);
  x2 = numer2(:,1); y2 = numer2(:,3);
  x3 = numer3(:,1); y3 = numer3(:,3);
  x4 = numer4(:,1); y4 = numer4(:,3);
  x5 = numer5(:,1); y5 = numer5(:,3);
  ex1 = f(x1);
  ex2 = f(x2);
  ex3 = f(x3);
  ex4 = f(x4);
  ex5 = f(x5);
  %plot(exact(:,1), exact(:,2), '-k', x1, y1, 'o', x2, y2, '*', x3, y3, '--', x4, ...
  %y4, '^', x5, y5, 'v');
  %cal_order
  nx = 20; h = 1/nx;
  err(1,1) = cal_norm(ex1-y1, numer1(:,2), h, 2);
  err(1,2) = cal_norm(ex1-y1, numer1(:,2), h, 1);
  err(1,3) = cal_norm(ex1-y1, numer1(:,2), h, 'inf');
  nx = 40; h = 1/nx;
  err(2,1) = cal_norm(ex2-y2, numer2(:,2), h, 2);
  err(2,2) = cal_norm(ex2-y2, numer2(:,2), h, 1);
  err(2,3) = cal_norm(ex2-y2, numer2(:,2), h, 'inf');
  nx = 80; h = 1/nx;
  err(3,1) = cal_norm(ex3-y3, numer3(:,2), h, 2);
  err(3,2) = cal_norm(ex3-y3, numer3(:,2), h, 1);
  err(3,3) = cal_norm(ex3-y3, numer3(:,2), h, 'inf');
  nx = 160; h = 1/nx;
  err(4,1) = cal_norm(ex4-y4, numer4(:,2), h, 2);
  err(4,2) = cal_norm(ex4-y4, numer4(:,2), h, 1);
  err(4,3) = cal_norm(ex4-y4, numer4(:,2), h, 'inf');
  nx = 320; h = 1/nx;
  err(5,1) = cal_norm(ex5-y5, numer5(:,2), h, 2);
  err(5,2) = cal_norm(ex5-y5, numer5(:,2), h, 1);
  err(5,3) = cal_norm(ex5-y5, numer5(:,2), h, 'inf');
  min_y = [min(y1);min(y2);min(y3);min(y4);min(y5)];
  order = [
  zeros(1,3);
  log2(err(1,:)./err(2,:));
  log2(err(2,:)./err(3,:));
  log2(err(3,:)./err(4,:));
  log2(err(4,:)./err(5,:))];
  N = [20;40;80;160;320];
  %diary table1.dat
  %diary on;
  for n = 1:5
    fprintf('%3d ', N(n));
    for i = 1:3
      fprintf('%.3e %.2f ', err(n,i), order(n,i));
    end
    fprintf('%.3e\n', min_y(n));
  end
  %diary off;

otherwise
  disp('Wrong choice of plot');
end

