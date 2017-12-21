K = 3;
limiter = 0;
numer1 = load(['example9_Nx20_K',num2str(K), '_PP',num2str(limiter),'.dat']);
numer2 = load(['example9_Nx40_K',num2str(K), '_PP',num2str(limiter),'.dat']);
numer3 = load(['example9_Nx80_K',num2str(K), '_PP',num2str(limiter),'.dat']);
numer4 = load(['example9_Nx160_K',num2str(K),'_PP',num2str(limiter),'.dat']);
numer5 = load(['example9_Nx320_K',num2str(K),'_PP',num2str(limiter),'.dat']);
x1 = numer1(:,1); y1 = numer1(:,2);
x2 = numer2(:,1); y2 = numer2(:,2);
x3 = numer3(:,1); y3 = numer3(:,2);
x4 = numer4(:,1); y4 = numer4(:,2);
x5 = numer5(:,1); y5 = numer5(:,2);
err = zeros(4,3);
err(1,1) = norm(y5(1:16:end)-y1,2)/sqrt(20);
err(2,1) = norm(y5(1:8:end)-y2,2)/sqrt(40);
err(3,1) = norm(y5(1:4:end)-y3,2)/sqrt(80);
err(4,1) = norm(y5(1:2:end)-y4,2)/sqrt(160);
err(1,2) = norm(y5(1:16:end)-y1,1)/20;
err(2,2) = norm(y5(1:8:end)-y2,1)/40;
err(3,2) = norm(y5(1:4:end)-y3,1)/80;
err(4,2) = norm(y5(1:2:end)-y4,1)/160;
err(1,3) = norm(y5(1:16:end)-y1,'inf');
err(2,3) = norm(y5(1:8:end)-y2,'inf');
err(3,3) = norm(y5(1:4:end)-y3,'inf');
err(4,3) = norm(y5(1:2:end)-y4,'inf');
err
log2(err(1,:)./err(2,:))
log2(err(2,:)./err(3,:))
log2(err(3,:)./err(4,:))

