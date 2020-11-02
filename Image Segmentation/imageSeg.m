

clc
clear all
Image=imread('shapes.jpg');
ref=im2double(Image);
getd = @(p)path(p,path);
getd('curve/toolbox_signal/');
getd('curve/toolbox_general/');
getd('curve/toolbox_graph/');
[sizeY, sizeX] = size(Image);
[Y,X] = meshgrid(1:sizeX,1:sizeY);
Xc = sizeY/2;
Yc = sizeX/2;
r = min(sizeX,sizeY)/2;
phi1 = sqrt( (X-Xc).^2 + (Y-Yc).^2 ) - r;  % surface. zero level line would be the circle with center (c1,c2) and radius r
clf; % clear current figure
u=phi1;

subplot(1,2,1);
plot_levelset(u,0,ref);
title('original')
drawnow;
dt=0.00001;
lambda = 0.001;
dx = 0.01;
dy = 0.01;
U = u;
alpha = 100;
[um,un] = size(U);
A = ones(um,un);
g = ref;
for n=1:5
    
    Umc=translateImage(U,-1,0); % Shifts down, ie, In place of U(i,j) we will have U(i-1,j)
    Upc=translateImage(U,1,0); % Shifts up, ie, In place of U(i,j) we will have   U(i+1,j)
    Ucm=translateImage(U,0,-1);% Shifts right, ie, In place of U(i,j) we will have U(i,j-1)
    Ucp=translateImage(U,0,1); % Shifts left, ie, In place of U(i,j) we will have U(i,j+1)
    Uxyp=translateImage(U,1,1); 
    Uxym=translateImage(U,-1,-1); 
    Uxypm=translateImage(U,1,-1); 
    Uxymp=translateImage(U,-1,1); 
    
    Ux = (Ucp - U)./(dx);
    Uy = (Upc - U)./(dy);
    Uxx = (Ucp + Ucm - 2*U)./(dx^2);        % (uxxuy2 -2uxyuxuy + uyyux2 ) / (ux2 + uy2)3/2 
    Uyy = (Upc + Umc - 2*U)./(dy^2);
    Uxy = (Uxym + Uxyp - Uxypm -Uxymp)./(4*dx*dy);
    k = ( (Uxx.*Uy).^2 - 2*(Uxy.*Uxy) + (Uyy.*Ux).^2 )./((Ux.^2 + Uy.^2).^1.5);
    dU = (Ux + Uy);
    %g = 1./(A + (dU.^2)*(1/(lambda^2)) );
    gcp=translateImage(g,0,1);
    gpc=translateImage(g,1,0);
    gx = (gcp-g)./(dx);
    gy = (gpc-g)./(dy);
    deltaX = (Ucp - Ucm)./(2*dx);
    deltaY = (Upc - Umc)./(2*dy);
    deltaYP = (Upc - U)./(dx);
    deltaYM = (U - Umc)./(dx);
    deltaXP = (Ucp - U)./(dy);
    deltaXM = (U - Ucm)./(dy);
    
    B = zeros(um,un);
    dPlus = ( max(deltaXM.*U,B).^2 + min(deltaXP.*U,B).^2 + max(deltaYM.*U,B).^2 + min(deltaYP.*U,B).^2 ).^(0.5);
    dMinus =( max(deltaXP.*U,B).^2 + min(deltaXM.*U,B).^2 + max(deltaYP.*U,B).^2 + min(deltaYM.*U,B).^2 ).^(0.5);
    t1 = (deltaX.*U).^2 + (deltaY.*U).^2;
    t2 = max(g,B).*dPlus + min(g,B).*dMinus;
    t3 = max(gx,B).*deltaXM.*U + min(gx,B).*deltaXP.*U;
    t4 = max(gy,B).*deltaYM.*U + min(gy,B).*deltaYP.*U;
    U = U + dt*(g.*k.*( (t1).^(0.5)) + alpha*(t2).*U + t3 + t4);
    subplot(1,2,2);
    plot_levelset(U,0,ref);
    drawnow;
end

