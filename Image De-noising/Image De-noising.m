clc

a=imread('cameraman.tif');
ref=im2double(a);
noisy=imnoise(ref,'gaussian',0.05);
U=noisy;

figure(1);
subplot(1,2,1);
imshow(noisy);
title('Noisy image')
drawnow;

dt=0.0025;
r=1/4;
dx=0.01;
disp(dx);
dy=dx;
loop=1/dt;

t=1;
%t=input('Enter time at which solution is needed '); % example t=0.1, 1 etc 
% Evolving upto 't'
 Upc=translateImage(U,1,0); %up
Umc=translateImage(U,-1,0);%down
Ucp=translateImage(U,0,1);%left
Ucm=translateImage(U,0,-1);%right
 UY=(Upc-Umc)/(2*dy);
 UX=(Ucp-Ucm)/(2*dx);
for i=1:100
    Upc=translateImage(U,1,0); %up
    Umc=translateImage(U,-1,0);%down
    Ucp=translateImage(U,0,1);%left
    Ucm=translateImage(U,0,-1);%right
    
    Uy=(Upc-Umc)/(2*dy);
    Ux=(Ucp-Ucm)/(2*dx);
    
    
    mod=((Ux.^2)+(Uy.^2));
    modf=((UX.^2)+(UY.^2));
    
    %D = 1./(exp(sqrt(mod)));
    %D=1./((mod+1));
    D = 1./(1 + modf);
    D1=translateImage(D,1,0); %up
    D2=translateImage(D,-1,0);%down
    D3=translateImage(D,0,1);%left
    D4=translateImage(D,0,-1);%right
    
    A=((D+D3).*(Ucp-U)) - ((D+D4).*(U-Ucm)) + ((D+D1).*(Upc-U)) -((D+D2).*(U-Umc));
    U=(0.5*r*A)+U;
    
    figure(1);
    subplot(1,2,2);
    imshow(U);
    title('De-noised image')
    drawnow;
end

    