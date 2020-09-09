% This program removes blur from an image using shock filter model proposed by Osher and Rudin
% diffusion filter u_t+(sign(u_x)sign(u_xx))u_x=0 with u(x,0)=cos(x)
clc
clear all


u=imread('cameraman.jpg');  % Blurred signal: smooth function
nitr=input('Enter no of iterations ');
U=im2double(u);
verbose=1;
if verbose
    figure(verbose);
    subplot(1,2,1);
    imshow(u)
    title('Blurred signal')
    drawnow;
end
for n=1:nitr
    
  
    Umc=translateImage(U,-1,0); % Shifts down, ie, In place of U(i,j) we will have U(i-1,j)
    Upc=translateImage(U,1,0); % Shifts up, ie, In place of U(i,j) we will have   U(i+1,j)
    Ucm=translateImage(U,0,-1);% Shifts right, ie, In place of U(i,j) we will have U(i,j-1)
    Ucp=translateImage(U,0,1); % Shifts left, ie, In place of U(i,j) we will have U(i,j+1)
    A=(U-Umc); 
    B=(Upc-U); 
    C=(U-Ucm); 
    D=(Ucp-U); 
    K=MIN(B,A);
    L=MIN(D,C);
    U=U-(((0.1)*sqrt((K.^2)+(L.^2))).*(sign((Umc+Upc+Ucm+Ucp-(4*U)))));
    
    
   if verbose
    figure(verbose);
    subplot(1,2,2);
    imshow(U)
    title(n)
    drawnow;
   end

end

function Y=MIN(a,b)
    c=sign(a.*b);
    c(c<=0)=0;
    Y=sign(a).*min((sign(a).*a),(sign(b).*b));
    Y=Y.*c;
end
