% This program removes noise from an image using linear second order
% diffusion filter u_t=u_xx+u_yy with homogeneous neumann boundary conditions
clc
clear all



a=imread('cameraman.tif');
ref=im2double(a);
%noisy=imnoise(ref,'gaussian',0.05);
U=ref;

t=input('Enter time at which solution is needed '); % example t=0.1, 1 etc 
dt=0.1;

nitr=t/dt; % No of iterations
r=dt;


verbose=1;
if verbose
    figure(verbose);
    subplot(1,3,1);
    imshow(ref);
    title('Original image')
    drawnow;
end

for n=1:nitr   % Evolving upto 't'

    Umc=translateImage(U,-1,0); % Shifts down, ie, In place of U(i,j) we will have U(i-1,j)
    Upc=translateImage(U,1,0); % Shifts up, ie, In place of U(i,j) we will have   U(i+1,j)
    Ucm=translateImage(U,0,-1);% Shifts right, ie, In place of U(i,j) we will have U(i,j-1)
    Ucp=translateImage(U,0,1); % Shifts left, ie, In place of U(i,j) we will have U(i,j+1)
    
    U=U+r*(Upc+Umc+Ucp+Ucm-4*U); % updating U
    
    if verbose
    figure(verbose);
    subplot(1,3,2);
    imshow(U);
    title('Blurred image')
    drawnow;
    end
    %pause(1) % pause the display for two seconds

end

    %Filtered image at given 't'
dx = 1;
dy = 1;


for n=1:200   % Evolving upto 't'

    Umc=translateImage(U,-1,0); % Shifts down, ie, In place of U(i,j) we will have U(i-1,j)
    Upc=translateImage(U,1,0); % Shifts up, ie, In place of U(i,j) we will have   U(i+1,j)
    Ucm=translateImage(U,0,-1);% Shifts right, ie, In place of U(i,j) we will have U(i,j-1)
    Ucp=translateImage(U,0,1); % Shifts left, ie, In place of U(i,j) we will have U(i,j+1)
    deltaYP = (Upc - U)./(dx);
    deltaYM = (U - Umc)./(dx);
    deltaXP = (Ucp - U)./(dy);
    deltaXM = (U - Ucm)./(dy);
    Ux = m( deltaXP , deltaXM );
    Uy = m( deltaYP , deltaYM );
    c1 = ( (Ux).^2 + (Uy).^2 ).^(0.5);
    c2 = (Upc+Umc+Ucp+Ucm-4*U)./(dx);
    %c3 = sin(c2);
    %c3 = c2./abs(c2);
    c3 = tanh(c2);
    %c3 = (abs(c2)).^(1/3);
    c3(isnan(c3))=0;
    U=U-dt*(c1.*c3); % updating U
    
    if verbose
    figure(verbose);
    subplot(1,3,3);
    imshow(U);
    title('De - Blurred image')
    drawnow;
    end
    %pause(1) % pause the display for two seconds

end

function O = m(a,b)
    O = a.*b;
    for i=1:256
        for j = 1:256
            if O(i,j)<=0
                O(i,j)=0;
            else
                s = abs(a(i,j))/a(i,j);
                O(i,j)=s*min(abs(a(i,j)),abs(b(i,j)));
            end
        end
    end
end
