function hAtomAngular()
%Function written by Hew, Justin Kin Jun
%Started on 4th November 2020, Completed on 6th
%Written for computational MATLAB assignment for PHYS2150: Methods in
%Physics I
%Computed spherical harmonics using legendre function as well as its norm
%squared value
%Considered that the legendre function does not work for -ve values, a
%different function was inputted in its stead for different m values.
%Any input of l,m,n that fits the conditions described below will display
%the following: a 3d surface of the angular probability density function,
%and the polar plot of the function against theta

l= input('value of l:');
m= input('value of m:');
n= input('value of n:');
if (l < n) && (-l <= m) && (m <= l)  %conditions as described in paper
    disp ('valid')
else
    disp ('invalid')
    return
end 

if isinteger(l) 
    disp('valid') %Angular Quantum Number should be an integer,isinteger function was used because there is problems with the mod function
elseif  mod(l,1) ~= 0
    disp('invalid') 
    return
end
%I have included some integer conditions for input values but this is not
%necessary. The Reshape function can be adjusted such that it only works for integers. 

%Establish a linear space with spherical coordinate system
if m > 0
dx = pi/240;
alt = -pi:dx:pi;
az = 0:dx:2*pi;
[phi,theta] = meshgrid(az,alt);
Plm= legendre(l,cos(theta));
P32 = reshape(Plm(m+1,:,:), size(phi));
a = (2*l+1)*factorial(l-m); 
b = 4*pi*factorial(l+m);
C = sqrt(a/b);
Y = abs(C .* P32 .* exp(1i*m*phi)).^2;
[Xm,Ym,Zm] = sph2cart(phi, theta,(Y));
title('Probability Density','interpreter','latex');
figure (1)
surf(Xm,Ym,Zm);

elseif m==0 && l~=0 %ran into some problems with reshape function for m=0, also l=0 requires different case
dx = pi/240;
alt = -pi:dx:pi;
az = 0:dx:2*pi;
[phi,theta] = meshgrid(az,alt);
Plm= legendre(l,cos(theta));
P32 = reshape(Plm(m+1,:,:), size(phi));
a = (2*l+1); 
b = 4*pi;
C = sqrt(a/b);
Y = abs(C .* P32).^2;
[Xm,Ym,Zm] = sph2cart(phi, theta, real(Y));
figure (1)
surf(Xm,Ym,Zm);
title('Probability Density','interpreter','latex');
elseif (m < 0)
m = -m;
dx = pi/240;
alt = -pi:dx:pi;
az = 0:dx:2*pi;
[phi,theta] = meshgrid(az,alt);
P= legendre(l,cos(theta));
a = (2*l+1)*factorial(l-m); 
b = 4*pi*factorial(l+m);
C = sqrt(a/b);
P32 = reshape(P(m+1,:,:), size(phi));
Y = abs((factorial(l-m)./factorial(l+m)).*C.*(-1).^m.*P32).^2; %plugging hint and using it
[Xm,Ym,Zm] = sph2cart(phi, theta,(Y));
figure (1)
surf(Xm,Ym,Zm);
title('Probability Density ','interpreter','latex');
elseif l==0 %ran into problems with l=0 so extended here
    dx = pi/240;
alt = -pi:dx:pi;
az = 0:dx:2*pi;
[phi,theta] = meshgrid(az,alt);
Plm= legendre(l,cos(theta));
P32 = reshape(Plm(), size(phi));
a = (2*l+1); 
b = 4*pi;
C = sqrt(a/b);
Y = abs(C .* P32).^2;
[Xm,Ym,Zm] = sph2cart(phi, theta, real(Y));
figure (1)
surf(Xm,Ym,Zm);
title('Probability Density','interpreter','latex');
end
figure (2)
polarplot(theta(:,2),Y(:,2)) %choose an arbitrary number on the surface
end

%This is an additional segment which shows the code for the evaluation of
%the legendre function for spherical harmonics; to answer part b of the
%homework sheet
%I tried to use a nested function but it seemed not to work, an anonymous
%function is used instead so a file is not created
%This function is not used as the spherical harmonics is already computed
%in the prior code
function SphHarmonic()

l= input('value of l:');
m= input('value of m:');
n= input('value of n:');
dx = pi/60;
alt = -pi:dx:pi;
az = 0:dx:2*pi;
[phi,theta] = meshgrid(az,alt);
s = @(m,n,l) legendre(l,cos(theta));
P32 = reshape(s(m+1,:,:), size(phi));
a = (2*l+1)*factorial(l-m); 
b = 4*pi*factorial(l+m);
C = sqrt(a/b);
Y = (C .* P32 .* exp(1i*m*phi));
[~,~,~] = sph2cart(phi, theta,real(Y));
title('Spherical Harmonics','interpreter','latex');
figure (3)
surf(Xm,Ym,Zm);
end
