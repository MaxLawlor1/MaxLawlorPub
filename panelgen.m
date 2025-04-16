function[xf, z]=panelgen(NACA,N,AoA)
%set numbers of NACA code as m, p, t
m=str2double(NACA(1))/100;
p=str2double(NACA(2))/10;
t=str2double(NACA(3:4))/100;

xa=1-0.5*(1-cos(2*pi*([0:1:N/2]/N)));
xa=sort(xa,2);
%find x and y below max thickness
x1=xa(xa<p);
yc1=(m/(p^2))*(2*p*x1-x1.^2);
dycdx1=(2*m/(p^2))*(p-x1);
%find x and y above max thickness
x2=xa(xa>p);
yc2=m/((1-p)^2)*((1-2*p)+2*p*x2-x2.^2);
dycdx2=2*m/((1-p)^2)*(p-x2);


%define theta
%theta=atan([unique(dycdx1');unique(dycdx2')])'
theta=atan([dycdx1 dycdx2]);

x=[x1 x2];
yc=[yc1 yc2];

%define yt, xU, xL, zU, zL
yt=5*t*(0.2969*x.^0.5-0.126*x-0.3516*x.^2+0.2843*x.^3-0.1036*x.^4);

xU=x-yt.*sin(theta);
xL=x+yt.*sin(theta);

zU=yc+yt.*cos(theta);
zL=yc-yt.*cos(theta);

%form arrays x and z
%form arrays x and z
xf=[flip(xL) xU(2:end)];
z=[flip(zL) zU(2:end)];

%plot for streamlines
figure(1)
plot(xf,z,'m','LineWidth',2);
xlim([-0.2,1.2]);
ylim([-0.7,0.7]);

%plot for velocity vectors
%foil=gca;
%figure(2)
%copyobj(foil,figure(2));
%hold on

bignumber=10^9;
xf(N+2)=bignumber;
z(N+2)=bignumber*tan(AoA);