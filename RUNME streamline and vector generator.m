clear
clc
clf
close all

%RUN THIS FILE to generate aerofoils, streamlines and vectors

%request inputs of free stream velocity and panel code
uinf=input('Free stream velocity = ');
NACA=num2str(input('Aerofoil Panel Code = '));

%initialise index for saving plots
i=2;

%if clause for NACA 2412 case
if NACA=='2412'
    for N=[50 100 200]
        cl=[];
        i=i+1;
        for AoA=0:10
            AoA=AoA*pi/180; %angle of attack in radians
            figure(1)
            
            %generate aerofoil and find mu and cl for each angle of attack
            [x,z]=panelgen('2412',N,AoA);
            
            [mu,cl(end+1)]=liftco(x,z,N,AoA,uinf);
        end
        
        %plot lift coefficient vs angle of attack
        figure(i)
        plot(0:10,cl,'LineWidth',1.5,'Color','r')
        hold on
        
        %load and plot the XFOIL data
        xfoil=readmatrix('xf-naca2412-il-1000000.txt','NumHeaderLines',12);
        xfoil=xfoil((xfoil(:,1)<=10)&(xfoil(:,1)>=0),:);
        xfoil=xfoil(:,1:2);
        plot(xfoil(:,1),xfoil(:,2),'LineWidth',1.5,'Color','k')
        
        %axis labels, title and legend
        words=['Cl vs alpha for NACA 2412 aerofoil with ',num2str(N),' panels'];
        title(words,'FontSize',16);
        xlabel('Angle of attack (deg)','FontSize',14);
        ylabel('Coefficient of lift','FontSize',14);
        legend('XFOIL data','Computed values','FontSize',12)
        hold off
        
        %save plots
        saveas(figure(i),[words,'.png']);
    end
else
    %ask for inputs if foil is not NACA 2412
    AoA=input('Angle of attack = ');
    AoA=AoA*pi/180; %angle of attack in radians
    N=input('number of panels = ');
    
    %generate aerofoil and calculate mu and lift coefficient
    figure(1)
    [x,z]=panelgen('2415',N,AoA);
    [mu,cl]=liftco(x,z,N,AoA,uinf);
end

%calculate velocity at each point
[x2,z2,U,V]=velocities(x,z,mu,uinf,N,AoA);

%plot streamlines
figure(1)
strm=streamslice(x2,z2,U,V);

%formatting of streamline plot
set(strm,'LineWidth',1)
xlim([-0.2,1.2]);
ylim([-0.7,0.7]);
manywords=['NACA ',NACA,' aerofoil with ',num2str(N),' panels'];
fewerwords=['Angle of attack = ',num2str(AoA*180/pi),' deg, free stream velocity = ',num2str(uinf),' m/s'];
title(['Streamlines of ',manywords],'FontSize',16);
subtitle(fewerwords,'FontSize',16);
xlabel('X coordinate relative to leading edge','FontSize',14);
ylabel('Y coordinate relative to leading edge','FontSize',14);

%save streamline plot
filename=['NACA ',NACA,' aerofoil with ',num2str(N),' panels ',num2str(AoA*180/pi),' deg'];
saveas(figure(1),['Streamlines of ',filename,'.png']);

%vector plot
figure(2)
plot(x(1:end-1),z(1:end-1),'m','LineWidth',2);
xlim([-0.2,1.2]);
ylim([-0.7,0.7]);
hold on
quiv=quiver(x2,z2,U,V);

%format vector plot
set(quiv,'LineWidth',1)
xlim([-0.2,1.2]);
ylim([-0.7,0.7]);
title(['Velocity vectors of ',manywords],'FontSize',16);
subtitle(fewerwords,'FontSize',16);
xlabel('X coordinate relative to leading edge','FontSize',14);
ylabel('Y coordinate relative to leading edge','FontSize',14);

%save vector plot
saveas(figure(2),['Vectors of ',filename,'.png']);
hold off