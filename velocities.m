function[x2,z2,U,V]=velocities(x,z,mu,uinf,N,AoA)

%create meshgrid of points for streamline plot
[x2,z2]=meshgrid(linspace(-0.2,1.2,200), linspace(-0.7,0.7,200));

%initialise velocity matrices
U=zeros(200);
V=zeros(200);

%remove points inside the aerofoil
for i=1:200
    for j=1:200
        for k=1:N+1    
            if inpolygon(x2(i,j),z2(i,j),x,z)
                    x2(i,j)=NaN;
                    z2(i,j)=NaN;
            end
            %find velocities at meshgrid points
            [u2(i,j),v2(i,j)]=cdoublet([x2(i,j),z2(i,j)],[x(k),z(k)],[x(k+1),z(k+1)]);
            U(i,j)=mu(k)*u2(i,j)+U(i,j);
            V(i,j)=mu(k)*v2(i,j)+V(i,j);
        end
    end
end

%finish finding velocities by adding free stream velocity influenced by
%angle of attack
U=U+uinf*cos(AoA);
V=V+uinf*sin(AoA);