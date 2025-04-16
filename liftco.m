function[mu,cl]=liftco(x,z,N,AoA,uinf)

%Creating p, a matrix of the coords of the centres of each panel
for i=1:N+1
    p(i,1)=(x(i+1)-x(i))/2+x(i);
    p(i,2)=(z(i+1)-z(i))/2+z(i);
end

%Loop to find u(i,j) and v(i,j)
for i=1:N+1
        beta(i)=mod(atan2((z(i+1)-z(i)),(x(i+1)-x(i))),2*pi); %find beta
    for j=1:N+1
            [u(i,j),v(i,j)]=cdoublet([p(i,1),p(i,2)],[x(j),z(j)],[x(j+1),z(j+1)]);
            A(i,j)=v(i,j)*cos(beta(i))-u(i,j)*sin(beta(i)); %create matrix A for linear system
    end
end

%Create and solve matrix system
B=-uinf*sin(AoA-beta);

%apply Kutta condition
B(end)=0;
A(N+1,:)=0;
A(N+1,[1 N+1])=1;
A(N+1,N)=-1;
mu=A\B';


%coefficient of lift
cl=-2*mu(N+1)/uinf;