function [tt,pp,rr,ss]=spatial_simulation(ap,ar,g)

% Takes parameters for cost of production (ap), cost of resistance (ar), and density (g)
% Performs reaction-diffusion simulation in a 2D grid using Crank-Nicholson
% Returns genotype abundances through the grid at time points

N=51;           % Length of grid
Tend = 210;     % Time of simulation
D = 0.0001;     % Diffusion constant

% Space and time intervals
dx = 1/(N-1);
dt = 0.1*(0.5*dx^2/D);
k=dt/(2*dx^2);
x = [0:dx:1]; y = [0:dx:1];

% Generates initial conditions
p=rand(N,N)<0.01;
r=rand(N,N)<0.01;
s=rand(N,N)<0.01;
p = reshape(p,N^2,1);
r = reshape(r,N^2,1);
s = reshape(s,N^2,1);
time = 0;

% Matrices to implement diffusion (Crank-Nicholson)
A1D = diffmat(N); Imat = eye(N);
A2D = kron(A1D,Imat) + kron(Imat, A1D);
Amat = eye(N^2) - D*k*A2D;

% Run reaction-diffusion simulation
Nt=Tend/dt;
n=0;
for i=1:Nt
   
   % Diffusion
   p = Amat\(p+D*k*A2D*p);
   r = Amat\(r+D*k*A2D*r);
   s = Amat\(s+D*k*A2D*s);
   
   % Reaction
   p= p + dt * (p.*(1-ap-p*(1-ap)-r*(1-ar)-s));
   r= r + dt * (r.*(1-ar-p*(1-ap)-r*(1-ar)-s));
   s= s + dt * (s.*(1-p*(1-ap+g)-r*(1-ar)-s));
   
   time = time+dt;
   
   % Record values
   if time/log(2)>n
       tt(n+1) = n;
       pp(n+1,:)=p;
       rr(n+1,:)=r;
       ss(n+1,:)=s;
       n=n+1;
   end
   
end

function [A] = diffmat(n)

A = zeros(n);

for i=1:n
    A(i,i)=-2.0;
end

for i=2:n
    A(i,i-1) = 1.0;
    A(i-1,i) = 1.0;
end

return

end

end