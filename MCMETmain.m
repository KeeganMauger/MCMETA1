% ELEC 4700
% Assignment 1: Monte-Carlo Modeling of Electron Transport
% Keegan Mauger

% Bounding Area: 200nm x 100nm

% rectangle('Position',[1 2 5 6])
% axis([0 10 0 10])

% Initialization
set(0,'DefaultFigureWindowStyle','docked')
set(0,'defaultaxesfontsize',20)
set(0,'defaultaxesfontname','Times New Roman')
set(0,'DefaultLineLineWidth', 2);

run WC
global C

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665;                      % metres (32.1740 ft) per s²
C.m_n = 0.26 * C.m_0;               % effective electron mass
C.am = 1.66053892e-27;              % atomic mass unit
C.t = 300;                          % temperature

temp = 300;

rectangle('Position',[0 0 200 100])
hold on
% Initialize particle positions

% i = 0;
% j = 0;
% N = 9;
% r = -1e-8;
% s = 1e-8;
% 
% for j=1:N
%     
%     for i = 1:N
%         plot(((j*10^-8)*2)+(r + (s - r).*rand(1,1)),((i*10^-8)+(r + (s - r).*rand(1,1))),'b.')
%     end
%     j = j+1;
% end
% k = 0 + (200-0).*rand(1,1);
% l = 0 + (100-0).*rand(1,1);

%--------------------------------------------------------------------------
% Initializing Positions
%--------------------------------------------------------------------------


N = 10;
i = 0;
j = 0;
k = 0;
l = 0;

for j=1:N
    for i = 1:N
        m(i,j) = 0 + (200-0).*rand(1,1);
        n(i,j) = 0 + (100-0).*rand(1,1);
    end
end

for k=1:N
    
    for l=1:N

        plot(m(l,k),n(l,k),'b.')
        hold on
    end
end

% Thermal Velocity and Direction

vth = sqrt(C.kb * temp / C.m_n);
vx = zeros(N);
vy = zeros(N);
theta = zeros(N);

r = 0;
s = 0;

for s=1:N;
    for r=1:N
        vx(r,s) = vth;
        vy(r,s) = vth;
        theta(r,s) = 0 + (359-0).*rand(1,1);
    end
end











