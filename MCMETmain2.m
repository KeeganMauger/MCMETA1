


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
C.g = 9.80665;                      % metres (32.1740 ft) per s�
C.m_n = 0.26 * C.m_0;               % effective electron mass
C.am = 1.66053892e-27;              % atomic mass unit
C.T = 300;
C.vth = sqrt(C.kb * C.T / C.m_n);


temp = C.T;

subplot(2,1,1);
rectangle('Position',[0 0 200e-9 100e-9])
hold on

%--------------------------------------------------------------------------
% Initializing Positions
%--------------------------------------------------------------------------


N = 100;        % Number of electrons
i = 0;
j = 0;

for i=1:N
    px(i) = 0 + (200e-9 - 0).*rand(1,1);
    py(i) = 0 + (100e-9 - 0).*rand(1,1);
    subplot(2,1,1);
    plot(px(i),py(i),'b.')
    hold on
end

%--------------------------------------------------------------------------
% Thermal Velocity and Direction
%--------------------------------------------------------------------------

vth = C.vth;

for j=1:N
    v0(j) = MaxBoltzDis();                                % Velocity of electron
    theta(j) = 0 + (360 - 0).*rand(1,1);        % Angle of electron
    if theta(j) == 360
        theta(j) = 0;
    end
    vx(j) = v0(j)*cos(theta(j));                % Velocity in x axis
    vy(j) = v0(j)*sin(theta(j));                % Velocity in y axis
end



t = 0;
T(1) = 0;
dt = 1e-14;     % time step
px_prev = 0;
py_prev = 0;
T_prev = 0;
P_scat = 0;

for t=2:100
    for k=1:N
        
        
        
        px_prev(k) = px(k);
        px(k) = px(k) + vx(k)*dt;
        py_prev(k) = py(k);
        py(k) = py(k) + vy(k)*dt;
        
        if py(k) >= 100e-9 || py(k) <= 0
            [theta(k),vx(k),vy(k)] = SpecRef(theta(k),vx(k),vy(k));
        end
        if px(k) >= 200e-9
            px(k) = 0;
            px_prev(k) = px(k);
        elseif px(k) <= 0
            px(k) = 200e-9;
            px_prev(k) = px(k);
        else
            px(k) = px(k);
        end
        
        v(k) = sqrt(vx(k)^2 + vy(k)^2);
        v2(k) = v(k).*v(k);
        
        subplot(2,1,1);
        plot([px_prev(k) px(k)],[py_prev(k) py(k)],'b')
        hold on
    end
    
    KE = 0.5 * C.m_n * mean(v2);
    KEtot = KE * N;
    T_prev = T;
    T(t) = KEtot / N / C.kb;
    subplot(2,1,2);
    plot(t, T, 'b.')
    hold on
    
    pause(0.1)
end

