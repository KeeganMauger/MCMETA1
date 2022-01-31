% ELEC 4700
% Assignment 1: Monte-Carlo Modeling of Electron Transport
% Keegan Mauger

% Bounding Area: 200nm x 100nm

% rectangle('Position',[1 2 5 6])
% axis([0 10 0 10])

set(0, 'DefaultFigureWindowStyle', 'docked')
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


rectangle('Position',[0 0 200e-9 100e-9])