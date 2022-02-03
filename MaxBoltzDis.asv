function [v] = MaxBoltzDis()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global C

n = C.vth + randn();
% if n < 0
    
% syms v_0
% eqn = (C.m_n/(2*pi*C.kb*C.T))^(1/2) * exp((-C.m_n*v_0^2)/(2*C.kb*C.T)) == n;
% v = isolate(eqn,v_0);
% v = sqrt(-...
%     ((2*C.kb*C.T)/(-C.m_n))*log(...
%     ((2^0.5)*n*(pi*C.kb*C.T)^0.5)/(C.m_n^0.5)));
v = sqrt((-2*C.T*C.kb*log(n) + C.T*C.kb*log(2*pi*C.T*C.kb) - C.T*C.kb*log(C.m_n))/C.m_n);
end

