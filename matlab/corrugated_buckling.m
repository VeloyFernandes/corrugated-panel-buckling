%% Buckling of flat, trapezoidal and round corrugated panels
% Corrected version of the script in the 2023 dissertation (Appendix).
% Loads are line loads in N per metre of loaded edge.
% x runs across the corrugations, y along them.
%
% Corrections to the round corrugation block:
%   REA22  used (Rl/2*R), which MATLAB reads as (Rl/2)*R; it must be Rl/Rc
%   REA66  used D66; the membrane shear term is A66
%   RED22  used the (wrong) REA22; it must use the sheet value A22
clear; clc;

%% Sheet and panel
E  = 2.1e11;  nu = 0.3;  t = 0.00635;
G  = E/(2*(1+nu));
a  = 1.016;   b  = 1.016;

A11 = E*t/(1-nu^2);        A22 = A11;   A12 = nu*A11;   A66 = G*t;
D11 = E*t^3/(12*(1-nu^2)); D22 = D11;   D12 = nu*D11;   D66 = G*t^3/12;

%% Trapezoidal corrugation: half depth f, half period c, web angle theta
f = 0.0127;  c = 0.0508;  theta = pi/4;
l  = 2*f/sin(theta) + c - 2*f/tan(theta);                 % half developed length
I1 = 4*f*cos(theta)^2/sin(theta) + 2*c - 4*f/tan(theta);
I2 = 4*f^3/(3*sin(theta)) + 2*f^2*(c - 2*f/tan(theta));
T = equivalent_plate(A11,A12,A22,A66,D11,D12,D22,D66,c,l,I1,I2);

%% Round corrugation: crest radius R, half straight length L
R = 0.0254;  L = 0.0254;
Rc = 2*R;                                                  % half period
Rl = pi*R + 2*L;                                           % half developed length
RI1 = pi*R;
RI2 = 4*L^3/3 + 2*pi*L^2*R + 8*L*R^2 + pi*R^3;
C = equivalent_plate(A11,A12,A22,A66,D11,D12,D22,D66,Rc,Rl,RI1,RI2);

%% Euler loads of the wide strut, N_cr = pi^2 D / (K L)^2
K = [1 0.5];  names = {'pinned','fixed'};
for k = 1:2
    fprintf('\nLoad across the corrugations, %s ends [N/m]\n', names{k});
    fprintf('  flat        %12.1f\n', pi^2*D11  /(K(k)*a)^2);
    fprintf('  trapezoidal %12.1f\n', pi^2*T.D11/(K(k)*a)^2);
    fprintf('  round       %12.1f\n', pi^2*C.D11/(K(k)*a)^2);
    fprintf('Load along the corrugations, %s ends [N/m]\n', names{k});
    fprintf('  flat        %12.1f\n', pi^2*D22  /(K(k)*b)^2);
    fprintf('  trapezoidal %12.1f\n', pi^2*T.D22/(K(k)*b)^2);
    fprintf('  round       %12.1f\n', pi^2*C.D22/(K(k)*b)^2);
end

%% Yield of the section under load along the corrugations
fy = 250e6;
fprintf('\nSection yields at [N/m]: flat %.0f, trapezoidal %.0f, round %.0f\n', ...
        fy*t, fy*t*l/c, fy*t*Rl/Rc);

function P = equivalent_plate(A11,A12,A22,A66,D11,D12,D22,D66,c,l,I1,I2)
% Xia, Friswell and Saavedra Flores (2012), Int. J. Solids Struct. 49, 1453-1462
P.A11 = 2*c/(I1/A11 + I2/D11);
P.A12 = A12/A11*P.A11;
P.A22 = A12/A11*P.A12 + (l/c)*(A11*A22 - A12^2)/A11;
P.A66 = (c/l)*A66;
P.D11 = (c/l)*D11;
P.D12 = D12/D11*P.D11;
P.D22 = (I2*A22 + I1*D22)/(2*c);
P.D66 = (l/c)*D66;
end
