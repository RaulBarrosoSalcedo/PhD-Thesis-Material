function [f,LM1,LM2,delta] = LM_cubo
% ==============================================================
% CALCULO DE L_M vs FRECUENCIA
%
% Cubo conductor de aluminio
%
% MATLAB R2009a
%
% --------------------------------------------------------------
% Primera expresion:
%
% LM1 = mu0*(b - 2*delta)^2/b
%
% Segunda expresion:
%
% LM2 = real{ mu0*b*(1-Suma1)/Suma2 }
%
% --------------------------------------------------------------
%
% Frecuencia:
%
%       1.5 kHz <= f <= 1.5 GHz
%
% Material:
%
%       sigma = 3.8e7 S/m
%       mu    = mu0
%
% Cubo:
%
%       b = 3 mm
%
% ==============================================================


clear;
clc;


% ==============================================================
% PARAMETROS FISICOS
% ==============================================================

sigma = 3.8e7;             % Conductividad [S/m]

mu0 = 4*pi*1e-7;            % Permeabilidad del vacio [H/m]

eps0 = 8.854187817e-12;     % Permitividad del vacio [F/m]

b = 3e-3;                   % Dimensiones del cubo [m]


% ==============================================================
% RANGO DE FRECUENCIA
% ==============================================================

fmin = 1.5e3;               % 1.5 kHz
fmax = 1.5e9;               % 1.5 GHz

Nf = 400;                   % Numero de puntos

% Distribucion logaritmica
f = logspace(log10(fmin),log10(fmax),Nf);


% ==============================================================
% NUMERO DE TERMINOS DE LA SUMATORIA
% ==============================================================

% La expresion teorica tiene:
%
% sum(n=1,infinito)
%
% Aqui usamos Nmode terminos.
%
% Posteriormente se recomienda comprobar la convergencia
% utilizando Nmode = 20, 40, 80, etc.

Nmode = 40;


% ==============================================================
% RESERVAR MEMORIA
% ==============================================================

LM1   = zeros(1,Nf);
LM2   = zeros(1,Nf);
delta = zeros(1,Nf);


% ==============================================================
% CALCULO PARA CADA FRECUENCIA
% ==============================================================

for ii = 1:Nf

    % ----------------------------------------------------------
    % Frecuencia angular
    % ----------------------------------------------------------

    omega = 2*pi*f(ii);


    % ----------------------------------------------------------
    % PROFUNDIDAD DE PENETRACION
    % ----------------------------------------------------------

    delta(ii) = sqrt( ...
        2/(omega*mu0*sigma) );


    % ==========================================================
    % PRIMERA EXPRESION DE L_M
    %
    % L_M = mu0*(b-2*delta)^2/b
    % ==========================================================

    LM1(ii) = ...
        mu0*(b-2*delta(ii))^2/b;


    % ==========================================================
    % SEGUNDA EXPRESION DE L_M
    % ==========================================================

    LM2(ii) = calcula_LM2( ...
        omega,b,sigma,mu0,Nmode);

end


% ==============================================================
% GRAFICA PRINCIPAL: L_M VS f
% ==============================================================

figure;

loglog(f,abs(LM1),'b-','LineWidth',2);
hold on;

loglog(f,abs(LM2),'r--','LineWidth',2);

grid on;

xlabel('Frecuencia f [Hz]');
ylabel('L_M [H]');

title('Inductancia equivalente L_M vs frecuencia');

legend('Primera expresion','Segunda expresion');


% ==============================================================
% TICKS DEL EJE X
% ==============================================================

set(gca,'XTick',[1e3 1e4 1e5 1e6 1e7 1e8 1e9]);


% ==============================================================
% GRAFICA DE DELTA
% ==============================================================

figure;

loglog(f,delta,'k-','LineWidth',2);

grid on;

xlabel('Frecuencia f [Hz]');
ylabel('\delta_c [m]');

title('Profundidad de penetracion');


% ==============================================================
% GRAFICA DE LA RELACION LM2/LM1
% ==============================================================

figure;

semilogx(f,LM2./LM1,'g-','LineWidth',2);

grid on;

xlabel('Frecuencia f [Hz]');
ylabel('L_{M,2}/L_{M,1}');

title('Relacion entre las dos expresiones de L_M');


% ==============================================================
% ERROR RELATIVO
% ==============================================================

error_relativo = ...
    abs(LM2-LM1)./abs(LM1)*100;


figure;

semilogx(f,error_relativo,'m-','LineWidth',2);

grid on;

xlabel('Frecuencia f [Hz]');
ylabel('Error relativo [%]');

title('Error relativo entre las dos expresiones');


% ==============================================================
% RESULTADOS NUMERICOS
% ==============================================================

fprintf('\n');
fprintf('======================================================\n');
fprintf(' RESULTADOS PARA L_M\n');
fprintf('======================================================\n');


fprintf('\n');
fprintf('Frecuencia inicial:\n');

fprintf('f     = %e Hz\n',f(1));

fprintf('delta = %e m\n',delta(1));

fprintf('LM1   = %e H\n',LM1(1));

fprintf('LM2   = %e H\n',LM2(1));


fprintf('\n');
fprintf('Frecuencia final:\n');

fprintf('f     = %e Hz\n',f(end));

fprintf('delta = %e m\n',delta(end));

fprintf('LM1   = %e H\n',LM1(end));

fprintf('LM2   = %e H\n',LM2(end));


fprintf('\n');

fprintf('Relacion LM2/LM1 inicial = %e\n', ...
        LM2(1)/LM1(1));

fprintf('Relacion LM2/LM1 final   = %e\n', ...
        LM2(end)/LM1(end));


fprintf('\n');

fprintf('Error relativo inicial = %e %%\n', ...
        error_relativo(1));

fprintf('Error relativo final   = %e %%\n', ...
        error_relativo(end));


fprintf('\n');
fprintf('======================================================\n');


end


% ################################################################
% ################################################################
% SUBFUNCION: SEGUNDA EXPRESION DE L_M
% ################################################################
% ################################################################

function LM = calcula_LM2( ...
    omega,b,sigma,mu0,Nmode)


% ==============================================================
% VECTOR DE INDICES
% ==============================================================

n = (1:Nmode)';


% ==============================================================
% k_n
%
% k_n = (2n-1)*pi/b
% ==============================================================

k = (2*n-1)*pi/b;


% ==============================================================
% kappa_n
%
% -k_n^2 + kappa_n^2 - j*omega*sigma*mu0 = 0
%
% Por tanto:
%
% kappa_n^2 = k_n^2 + j*omega*sigma*mu0
%
% ==============================================================

kappa = sqrt(k.^2 + j*omega*sigma*mu0);


% ==============================================================
% SUMA 1
%
% Suma1 =
%
% SUM {
%
% 32*tanh(kappa_n*b/2)
% ----------------------------------------
% pi^2*(2n-1)^2*kappa_n*b
%
% }
%
% ==============================================================

Suma1 = 0;


for nn = 1:Nmode

    kap = kappa(nn);

    x = kap*b/2;


    % ----------------------------------------------------------
    % tanh(x)
    %
    % Para Re(x) grande:
    %
    % tanh(x) -> 1
    %
    % Evitamos calcular valores innecesarios.
    % ----------------------------------------------------------

    if real(x) > 20

        th = 1;

    else

        th = tanh(x);

    end


    % ----------------------------------------------------------
    % Termino de la sumatoria
    % ----------------------------------------------------------

    termino = ...
        32*th / ...
        (pi^2*(2*nn-1)^2*kap*b);


    Suma1 = Suma1 + termino;

end


% ==============================================================
% SUMA 2
%
% Suma2 =
%
% SUM {
%
% 4*(-1)^(n+1)
% ---------------------------
% (2n-1)*pi
%
% * [1 - 2/cosh(kappa_n*b/2)]
%
% }
%
% ==============================================================

Suma2 = 0;


for nn = 1:Nmode

    kap = kappa(nn);

    x = kap*b/2;


    % ----------------------------------------------------------
    % Calculo estable de:
    %
    % 1 - 2/cosh(x)
    %
    % ----------------------------------------------------------

    if real(x) > 20

        factor = 1;

    else

        factor = 1 - 2/cosh(x);

    end


    % ----------------------------------------------------------
    % Termino de la sumatoria
    % ----------------------------------------------------------

    termino = ...
        4*(-1)^(nn+1)/((2*nn-1)*pi) ...
        * factor;


    Suma2 = Suma2 + termino;

end


% ==============================================================
% EXPRESION FINAL
%
% LM =
%
% Re{
%
% mu0*b*(1-Suma1)
% ----------------
% Suma2
%
% }
%
% ==============================================================

LM = real( ...
    mu0*b*(1-Suma1)/Suma2 );



end

plot(LM2(1,1:100), 'DisplayName', 'LM2(1,1:100)', 'YDataSource', 'LM2(1,1:100)'); figure(gcf)
hold on; plot(LM1(1,1:100), 'DisplayName', 'LM2(1,1:100)', 'YDataSource', 'LM2(1,1:100)'); figure(gcf)
