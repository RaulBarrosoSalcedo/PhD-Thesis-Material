%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RM_vs_f.m
%
% Comparacion entre:
%  1) Aproximacion por profundidad de penetracion
%  2) Expresion modal completa
%
% MATLAB 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

%% Constantes

mu0     = 4*pi*1e-7;
eps0    = 8.854187817e-12;

sigma_c = 3.8e7;      % S/m
mu_c    = mu0;
eps_c   = eps0;

H0      = 1;

b       = 3e-3;       % 3 mm

%% Frecuencia

Nf = 150;

f = logspace(log10(1.5e3),log10(1.5e9),Nf);

%% Numero de modos

N = 50;

%% Inicializacion

RM_skin  = zeros(size(f));
RM_modal = zeros(size(f));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bucle principal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ifreq = 1:length(f)

    freq  = f(ifreq);

    omega = 2*pi*freq;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CASO 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    delta_c = sqrt(2/(omega*mu_c*sigma_c));

    RM_skin(ifreq) = 4/(delta_c*sigma_c);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CASO 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    k      = zeros(N,1);
    kappa  = zeros(N,1);
    Aprime = zeros(N,1);

    for m=1:N

        k(m) = (2*m-1)*pi/b;

        kappa(m) = sqrt( k(m)^2 + 1i*omega*mu0*sigma_c );

        Aprime(m)= ...
            4*(-1)^(m+1)*H0/...
            ((2*m-1)*pi*cosh(kappa(m)*b/2));

    end

    %----------------------------------------------------------------------
    % Numerador
    %----------------------------------------------------------------------

    Numerador = 0;

    for m=1:N

        for n=1:N

            km  = k(m);
            kn  = k(n);

            kapm = kappa(m);
            kapn = kappa(n);

            kapmc = conj(kapm);
            kapnc = conj(kapn);

            
            if (m==n)

                ctmp = ...
                    b*( ...
                    (km^2-abs(kapm)^2)* ...
                    sinh(0.5*b*(kapm-conj(kapm))) / ...
                    (kapm-conj(kapm)) ...
                    + ...
                    (km^2+abs(kapm)^2)* ...
                    sinh(0.5*b*(kapm+conj(kapm))) / ...
                    (kapm+conj(kapm)) ) ...
                    ...
                    - ...
                    2*km^2*(conj(kapm)^2)* ...
                    ( cosh(b*kapm/2)/(km^2+kapm^2) )* ...
                    ( cosh(b*kapmc/2)/(km^2+kapm^2) );

            else

                term1 = ...
                    -8*(-1)^(m+n)* ...
                    ( km*(kapnc^2)*kn + kn*(kapm^2)*km )* ...
                    ( cosh(b*kapm/2)/(kn^2+kapm^2) )* ...
                    ( cosh(b*kapnc/2)/(km^2+kapnc^2) );

                term2 = ...
                    2*km*kn* ...
                    ( (-1)^(m-n)/(km-kn) ...
                    -(-1)^(m+n)/(km+kn) )* ...
                    ( ...
                    sinh(b*(kapm-kapnc)/2)/(kapm-kapnc) ...
                    + ...
                    sinh(b*(kapm+kapnc)/2)/(kapm+kapnc) );

                term3 = ...
                    2*kapm*kapnc* ...
                    ( (-1)^(m-n)/(km-kn) ...
                    +(-1)^(m+n)/(km+kn) )* ...
                    ( ...
                    -sinh(b*(kapm-kapnc)/2)/(kapm-kapnc) ...
                    + ...
                    sinh(b*(kapm+kapnc)/2)/(kapm+kapnc) );

                ctmp = term1 + 0*term2 + 0*term3;

            end

            Numerador = Numerador + ...
                Aprime(m)*conj(Aprime(n))*ctmp;

        end
    end

    Numerador = (b/sigma_c)*Numerador;

    %----------------------------------------------------------------------
    % Denominador
    %----------------------------------------------------------------------

    IM = 0;

    for n=1:N

        IM = IM + ...
            4*(-1)^(n+1)*H0*b/((2*n-1)*pi)* ...
            (1 - 2/cosh(kappa(n)*b/2));

    end

    Denominador = abs(IM)^2;

    RM_modal(ifreq) = real(Numerador/Denominador);

end

ft = f';
RM_modalt = RM_modal';
RM_skint = RM_skin';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graficas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

loglog(f,RM_skin,'b','LineWidth',2)
hold on

loglog(f,RM_modal,'r','LineWidth',2)

grid on

xlabel('Frecuencia (Hz)')
ylabel('R_M (\Omega)')

title('R_M vs frecuencia')

legend('Aproximacion profundidad penetracion',...
       'Expresion modal completa')
