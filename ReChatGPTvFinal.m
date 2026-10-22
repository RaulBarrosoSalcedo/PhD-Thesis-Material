clear all;
close all;
clc;

%=========================================================
% CONSTANTES
%=========================================================
mu0     = 4*pi*1e-7;
eps0    = 8.854187817e-12;

sigma_c = 3.8e7;      % S/m
b       = 3e-3;       % 3 mm
E0      = 1;

% Numero de terminos de las series
N = 50;

%=========================================================
% VECTOR DE FRECUENCIAS
%=========================================================
f = logspace(log10(1.5e1),log10(2e9),50);

Re_aprox = zeros(size(f));
Re_serie = zeros(size(f));

%=========================================================
% BARRIDO EN FRECUENCIA
%=========================================================
for kf = 1:length(f)

    omega = 2*pi*f(kf);

    %-----------------------------------------------------
    % CASO 1: FORMULA APROXIMADA
    %-----------------------------------------------------
%    Re_aprox(kf) = sqrt(mu0*omega/(32*sigma_c));
    Re_aprox(kf) = sqrt(pi*mu0*f(kf)/(16*sigma_c));
    %-----------------------------------------------------
    % CASO 2: SERIE COMPLETA
    %-----------------------------------------------------

    numerador = 0;
    sumaDen   = 0;

    kvec      = zeros(N,1);
    kappavec  = zeros(N,1);
    Avec      = zeros(N,1);

    %=====================================================
    % COEFICIENTES km, kappam, Am
    %=====================================================
    for m = 1:N

        km = (2*m-1)*pi/b;

        % ecuacion:
        % -km^2 + kappa^2 - j*w*sigma*mu0 = 0
        %
        % kappa^2 = km^2 + j*w*sigma*mu0

        kappa = sqrt(km^2 + 1i*omega*sigma_c*mu0);

        kvec(m)     = km;
        kappavec(m) = kappa;

        % Se asume Ec = E0 = 1
        Avec(m) = 4*(-1)^(m+1)*E0 / ...
                 ((2*m-1)*pi*cosh(kappa*b/2));

        sumaDen = sumaDen + ...
            tanh(kappa*b/2)/((2*m-1)*kappa*km);

    end

    %=====================================================
    % SUMA DOBLE
    %=====================================================
    for m = 1:N
m = m
kf = kf
        km      = kvec(m);
        kappam  = kappavec(m);
        Am      = Avec(m);

        for n = 1:N

            kn      = kvec(n);
            kappan  = kappavec(n);
            An      = Avec(n);

            kappanC = conj(kappan);

            %---------------------------------------------
            % CALCULO DE c_mn
            %---------------------------------------------
            if m == n

                d1 = kappam - kappanC;

                if abs(d1) < 1e-12
                    T1a = b/2;
                else
                    T1a = sinh(b*d1/2)/d1;
                end

                T1b = sinh(b*(kappam+kappanC)/2) / ...
                      (kappam+kappanC);

                T1 = b*(T1a + T1b);

                T2 = 8*km^2 * ...
                    ( cosh(kappam*b/2) / ...
                      (km^2+kappam^2) ) * ...
                    ( cosh(kappanC*b/2) / ...
                      (km^2+kappanC^2) );

                cmn = T1 + T2;

            else % Ojo porque en la siguiente línea està el error que detecta la AI

                T1 = 0*2*( ...
                    (-1)^(m-n)/(km-kn) + ...
                    (-1)^(m+n)/(km+kn) );

                T2 = ...
                    sinh(b*(kappam-kappanC)/2) / ...
                    (kappam-kappanC) ...
                    + ...
                    sinh(b*(kappam+kappanC)/2) / ...
                    (kappam+kappanC);

                T3 = 8*(-1)^(m+n)*km*kn * ...
                    ( cosh(kappam*b/2) / ...
                      (kn^2+kappam^2) ) * ...
                    ( cosh(kappanC*b/2) / ...
                      (km^2+kappanC^2) );

                cmn = T1*T2 + T3;

            end

            numerador = numerador + Am*conj(An)*cmn;
            cnmAnalitico(m,n) = cmn;

%           F = @(y,z)((cosh(kappam*y)*cos(km*z)+cosh(kappam*z)*cos(km*y))*(cosh(kappanC*y)*cos(kn*z)+cosh(kappanC*z)*cos(kn*y)));
            %F1 = @(y,z)( cosh(kappam*y).*cos(km*z)*cosh(kappam*z).*cos(km*y) );
            %F2 = @(y,z)( cosh(kappam*y).*cos(km*z)*cosh(kappanC*z).*cos(kn*y) );            
            %F3 = @(y,z)( cos(km*y).*cosh(kappam*z)*cos(kn*z).*cosh(kappanC*y) );             
            %F4 = @(y,z)( cos(km*y).*cosh(kappam*z)*cosh(kappanC*z).*cos(kn*y) );
            
            %cmn2 = dblquad(F1,-b/2,b/2,-b/2,b/2);
            %cmn3 = dblquad(F2,-b/2,b/2,-b/2,b/2);
            %cmn4 = dblquad(F3,-b/2,b/2,-b/2,b/2);
            %cmn5 = dblquad(F4,-b/2,b/2,-b/2,b/2);
            %cnmNumerico(m,n) = cmn2+cmn3+cmn4+cmn5;            

            

        end
    end

    numerador = (b*(pi^2))/(1024*sigma_c)*numerador;

    denominador = abs(E0*sumaDen)^2;

    Re_serie(kf) = real(numerador/denominador);

end

ft = f';
Re_seriet = Re_serie';
Re_aproxt = Re_aprox';

%=========================================================
% GRAFICA
%=========================================================
figure;

loglog(f,Re_aprox,'b','LineWidth',2);
hold on;

loglog(f,abs(Re_serie),'r--','LineWidth',2);

grid on;

xlabel('Frecuencia (Hz)');
ylabel('R_e (\Omega)');

title('R_e vs Frecuencia');

legend('Aproximacion','Serie completa');

set(gca,'FontSize',12);