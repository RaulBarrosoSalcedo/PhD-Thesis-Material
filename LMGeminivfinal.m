% =========================================================================
% SCRIPT: Comparación de Expresiones para Lm vs Frecuencia
% Compatible con MATLAB R2009a (Sin funciones locales, un solo archivo)
% =========================================================================

clear all;
close all;
clc;

% --- CONSTANTES DEL PROBLEMA ---
H_0 = 1.0;                  % Amplitud de campo magnético (H_0 = 1)
b = 3e-3;                   % Dimensiones del cubo de aluminio (3 mm)
sigma_c = 3.8e7;            % Conductividad eléctrica del aluminio (S/m)
mu_0 = 4 * pi * 1e-7;       % Permeabilidad magnética del vacío (H/m)

% --- CONFIGURACIÓN DEL RANGO DE FRECUENCIAS ---
f_min = 1.5e3;              % 1.5 KHz
f_max = 1.5e9;              % 1.5 GHz
puntos_f = 200;             % Cantidad de puntos para una curva suave
f = logspace(log10(f_min), log10(f_max), puntos_f);

% --- ALMACENAMIENTO DE RESULTADOS ---
Lm1 = zeros(1, puntos_f);   % Primera expresión (Aproximada)
Lm2 = zeros(1, puntos_f);   % Segunda expresión (Exacta por serie)

% --- NÚMERO DE TÉRMINOS PARA LA CONVERGENCIA ---
N_terms = 50;               % Al ser sumas simples, podemos usar más términos

% --- BUCLE PRINCIPAL DE FRECUENCIAS ---
for idx = 1:puntos_f
    f_val = f(idx);
    omega = 2 * pi * f_val;
    
    % --- 1. CÁLCULO DE LA PRIMERA EXPRESIÓN (APROXIMADA) ---
    delta_c = sqrt(2 / (omega * mu_0 * sigma_c));
    
    % Se asume que b - 2*delta_c no debe ser negativo a bajas frecuencias
    % donde delta_c es mayor que el radio físico del cubo.
    b_efectivo = b - 2 * delta_c;
    if b_efectivo < 0
        b_efectivo = 0; % Límite físico cuando el campo penetra por completo
    end
    Lm1(idx) = (mu_0 * (b_efectivo)^2) / b;
    
    % --- 2. CÁLCULO DE LA SEGUNDA EXPRESIÓN (SERIE ANALÍTICA) ---
    n_arr = 1:N_terms;
    k = (2 * n_arr - 1) * pi / b;
    kappa = sqrt(k.^2 + 1i * omega * sigma_c * mu_0);
    
    % Inicialización de las sumas para la frecuencia actual
    suma_num = 0;
    suma_den = 0;
    
    for n = 1:N_terms
        kn = k(n);
        kpan = kappa(n);
        
        % Cálculo del término tanh(x) de forma segura
        arg_hyper = kpan * b / 2;
        val_tanh = tanh(arg_hyper);
        
        % Tratamiento seguro para el cosh a altas frecuencias para evitar NaN
        % Si la parte real de arg_hyper es muy grande, 2/cosh tiende a 0.
        if real(arg_hyper) > 30
            val_cosh_inv = 0;
        else
            val_cosh_inv = 2 / cosh(arg_hyper);
        end
        
        % Término de la sumatoria del Numerador
        term_num = 32 * val_tanh / ((pi^2) * ((2*n-1)^2) * kpan * b);
        suma_num = suma_num + term_num;
        
        % Término de la sumatoria del Denominador
        term_den = (4 * (-1)^(n+1) / ((2*n-1) * pi)) * (1 - val_cosh_inv);
        suma_den = suma_den + term_den;
    end
    
    % División final y extracción de la parte real según la fórmula
    L_M_compleja = mu_0 * b * (1 - suma_num) / suma_den;
    Lm2(idx) = real(L_M_compleja);
end

% --- GENERACIÓN DE LA GRÁFICA (SEMILOGX) ---
% Como la inductancia cae hacia cero a altas frecuencias, la escala logarítmica 
% en Y podría ocultar los valores negativos o cercanos a cero. Usamos semilogx.
figure('Color', [1 1 1]);
semilogx(f, Lm1, 'r--', 'LineWidth', 2); hold on;
semilogx(f, Lm2, 'b-', 'LineWidth', 2);
grid on;

% Etiquetas y formato
xlabel('Frecuencia (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Inductancia Equivalente L_M (H)', 'FontSize', 12, 'FontWeight', 'bold');
title('Comparación de Expresiones de L_M vs Frecuencia', 'FontSize', 13, 'FontWeight', 'bold');
legend('Aproximación del Cubo (L_{M1})', 'Modelo Analítico de Series (L_{M2})', 'Location', 'SouthWest');

% Ajustes del formato del eje
set(gca, 'FontSize', 11);
