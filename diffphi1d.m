%============================================================
% RISOLUZIONE EQUAZIONE DI DIFFUSIONE NON LINEARE (1D)
% ============================================================

clear all;
global c0_global phi0_global

% 1. DEFINIZIONI

h = 0.95e-3;                % altezza del campione (m)
z = linspace(0, h, 1000);   % griglia spaziale
t = linspace(0, 2000, 400); % tempi (secondi)
c0 = 0.4;
c0_global = 0.4;            %concentrazione iniziale
phi0_global = 0.01;


% 2. RISOLUZIONE

m = 0;                      % geometria cartesiana (1D)
sol = pdepe(m, @pdefun, @icfun, @bcfun, z, t);
c = sol(:,:,1);             % glicerolo
phi = sol(:,:,2);           % colloide



% 3. PLOT RISULTATI

times_to_plot = [1 4 25 194];    %Dati vari

rho0 = 997;
b = 230;
rhoL = 2200;  
c0 = 0.4;
phi0 = 0.01;

z_norm = z / h;



figure; hold on;                   %Plot glicerolo (t=0,25,125,970 s)

for i = times_to_plot
    plot(z_norm, c(i,:) / c0, 'LineWidth', 2);
end

xlabel('z/h');
ylabel('c / c_0');
title('Glicerolo');
grid on;



figure; hold on;                   %Plot Colloide (t=0,25,125,970 s)

for i = times_to_plot
    plot(z_norm, phi(i,:) / phi0, 'LineWidth', 2);
end

xlabel('z/h');
ylabel('\phi / \phi_0');
title('Colloide');
grid on;



figure; hold on;                   %Plot densità1 (t=0,25,125,970 s)

for i = times_to_plot
    
    c_i = c(i,:);
    phi_i = phi(i,:);
   
    rho = rho0 + (rhoL - rho0).*phi_i + b*(1 - phi_i).*c_i;
    
    plot(z_norm, rho/rho0 - 1, 'LineWidth', 2);
end

xlabel('z/h');
ylabel('\rho/\rho_0 - 1');
title('Profilo di densità');
grid on;



times_to_plot = [127 140 194];

figure; hold on;                   % Plot densità2 (t=633, 700, 970 s) 

for i = times_to_plot
    
    c_i = c(i,:);
    phi_i = phi(i,:);
   
    rho = rho0 + (rhoL - rho0).*phi_i + b*(1 - phi_i).*c_i;
    y = rho/rho0 - 1;

    % --- TROVA PICCO ---
    [pks, locs] = findpeaks(y, z_norm);

    if ~isempty(pks)
        [peakValue, idx] = max(pks);     % picco principale
        peakLocation = locs(idx);

        % segna il picco sul grafico
        plot(peakLocation, peakValue, 'ro', 'MarkerFaceColor','r');
    end

    % plot curva
    plot(z_norm, y, 'LineWidth', 2);
end

xlabel('z/h');
ylabel('\rho/\rho_0 - 1');
title('Profilo di densità');
grid on;


%============================================================
% FIGURA 5 - RAYLEIGH (VERSIONE CORRETTA)
%============================================================

% --- STOKES-EINSTEIN ---
kB = 1.380649e-23;      
T  = 298.15;            % 25 gradi + 273.15 
r  = 11e-9;              
g = 9.81;

Ra_crit = 27*pi^4/4;


c_values = [0.20 0.28 0.40];  
phi_values = [0.01 0.02 0.04 0.08];

figure;

for p = 1:length(phi_values)

    phi0_global = phi_values(p);

    subplot(2,2,p)
    hold on

    for k = 1:length(c_values)

        c0_global = c_values(k);

        % RISOLUZIONE PDE
        sol = pdepe(m, @pdefun, @icfun, @bcfun, z, t);
        c = sol(:,:,1);
        phi = sol(:,:,2);

        Ra = zeros(length(t),1);

        for i = 1:length(t)

            c_i = c(i,:);
            phi_i = phi(i,:);

            % DENSITÀ
            rho = rho0 + (rhoL - rho0).*phi_i + b*(1 - phi_i).*c_i;

            % ===== REGIONE INSTABILE =====
            drhodz = gradient(rho, z);

            unstable = drhodz > 0;  % o <0 a seconda del caso

            if any(unstable)
              z_unst = z(unstable);
              rho_unst = rho(unstable);

              Dz = max(z_unst) - min(z_unst);
              Dr = max(rho_unst) - min(rho_unst);
            else
              Dz = 0;
              Dr = 0;
            end

            % Rayleigh (stokes-einstein)
            if Dz > 0 && Dr > 0
                Ra(i) = g * Dr * Dz^3 * (6*pi*r) / (kB*T);
            else
               Ra(i) = 0;
            end
        end

        % NORMALIZZAZIONE
        Ra_norm = Ra / Ra_crit;

        % PLOT
        semilogy(t, Ra_norm, 'LineWidth', 2, 'DisplayName', ['c_0 = ' num2str(c0_global)]);
    end

    yline(1,'--k')

    xlabel('t (s)')
    ylabel('Ra_S / Ra_S^*')
    title(['\phi_0 = ' num2str(phi0_global)])
    ylim([0 1.5])
    grid on

    if p == 1
        legend show
    end
end



% ============================================================
% 4. FUNZIONI
% ============================================================

function [c,f,s] = pdefun(z,t,u,dudz)

    c = zeros(2,1);     % PREALLOCAZIONE (FONDAMENTALE ALICE)
    f = zeros(2,1);
    s = zeros(2,1);
    
    cgli = u(1);        %variabili
    phi = u(2);

    dcglidz = dudz(1);
    dphidz = dudz(2);

    
    rho0 = 997;         %parametri
    b = 230;
    D0 = 1.025e-9;
    eps = -1.308e-9;
    a = 125;
    Kb=1.38e-23;
    
    T=25+273.15; %K
    Tc=T-273.15;  %C
    chi=0.705;
    csi=2;
    gamma=1-u(1)+(chi*csi*u(1).*(1-u(1)))./(chi*u(1)+csi*(1-u(1)));

    eta_g=12100*exp((-1233+Tc)*Tc/(9900+70*Tc));
    eta_w=1.79*exp((-1230-Tc)*Tc/(36100+360*Tc));%0.001 Ns/m2
    eta=eta_w.^gamma*eta_g.^(1-gamma)*0.001;% Ns/m2

    r=11e-9;%m
    Dc=Kb*T./(6*pi*eta*r);%cm2/s


    
    Ds = D0 + eps*cgli;  %coefficienti
    rho = rho0 + b*cgli;

    
    c(1) = rho;              %equaz. glicerolo
    f(1) = Ds*(cgli*b*dcglidz+rho*dcglidz);   %Ds * rho* dcglidz;
    s(1) = 0;

    
    c(2) = 1;                %equaz. colloide
    f(2) = Dc * dphidz + (a/rho0) * Dc * phi * (cgli*b*dcglidz+rho*dcglidz); %Dc * dphidz + (a/rho0) * Dc * phi * rho*dcglidz;
    s(2) = 0;

end




function u0 = icfun(z)                      %condiz iniziali

    global c0_global phi0_global

    h = 0.95e-3;
    
    if z < h/2
        c0 = c0_global;
        phi0 = phi0_global;
    else
        c0 = 0;
        phi0 = 0;
    end

    u0 = [c0; phi0];   % <-- COLONNA (fondamentale!)
end





function [pl,ql,pr,qr] = bcfun(zl,ul,zr,ur,t) % condiz al contorno (no flusso alle pareti)

    pl = [0; 0];
    ql = [1; 1];

    pr = [0; 0];
    qr = [1; 1];
end






