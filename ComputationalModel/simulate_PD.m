function [Y,F,S] = simulate_PD(varargin)
clc;
%close all;
%clear all;

% the parameter descriptions of the original model as per Nevado Holgado
% 2010 with additional features
% note we add autoinhibition in cortical layers in accordance with the CMC
% paper of Rosalyn Moran. Additional additions are a model of the GPi. The
% time constants of the GPi are taken from the literature.

Y = [];
F = [];
S = [];

if isempty(varargin)
    p.D        = 1;
    D = p.D;
    p.dTsg     = 6e-3;     %D*6e-3;   % ms
    p.dTgs     = 6e-3;     %D*6e-3;   % ms
    p.dTgg     = D*4e-3;   % ms
    p.dTcs     = D*5.5e-3;
    p.dTsc     = D*20e-3;  %21.5e-3;
    p.dTsgpi   = 6e-3;
    p.dTgpegpi = 6e-3;
    p.dTgpigpi = 4*1e-3;
    
    p.Ts   = 12.8*1e-3;   % ms %12.8
    p.Tg   = 20*1e-3;     % ms %20
    p.Te   = 11e-3;       %10e-3;  % ms
    p.Ti   = 11e-3;       %10e-3;  % ms
    p.Tgpi = 14e-3;       % see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2585404/#!po=67.1053 for reference
    
    p.Ctx  = 172.18;      %27;       % spm/s
    p.Str  = 8.46;        %  2;      % sp/s
    p.Ms   = 300;         % sp/s
    p.Bs   = 10;          %17;     % sp/s
    p.Mg   = 400;         % sp/s
    p.Bg   = 20;          %75;     % sp/s
    p.Mgi  = 300;
    p.Bgi  = 18;
    p.Mi   = 205;         % normal range 200-330 sp/s
    p.Bi   = 9.87;        % normal range 0-20 sp/s
    p.Me   = 75;          % normal range 50-80 sp/s
    p.Be   = 17.85;       % normal range 0-20 sp/s
    
    % for the cortical populations
    p.dTcs  = D*5.5e-3;   % ms - allowed range 1-10 ms
    p.dTcc  = D*6e-3;     %5e-3%1e-3;%4.65e-3;  % 2;     % ms - allowed range 1-10 ms
    
    p.K = 1;
    % disease
    p.Wsg  = 4.87;  %19 + p.K*(20-19);
    p.Wgs  = 1.33;  %1.12 + p.K*(10.7-1.12);
    p.Wgg  = 0.53;  %0.53;%6.6 + p.K*(12.3 - 6.6);
    p.Wcs  = 7.98;     %2.42 + p.K*(9.2 - 2.42);
    
    % cortical parameters
    p.Wsc     = 0; %0.4;%0.4;%0.4;%0.4%8.93;
    p.Wcc     = 4;
    p.Wsgpi   = 3.7;
    p.Wgpegpi = 1;
    
    % stochastic component - for integration and for summation
    p.S       = 40;
    p.SN      = 0;
    
    % generate plot
    p.Plot    = 1;
    
    % plot phase portrait
    p.PhaseP  = 0;
    
    % integration method
    p.Method  = 'S';
    
    % Time duration
    p.length  = 5;
else
    p = varargin{1};
end


% integrate model via deterministic or stochastic approach
IC              = zeros(5,1);

if     strcmp(p.Method,'D')
    %ops_DDE         = ddeset('RelTol',1e-5,'AbsTol',1e-5,'MaxStep',5*1e-3);
    ops_DDE         = ddeset('MaxStep',1e-3);
    sol             = dde23(@ddefun, [p.dTsg,p.dTgs,p.dTgg, p.dTcs, p.dTsc, p.dTcc, p.dTsgpi, p.dTgpegpi, p.dTgpigpi]', IC, [0 p.length], ops_DDE,p);
elseif strcmp(p.Method,'S')
    sol = Euler_Maruyama(p,IC,p.length,1e-3);
end
sol.y           = sol.y + p.SN*randn(size(sol.y,1),size(sol.y,2));

if p.Plot
    if p.PhaseP
        im = 3;
    else
        im = 2;
    end
    
    figure;subplot(im,1,1);
    plot(sol.x,sol.y(1,:),'b','LineWidth',2);hold on;
    plot(sol.x,sol.y(2,:),'g','LineWidth',2);hold on;
    plot(sol.x,sol.y(3,:),'r','LineWidth',2);hold on;
    plot(sol.x,sol.y(4,:),'k','LineWidth',2);hold on;
    plot(sol.x,sol.y(5,:),'m','LineWidth',2);hold on;
    axis([0 1 0 max(sol.y(:))+20]);
    xlabel('Time (s)');ylabel('Firing Rate (Hz)');
    set(gca,'FontSize',12);
    %title('DDE 23');
    
    
    [S,~,f] = ft_specest_mtmfft(sol.y(:,sol.x>0.5),sol.x(sol.x>0.5),'taper','dpss','tapsmofrq',1,'freqoi',[0:1: 50]);
    s = squeeze(mean(abs(S),1));
    
    subplot(im,1,2);
    plot(f,s(1,:),'b','LineWidth',2);hold on;
    plot(f,s(2,:),'g','LineWidth',2);hold on;
    plot(f,s(3,:),'r','LineWidth',2);hold on;
    plot(f,s(4,:),'k','LineWidth',2);hold on;
    plot(f,s(5,:),'m','LineWidth',2);hold on,
    xlabel('Frequency (Hz)');ylabel('Power (AU)');
    h = legend({'STN','GPe','GPi','Cortex E','Cortex I'});
    set(h,'box','off','LineWidth',2);
    if isfield(p,'axis')
        axis([0 p.axis(1) 0 p.axis(2)]);
    end
    set(gca,'FontSize',12);
    
    if p.PhaseP
        subplot(im,1,3)
        scatter(sol.y(1,sol.x>1.5),sol.y(2,sol.x>1.5),10,sol.y(4,sol.x>1.5),'filled')
        xlabel('STN rate (Hz)');
        ylabel('GPe rate (Hz)');
        box on;
        set(gca,'FontSize',12);
        try
        ylim(p.ylim)
        end
    end
end
Y.y = sol.y;
Y.x = sol.x;
if p.Plot
    F.f = f;
    F.s = s;
end


end


function y = sigmoid (x,max,min)
y = max./(1+(exp(-4*x*inv(max))*(max-min)/min));
end

function dydt = ddefun(t,Y,Z,p)
% the integral function
% Y(1) = STN
% Y(2) = GP
% Y(3) = GPi
% Y(4) = Cortex - E
% Y(5) = Cortex - I

% Z(:,1) = dTsg
% Z(:,2) = dTgs
% Z(:,3) = dTgg
% Z(:,4) = dTcs
% Z(:,5) = dTsc  % Note that this is now an inhbitory effect of the GPi
% Z(:,6) = dTcc
% Z(:,7) = dTsGPi
% Z(:,8) = dTGPeGPi
% Z(:,9) = 

% remember that Z is coded Z(i,j) such that i is the giving population and
% j is the delay

% [p.dTsg,p.dTgs,p.dTgg, p.dTcs, p.dTsc, p.dTcc, p.dTsgpi, p.dTgpegpi, p.dTgpigpi]

dydt = zeros(5,1);
S       =   p.S;
dydt(1) =   (sigmoid(-p.Wgs*Z(2,2)   + p.Wcs*Z(4,4),p.Ms,p.Bs)                                 - Y(1) + S*randn)/p.Ts;   % STN
dydt(2) =   (sigmoid( p.Wsg*Z(1,1)   - p.Wgg*Z(2,3) - p.Str,p.Mg,p.Bg)                         - Y(2) + S*randn)/p.Tg;   % GPe
dydt(3) =   (sigmoid( p.Wsgpi*Z(1,7) - p.Wgpegpi*Z(2,8) - p.Wgg*Z(3,9) -p.Str,p.Mgi,p.Bgi)     - Y(3) + S*randn)/p.Tgpi; % GPi          
dydt(4) =   (sigmoid(-p.Wsc*Z(3,5)   - p.Wcc*Z(5,6) + p.Ctx,p.Me,p.Be)                         - Y(4) + S*randn)/p.Te;   % Cortex E
dydt(5) =   (sigmoid(p.Wcc*Z(4,6)    - p.Wcc/1*Z(5,6),p.Mi,p.Bi)                               - Y(5) + S*randn)/p.Ti;   % Cortex I

end

function sol = Euler_Maruyama(p,IC,T,dt)

% solves a bespoke SDDE of the following form using the Euler-Maruyama
% method which is a truncation of the Ito-Taylor expansion. 

% X(j) = X(j-1) + f(X(j-1),X(j-Nd)*dT + K*(W(j) - W(j-1)

% https://epubs.siam.org/doi/pdf/10.1137/S0036144500378302
% https://www.sciencedirect.com/science/article/pii/S037704270500124X#bib1

% y(1) = STN
% y(2) = GP
% y(3) = GPi
% y(4) = Cortex - E
% y(5) = Cortex - I


Nt         = T/dt;                         % discretise integration steps
max_delays = max([p.dTsg,p.dTgs,p.dTgg, p.dTcs, p.dTsc, p.dTcc, p.dTsgpi, p.dTgpegpi, p.dTgpigpi]);

% express delays in terms of time steps
dTsg       = floor(p.dTsg/dt);
dTgs       = floor(p.dTgs/dt);
dTgg       = floor(p.dTgg/dt);
dTcs       = floor(p.dTcs/dt);
dTsc       = floor(p.dTsc/dt);
dTcc       = floor(p.dTcc/dt);
dTsgpi     = floor(p.dTsgpi/dt);
dTgpegpi   = floor(p.dTgpegpi/dt);
dTgpigpi   = floor(p.dTgpigpi/dt);

Ndt        = max_delays/dt;

y          = zeros(numel(IC),Nt);          % set the output matrix
dW         = sqrt(dt).*randn(size(y));     % brownian motion 
W          = cumsum(dW,2);                 % discretised brownian path
yA         = [zeros(numel(IC),Ndt) y];     % augment for delays

for i = Ndt+2:size(yA,2)
    
    % deterministic portion  - equation of motion multiplied by step size dfdt*dt
    dfdt_dt(1) = dt * (sigmoid(-p.Wgs      *yA(2,i-dTgs)       + p.Wcs     * yA(4,i-dTcs),p.Ms,p.Bs)    - yA(1,i-1))/p.Ts;
    dfdt_dt(2) = dt * (sigmoid( p.Wsg      *yA(1,i-dTsg)       - p.Wgg     * yA(2,i-dTgg)              - p.Str,p.Mg,p.Bg)    - yA(2,i-1))/p.Tg;
    dfdt_dt(3) = dt * (sigmoid( p.Wsgpi    *yA(1,i-dTsgpi)     - p.Wgpegpi * yA(2,i-dTgpegpi)          - p.Wgg*yA(3,i-dTgpigpi) -p.Str,p.Mgi,p.Bgi) - yA(3,i-1))/p.Tgpi;
    dfdt_dt(4) = dt * (sigmoid(-p.Wsc      *yA(3,i-dTsc)       - p.Wcc     * yA(5,i-dTcc)              + p.Ctx                       ,p.Me,p.Be)    - yA(4,i-1))/p.Te;
    dfdt_dt(5) = dt * (sigmoid( p.Wcc      *yA(4,i-dTcc)       - p.Wcc     * yA(5,i-dTcc)              ,p.Mi,p.Be)    - yA(5,i-1))/p.Ti; 
    
    if 1
        % stochastic portion - Brownian motion (signal dependence)
        g_dw(1)    = yA(1,i-1) * p.S * dW(1,i-Ndt-1);
        g_dw(2)    = yA(2,i-1) * p.S * dW(2,i-Ndt-1);
        g_dw(3)    = yA(3,i-1) * p.S * dW(3,i-Ndt-1);
        g_dw(4)    = yA(4,i-1) * p.S * dW(4,i-Ndt-1);
        g_dw(5)    = yA(5,i-1) * p.S * dW(5,i-Ndt-1);
    else
        g_dw(1)    = p.S * dW(1,i-Ndt-1);
        g_dw(2)    = p.S * dW(2,i-Ndt-1);
        g_dw(3)    = p.S * dW(3,i-Ndt-1);
        g_dw(4)    = p.S * dW(4,i-Ndt-1);
        g_dw(5)    = p.S * dW(5,i-Ndt-1);
    end
    % EM update
    yA(1,i)    = yA(1,i-1) + dfdt_dt(1) + g_dw(1);
    yA(2,i)    = yA(2,i-1) + dfdt_dt(2) + g_dw(2);
    yA(3,i)    = yA(3,i-1) + dfdt_dt(3) + g_dw(3);
    yA(4,i)    = yA(4,i-1) + dfdt_dt(4) + g_dw(4);
    yA(5,i)    = yA(5,i-1) + dfdt_dt(5) + g_dw(5);
    
end
yA = yA(:,Ndt+1:end);
sol.y = yA;
sol.x = 0:dt:T-dt;


end
