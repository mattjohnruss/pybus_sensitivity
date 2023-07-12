%% parameter fitting and plot 

% experimental data means
% data.TGFB = 0.22175;
data.ASM = [0.075350741519162; 0.134566093663463; 0.115725487266925;...
    0.092506937555619; 0.104128617681256; 0.102772942714612; 0.078520645547573];  
data.ECM = [0.108408201555051; 0.165391819652303; 0.159142424122995; ...
    0.171816113229631; 0.176121994229251; 0.190120346184978; 0.140735916818001];

% controls
% con.TGFB_34 = 0.185083;
con.ASM_34 = data.ASM(1);
con.ASM_41 = data.ASM(end);
con.ECM_34 = data.ECM(1);
con.ECM_41 = data.ECM(end);

% normalised OVA relative to day 34 controls
% norm.TGFB = data.TGFB./con.TGFB_34;
norm.ASM = data.ASM(1:end-1)./con.ASM_34;
norm.ECM = data.ECM(1:end-1)./con.ECM_34;

% scale time and tspan for simulation plots
K.ts = 1000;
K.phia = 10;
tspan.total = linspace(0,(100+K.ts)*K.phia,(100+K.ts)*K.phia + 1);

% time points from data for fitting function
tspan.TGFB = (34+K.ts)*K.phia;
tspan.ASM_ECM = ([0; 34; 35; 37; 39; 41]+K.ts)*K.phia;
% T = [tspan.TGFB; tspan.ASM_ECM; tspan.ASM_ECM]; 
T = [tspan.ASM_ECM; tspan.ASM_ECM]; 

% data column vector for lsqcurvefit and ODE solver
% D = [norm.TGFB; norm.ASM; norm.ECM];
D = [norm.ASM; norm.ECM];

% intial conditions 
a0 = 0.19;
p0 = 0.005;
c0 = 0.07;
m0 = 0.12;
ICS = [a0 p0 c0 m0];

% initial parameter guesses
%k0 = ones(1,6);
k0 = [0.157, 1.165, 1.384, 2.314, 0.1325, 1.17, 0.3275, 0.1, 2.64, 30];

% k0 = [0.157, 1.165, 1.384, 2.314, 0.1325, 1.17, 0.3275, 0.1, 2.64, 30,...
%     1, 17.3, 0.165, 0.063, 0.001, 1, 0.001, 0.003, 1, 0.114, 3.7, 0.001, 1.8, 0.00138];

% k1 = k(1);%0.157 ka
% k2 = k(2);%1.165 kac
% k3 = k(3);%1.384 nac
% k4 = k(4);%2.314 kb
% k5 = k(5); %0.1325; % kp
% k10 = k(6); %1.17; % kpc
% k11 = k(7); %0.3275; % phic
% k17 = k(8);%0.1 kpm
% ks = k(9);%2.64 ks
% ku = k(10);%30; v

% bounds on parameters
lb = [0.1, 1, 0.1, 1.5, 0.05, 0.5, 0.1, 0.05, 2, 10];
ub = [1, 3, 3, 4, 0.3, 2, 1, 0.15, 6, 100];
% lb = [0.1, 1, 0.1, 1.5, 0.05, 0.5, 0.1, 0.05, 2, 10,...
%     0.5, 15, 0.1, 0.05, 0.0005, 0.5, 0.0005, 0.0001, 0.5, 0.1, 2, 0.0005, 1, 0.0001];
% ub = [1, 3, 3, 4, 0.3, 2, 1, 0.15, 6, 100,...
%     2, 20, 0.2, 0.1, 0.01, 2, 0.01, 0.01, 2, 0.3, 5, 0.01, 3, 0.01];

%% fit parameters

% output optimum parameters with SSE
[k_optim,SSE]=lsqcurvefit(@(k,t)objective(k,tspan,ICS,K),k0,T,D,lb,ub);

% solve ODEs with optimised parameters
sols = solveODE(k_optim,tspan.total,ICS,K);

% sols normalised using steady state
TGFB = sols(:,1)./sols((K.ts/2)*K.phia,1);
P = sols(:,2)./sols((K.ts/2)*K.phia,2);
C = sols(:,3)./sols((K.ts/2)*K.phia,3);
ECM = sols(:,4)./sols((K.ts/2)*K.phia,4);

% pool ASM populations
ASM = (sols(:,2) + sols(:,3))./(sols((K.ts/2)*K.phia,2)+sols((K.ts/2)*K.phia,3));

%% plot solutions

% colours
% colours = [0.6902 0.7686 0.8706; 0.2745 0.5098 0.7059; 0.8471 0.7490 0.8471;...
%     0.8667 0.6275 0.8667; 0.5647 0.9333 0.5647; 0.2353 0.7020 0.4431];
colours = [0.2745 0.5098 0.7059; 0.8471 0.7490 0.8471;...
    0.8667 0.6275 0.8667; 0.5647 0.9333 0.5647; 0.2353 0.7020 0.4431];

% plot data points and optimised parameter set solutions
figure()
% plot(tspan.TGFB,norm.TGFB,'*','markersize',10,'linewidth',3)
hold on
plot(tspan.total,TGFB,'linewidth',3);
plot(tspan.ASM_ECM,norm.ASM,'*','markersize',10,'linewidth',3)
plot(tspan.total,ASM,'linewidth',3);
plot(tspan.ASM_ECM,norm.ECM,'*','markersize',10,'linewidth',3)
plot(tspan.total,ECM,'linewidth',3);
% xlim([0 max(tspan.total)])
% xlim([K.ts*K.phia (45+K.ts)*K.phia])
xlim([9700 10700])
ylim([0.9 2]);
xlabel('Time (days)', 'fontweight', 'bold')
ylabel('Constiuent density', 'fontweight', 'bold')
% l = legend('TGFB data','TGFB sols', 'ASM data', 'ASM sols', 'ECM data',...
%     'ECM sols', 'location','southoutside');
l = legend('TGFB sols', 'ASM data', 'ASM sols', 'ECM data',...
    'ECM sols', 'location','southoutside');
l.NumColumns = 3;  
set(gca, 'fontsize', 26, 'FontName','Calibri')
set(gcf, 'Position',  [100, 100, 600, 600])
colororder(colours);

%% display SSE
disp(SSE)

%% objective function

function output = objective(k,tspan,ICS,K)

% solve ODES with current parameter set, evaluate at mutual timepoints
sols_OVA = solveODE(k,tspan.total,ICS,K);

% pick out timepoints to match data sets normalised using steady state
% indxTGFB = zeros(size(tspan.TGFB));
indxASM_ECM = zeros(size(tspan.ASM_ECM));

% for i = 1:length(tspan.TGFB)
%     indxTGFB(i) = find(tspan.total == tspan.TGFB(i));
% end
for j = 1:length(tspan.ASM_ECM)
    indxASM_ECM(j) = find(tspan.total == tspan.ASM_ECM(j));
end

% sols_TGFB = sols_OVA(indxTGFB,1)./sols_OVA((K.ts-1)*K.phia,1);
sols_ECM = sols_OVA(indxASM_ECM,4)./sols_OVA((K.ts - 1)*K.phia,4);

% pool ASM populations
sols_ASM = (sols_OVA(indxASM_ECM,2) + sols_OVA(indxASM_ECM,3))./(sols_OVA((K.ts - 1)*K.phia,2)+sols_OVA((K.ts - 1)*K.phia,3));
    
% output as colum vector to be consistent with lsqcurvefit input T and data
% output = [sols_TGFB; sols_ASM; sols_ECM];
output = [sols_ASM; sols_ECM];

end

%% ODE function

function outputODE = solveODE(k,tspan,ICS,K)
    
    % solve ODEs 
    [~,outputODE] = ode45(@(t,y)DifEq(t,y,K),tspan,ICS,K);
   
    % ODE system: rate of change in airway wall constituents
    function dy = DifEq(t,y,K) 
        dydt = zeros(4,1); % Size of 3 ODEs
        
        % parameters
        k1 = k(1);%0.157 ka
        k2 = k(2);%1.165 kac
        k3 = k(3);%1.384 nac
        k4 = k(4);%2.314 kb
        k5 = k(5); %0.1325; % kp
        k6 = 1; % pmax k(11);%
        k7 = 17.3; % kap k(12);%
        k8 = 0.165; % nap k(13);%
        k9 = 0.063; % kcp k(14);%
        k10 = k(6); %1.17; % kpc
        k11 = k(7); %0.3275; % phic
        k12 = 0.001; % phicm k(15);%
        k13 = 1; % ncm k(16);%
        k14 = 0.001; % kcm k(17);%
        k15 = 0.003; % kacm k(18);%
        k16 = 1; % nacm k(19);%
        k17 = k(8);%0.1 kpm
        k18 = 0.114; % kapm k(20);%
        k19 = 3.7; % napm k(21);%
        k20 = 0.001; % ke k(22);%
        k21 = 1.8; % ne k(23);%
        k22 = 0.00138; % phim k(24);%
        ks = k(9);%2.64 ks
        ku = k(10);%30; v
        rhom = 1;
        rhoc = 1;
        rhop = 1;
        rhoa = 0.01;
        
        % variables
        a = y(1); % active TGFB
        p = y(2); % proliferative ASM
        c = y(3); % contractile ASM
        m = y(4); % ECM
        
        % stimulus
        % sensitisation 0, 10 then challenges days 17-33
        ti_vec = ([17 18 19 20 21 22 24 26 28 33]+K.ts)*K.phia; 
        k_stim = 0;
            for i = 1:length(ti_vec)
                ti = ti_vec(i);
                k_stim = k_stim + (ks/(2*pi*ku)^0.5)*exp((-(t-ti).^2)/(2*ku^2));
            end
        
        % ODES
        dydt(1) = k1 + (k_stim + k2*a*c/(k3 + a*c))*c*m*rhom/rhoa ...
            - k4*((rhop/rhoc)*p + c)*a - a;
        dydt(2) = k5*(1 - p/k6)*(1 + k7*a*p/(k8 + a*p))*p + (k9*rhoc/rhop)*c - k10*p;
        dydt(3) = (k10*rhoc/rhop)*p - k9*c - (k11*c - k12*m/(k13 + m))*c;
        dydt(4) = (k14 + k15*a*c/(k16 + a*c))*rhoc*c/rhom ...
            + (k17 + k18*a*p/(k19 + a*p))*rhop*p/rhom + k20*a/(k21 + a) - k22*m;
        
        dy = dydt;

    end
        
end