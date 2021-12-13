%% File that runs model Based on size non specific liebig full
close all
clear
res       = 100;                     % model resolution
 
%% Environmental timeseries
L  = 70; % W/m2 at 10 meters depth
N  = 40; % µg N/l
Si = 30; % µg Si/l [0.00494 0.00451 0.00367 0.00238 0.00083 0.00044 0.00040 0.00073 0.00167 0.00290 0.00403 0.00427]; % mol Si/m3
G  = 15; % µg C/l
 
% trait space
c         = logspace(-7,1,res); % carbon size range
c         = repmat(c,res,1);    % matrix of radius (v,c)
vfrac     = linspace(0,1,res);  % range of vacuole fraction
vfrac     = repmat(vfrac',[1 res]);
 
%
ccb_g     = 0.125*10^-6;        % µg C µm-3 (ONLY BIOMASS)
Vb        = c./ccb_g;           % biovolume
r0        = (Vb./(4/3*pi)).^(1/3); % equvalent spherical radius of biovolume (µm)
m         = 8*10^-3;            % membrane thickness (µm)
% equvalent spherical radius (including outer membrane, but no shell)
ESR       = (m^2./(vfrac - 1) + m^2./(vfrac - 1).^2)./((((3*m^3)./(2*(vfrac - 1).^2) + m^3./(vfrac - 1).^3 + (m^3 + r0.^3)./(2*vfrac - 2)).^2 - (m^2./(vfrac - 1) + m^2./(vfrac - 1).^2).^3).^(1/2) - (3*m^3)./(2*(vfrac - 1).^2) - m^3./(vfrac - 1).^3 - (m^3 + r0.^3)./(2*vfrac - 2)).^(1/3) - m./(vfrac - 1) + ((((3*m^3)./(2*(vfrac - 1).^2) + m^3./(vfrac - 1).^3 + (m^3 + r0.^3)./(2*(vfrac - 1))).^2 - (m^2./(vfrac - 1) + m^2./(vfrac - 1).^2).^3).^(1/2) - m^3./(vfrac - 1).^3 - (m^3 + r0.^3)./(2*(vfrac - 1)) - (3*m^3)./(2*(vfrac - 1).^2)).^(1/3);
ESR(1,:)  = r0(1,:)+m;
V0        = 4/3*pi*ESR.^3;      % total volume (outer membrane, no shell)
%V0(1,:)   = 0;
Vmem      = 4/3*pi*(ESR).^3-4/3*pi*(ESR-m).^3; % outer (plasmalemma) membrane volume
Vv        = V0 - Vb - Vmem;     % Vacuole volume (inkluding inner totnoplast membrane)
Vv(1,:)   = 0;
rv        = (Vv./(4/3*pi)).^(1/3); % radius of vacuole
t         = (10.^(0.24*log10(V0)+0.21))*10^-3; % Shell thickness (µm)
Vs        = (4/3*pi*(ESR+t).^3)-4/3*pi*(ESR).^3;% shell volume (µm3)
RNC       = 1/7;                % (mol N (mol C)-1 )
RCN       = 7;                  % (mol C (mol N)-1 )
RCN_g     = RCN*12/14;          % (g C (g N)-1 )
RNC_g     = RNC*14/12;          % (g N (g C)-1 )
Ncyto     = c*RNC_g;            % µg N/cell cytoplasm (no tonoplast)
Vtono     = 4/3*pi*(rv).^3-4/3*pi*(rv-m).^3; % volume of (inner) tonoplast membrane
 
% densities
rhoW      = 1027;               % Density of seawater (kg/m3)
rhoB      = 1200;               % Density of biomass (kg/m3) SA
rhoS      = 2600;               % Density of shell (kg/m3) (max) SA
rhoV      = 1010;               % Density of vacuole (kg/m3)
% diatom overall density (r,v)
rhoD      = rhoS*Vs./(V0+Vs)+rhoV*Vv./(V0+Vs)+rhoB*(Vb+Vtono+Vmem)./(V0+Vs);
 
rhoCtono  = 600*10^9/10^18;     % carbon density of tonoplast µg C/µm3
rhoCmoltono = rhoCtono/12;      % µmol/µm3
ctono     = Vtono*rhoCtono;     % carbon in tonoplast µg C
cmem      = Vmem*rhoCtono;      % carbon in plasmalemma µg C
Cmoltono  = Vtono*rhoCmoltono;  % carbon in tonoplast µmol C
Cmolmem   = Vmem*rhoCmoltono;   % carbon in plasmalemma µmol C
CNtono    = 48;                 % CN ratio in membrane material (mol/mol)
Nmoltono  = Cmoltono./CNtono;   % Nitrogen in tonoplast µmol N
Nmolmem   = Cmolmem./CNtono;    % Nitrogen in plasmalemma µmol N
Ntono     = Nmoltono*14;        % Nitrogen in tonoplast µg N
Nmem      = Nmolmem*14;         % Nitrogen in plasmalemma µg N
Ntot      = Ncyto+Ntono+Nmem;   % Total nitrogen in cell (µg N cell-1)
Ctot      = c+ctono+cmem;       % Total carbon in cell (µg N cell-1)
membranecarb = (ctono+cmem)./Ctot; % fraction of carbon in membrane
CNtot     = Ctot./Ntot;         % effective CN ratio µg/µg
CNtot_mol = CNtot*14/12;        % mol/mol
ccb       = ccb_g/12*10^-6;     % (mol C µm-3) (ONLY BIOMASS)
ccbtot_mol= (Ctot)./(Vb+Vtono+Vmem)/12*10^-6; % mol C (µm3 biomass and membrane)-1
gammaSi   = rhoS*10^3/60.08*10^-18;% mole Si (µm3 shell)-1
K         = gammaSi./(ccbtot_mol); % mol Si per µm3 shell per mol C per µm3 biomass (GENOVEVEJ)
RSiC      = K.*Vs./(Vb+Vtono+Vmem);% mol Si (mol C)-1
RCSi      = 1./RSiC;            % mol C (mol Si)-1
RCSi_g    = RCSi*28/12;         % g C (g Si)-1
jmax      = 1;                  % maximum synthesis rate (day)-1
Jmax      = jmax*(c);           % maximum synthesis (µg C (day)-1)
 
 
% grazing
fG        = 0.09;
m0        = 0.05;
 
% affinities
 
r_meter   = (ESR)*10^-6;        % radius in meters
DN        = 1.88*10^-5*10^-4*60*60*24;% Nitrogen diffusion (m2/day) (1.62*10^-4)
DSi       = 1.88*10^-5*10^-4*60*60*24;% Silicon diffusion (m2/day) (1.62*10^-4)
cL        = 0.01*(ccb_g)^(2/3); % µg C /day /(W/m2)/µm3^(2/3))
AL        = cL.*(V0).^(2/3);    % µg C/day/(W/m2) (r,v)
ALspec    = AL./c;              % specific affiliny for light (1/day/(W/m2))
AN        = 4*pi*DN*(r_meter)*10^3;% Nitrogen affinity (l/day)
ASi       = 4*pi*DSi*(r_meter)*10^3;% Nitrogen affinity (l/day)
 
resp0     = 0.05;               % respiration rate (day-1)
Jr0       = resp0.*(Ctot);      % µg C/day
 
% costs
betaL     = 0.35;               % Cost of photosynthetic uptake and synthesis (µg C (µg C)-1)
betaN     = 3;                  % Cost of nitrogen uptake and synthesis (µg C (µg N)-1) SAT NED FRA 3.5
betaSi    = 1.05;               % Cost of silicification (µg C (µg Si)-1)
 
 
% sinking mortality
grav     = 9.80665;             % tyngdeaccelerationen (m s-2)
nu       = 1.07e-3;             % dynamic viscosity of water (Pa s)
H        = 20;                  % Depth of euphotic zone.
% sinking speed (m s-1)
w        = 2*grav*(10^-6*(V0./(4/3*pi)).^(1/3)).^2/(9*nu).*(rhoD-rhoW);
ms       = w./H.*60*60*24;      % sinkin mortality (day^-1)
ms(ms<0) = 0;
 
%% Find optimal growth in trait space
for n    = 1:length(L)
    
    % Division rate
    JN(:,:,n)       = AN*N(n);  % µg C/day
    JSi(:,:,n)      = ASi*Si(n);% µg Si/day
    JL(:,:,n)       = AL*L(n);  % µg N/day
    
    % pC = 1
    pN1(:,:,n)      = max(0,min(-(Jr0.*RCSi_g - JL(:,:,n).*RCSi_g + JL(:,:,n).*RCSi_g*betaL)./(JN(:,:,n).*(CNtot*betaSi + RCSi_g*betaN + CNtot.*RCSi_g)),1));
    pSi1(:,:,n)     = max(0,min(-(CNtot.*Jr0 - CNtot.*JL(:,:,n) + CNtot.*JL(:,:,n)*betaL)./(JSi(:,:,n).*(CNtot*betaSi + RCSi_g*betaN + CNtot.*RCSi_g)),1));    
    % pN = 1
    pSi2(:,:,n)     = max(0,min((CNtot.*JN(:,:,n))./(JSi(:,:,n).*RCSi_g),1));
    pC2(:,:,n)      = max(0,(Jr0.*RCSi_g + CNtot.*JN(:,:,n).*RCSi_g + CNtot.*JN(:,:,n)*betaSi + JN(:,:,n).*RCSi_g*betaN)./(RCSi_g.*(JL(:,:,n) - JL(:,:,n)*betaL)));
    % pSi = 1
    pN3(:,:,n)      = max(0,min((JSi(:,:,n).*RCSi_g)./(JN(:,:,n).*CNtot),1));
    pC3(:,:,n)      = max(0,(CNtot.*Jr0 + CNtot.*JSi(:,:,n).*RCSi_g + CNtot.*JSi(:,:,n)*betaSi + JSi(:,:,n).*RCSi_g*betaN)./(CNtot.*(JL(:,:,n) - JL(:,:,n).*betaL)));
 
    % loop through all carbon (v) and vacuole (rad) sizes
    for rad=1:res
        for v=1:res
            if pN1(rad,v,n)<1 && pSi1(rad,v,n)<1
                pN(rad,v,n)     = pN1(rad,v,n);
                pSi(rad,v,n)    = pSi1(rad,v,n);
                pC(rad,v,n)     = 1;
                
                limmat(rad,v,n) =1;
 
            elseif pSi2(rad,v,n)<1 && pC2(rad,v,n)<1
                pSi(rad,v,n)    = pSi2(rad,v,n);
                pN(rad,v,n)     = 1;
                pC(rad,v,n)     = pC2(rad,v,n);
                
                limmat(rad,v,n) =2;
 
            elseif pN3(rad,v,n)<1 && pC3(rad,v,n)<1
                pSi(rad,v,n)    = 1;
                pN(rad,v,n)     = pN3(rad,v,n);
                pC(rad,v,n)     = pC3(rad,v,n);
 
                limmat(rad,v,n) =3;
 
            end
            pC(res,:,n)         = NaN;
            pSi(res,:,n)        = NaN;
            pN(res,:,n)         = NaN;
            Jlieb(rad,v,n)      = max(min([JL(rad,v,n)-Jr0(rad,v)-betaL*JL(rad,v,n)-pN(rad,v,n)*betaN*JN(rad,v,n)-pSi(rad,v,n)*betaSi*JSi(rad,v,n) pN(rad,v,n)*JN(rad,v,n)*CNtot(rad,v) pSi(rad,v,n)*JSi(rad,v,n).*RCSi_g(rad,v)]),0);
 
            if JL(rad,v,n)-Jr0(rad,v)-betaL*JL(rad,v,n)-pN(rad,v,n)*betaN*JN(rad,v,n)-pSi(rad,v,n)*betaSi*JSi(rad,v,n)<0
                limmat(rad,v,n) =0;
            end
        end
    end
    Jtot(:,:,n)                 = Jmax.*Jlieb(:,:,n)./(Jlieb(:,:,n)+Jmax); % µg C/day
    dgrowth(:,:,n)              = Jtot(:,:,n)./(Ctot);
    
    % Predation mortality
    mu(:,:,n)                   = m0+fG*G(n)*(V0+Vs).^(-1/4);
    
    % Net growth in trait space
    growth(:,:,n)               = dgrowth(:,:,n) - mu(:,:,n) - ms;
    % calculate all optimal parameters
    growth1               = growth(:);
    [g_opt(n),I(n)]       = max(growth1);          % save optimal growth with index
    [I_r(n), I_v(n)]      = ind2sub(size(growth(:,:,1)),I(n)); 
    v_opt(n)              = vfrac(I_r(n),I_v(n));  % optimal vacuole fraction
    r_opt(n)              = ESR(I_r(n),I_v(n));    % optimal size (radius)
    C_opt(:,n)            = c(I_r(n),:);           % cytoplasmic carbon content for optimal vacuole size
    C_opt_p(:,n)          = c(I_r(n),I_v(n));      % optimal cytoplasmic carbon content
    Ctotopt(n)            = Ctot(I_r(n),I_v(n));   % optimal total carbon content(µg C/µm3)
    dcost_opt(n)          = dgrowth(I_r(n),I_v(n),n);  % division rate at optimal growth
    ms_opt(n)             = ms(I_r(n),I_v(n));     % sinking rate at optimal growth
    mu_opt(n)             = mu(I_r(n),I_v(n),n);   % grazing loss at optimal growth
    % upper and lower 95% of optimal parameters
    growth2(:,:,n)              = growth/max(max(growth)); % scaled optimal growth
    [I_r95, I_v95]              = find(round(growth2(:,:,n),2)==0.95);
    if isempty(I_r95) || isempty(I_v95)
        Ir_max(n)               = res;
        Ir_min(n)               = 1;
        Iv_max(n)               = res;
        Iv_min(n)               = 1;
    else
        Ir_max(n)               = max(I_r95);
        Ir_min(n)               = min(I_r95);
        Iv_max(n)               = max(I_v95);
        Iv_min(n)               = min(I_v95);
    end
    Vpercent_opt95u(n)   = vfrac(Ir_max(n),Iv_max(n));     % vacuole percent at upper 95% of optimal growth
    r_opt95u(n)          = ESR(Ir_max(n),Iv_max(n));       % size (radius) at upper 95% of optimal growth
    C_opt95u(n)          = Ctot(Ir_max(n),Iv_max(n));      % carbon content at upper 95% of optimal growth
    Vpercent_opt95l(n)   = vfrac(Ir_min(n),Iv_min(n));     % vacuole percent at lower 95% of optimal growth
    r_opt95l(n)          = ESR(Ir_min(n),Iv_min(n));       % size (radius) at lower 95% of optimal growth
    C_opt95l(n)          = Ctot(Ir_min(n),Iv_min(n));      % carbon content at lower 95% of optimal growth
end
 
%% plots
set(0,'DefaultFigurePosition', [100 100 900 600]);
set(0,'DefaultAxesFontSize','default')
set(0,'DefaultLineLineWidth','default')
set(0,'DefaultFigureColor','w')
monthstr = {'Light limited', 'Non limited', 'Nut/grazing limited'};
color = [0 0 255;30 144 255];
 
fig=figure(1);
fig.Units = 'inches';
p = [7 3 6 5];
set(fig,'pos',p)
lpos2 = 0.13;
lpos22 = 0.06;
lpos1 = 0.08;
lofset = 0.1;
name1 = {'a','b','c','d','e','f'};
for n = 1:length(L)
    h = subplot(3,3,n);
    p = get(h,'pos');
    p(3:4) = p(3:4)*1.2;
    p(2)= p(2)*0.9;
    set(h,'pos',p)
    scatter(log10(C_opt_p(n)*10^6),v_opt(n)*100,g_opt(n)*100, 'filled','d','markerfacecolor','k','markeredgecolor','k'),
    box on
    hold on
    contour(log10(c(1,:)*10^6),vfrac(:,1)*100,growth2(:,:,n),[0.99 0.99],'linecolor',[0 0 255/255],'linewidth', 1.2);
    hold on
    contour(log10(c(1,:)*10^6),vfrac(:,1)*100,growth2(:,:,n), [0.90 0.90],'LineColor',[30 144 255]./255,'linewidth', 1.2);
    hold on
    contour(log10(c(1,:)*10^6),vfrac(:,1)*100,growth2(:,:,n), [0.8 0.8],'LineColor',[130 200 255]./255,'linewidth', 1.2);
    hold on
    contour(log10(c(1,:)*10^6),vfrac(:,1)*100,growth(:,:,n),[0 0],'k','linewidth', 1.2);
    hold on
    plot(log10(c(1,:)*10^6),ones(1,res)*v_opt(n)*100,'--k')% 'color',[0.1 0.1 0,1],'linestyle',)
    hold on
    contour(log10(c(1,:)*10^6),vfrac(:,1)*100,RSiC,[0.052 0.55],'linecolor',[0 1 0],'linewidth',0.7);
    hold on
    scatter(log10(C_opt_p(n)*10^6),v_opt(n)*100,g_opt(n)*100, 'filled','d','markerfacecolor','k','markeredgecolor','k'),
    box on
    if n == 1
        ylabel('Vacuole (%)')
    end
    title(monthstr(2))
    if n==1 || n==2 || n == 3
        set(gca, 'XTickLabel', []);
    end
    if n==2 || n == 3
        set(gca, 'YTickLabel', []);
    end
    if n== 3
        lh = legend('net growth','99%','95%','90%', 'growth = 0');
        set(lh,'orientation','horizontal')
        legend boxoff
        pos = get(lh,'pos');
        pos(1) = 0.12;
        pos(2) = lpos2;
        set(lh,'pos',pos)
    end
    %     h = subplot(3,3,n+3);
    %     p = get(h,'pos');
    %     p(3:4) = p(3:4)*1;
    %     set(h,'pos',p)
    text(-0.8,93,name1(n),'fontweight','bold')
    h = subplot(3,3,n+3);
    p = get(h,'pos');
    p(3:4) = p(3:4)*1.2;
    p(2)= p(2)*0.8;
    set(h,'pos',p)
    plot(log10(c(1,:)*10^6),growth(I_r(n),:,n),'k', 'linewidth',1.5), hold on
    plot(log10(c(1,:)*10^6),-ms(I_r(n),:),'m','linewidth',1.2), hold on
    plot(log10(c(1,:)*10^6),-mu(I_r(n),:,n),'r','linewidth',1.2), hold on
    plot(log10(c(1,:)*10^6),dgrowth(I_r(n),:,n),'b','linewidth',1.2), hold on
    gg = growth(I_r(n),:,n);
    indx = find(gg>0);
    gg = gg(indx);
    tt = log10(c(1,:)*10^6);
    tt = tt(indx);
    gg0 = zeros(size(gg));
    patch([tt(:); flipud(tt(:))],[gg(:); flipud(gg0(:))],[0.8 0.8 0.8], 'EdgeColor', [0.8 0.8 0.8],'Linewidth',1.5), hold on
    axis tight
    box on
    plot(log10(c(1,:)*10^6),growth(I_r(n),:,n),'k', 'linewidth',1.5), hold on
    plot(log10(c(1,:)*10^6),-ms(I_r(n),:),'m','linewidth',1), hold on
    plot(log10(c(1,:)*10^6),-mu(I_r(n),:,n),'r','linewidth',1), hold on
    plot(log10(c(1,:)*10^6),dgrowth(I_r(n),:,n),'b','linewidth',1), hold on
    ylim([-1.5 1.5])
    if n+3==5
        xlabel('Cell carbon (pg)')
    end
    if n+3==4
        ylabel('Rate (day^{-1})')
    end
    if n==2 || n == 3
        set(gca, 'YTickLabel', []);
    end
    if n== 1
        lh = legend('Net growth','Sinking loss','Grazing loss','Division rate');
        set(lh,'orientation','horizontal')
        legend boxoff
        pos = get(lh,'pos');
        pos(2) = lpos22;
        pos(1) = 0.12;
        set(lh,'pos',pos)
    end
    text(-0.8,1.35,name1(n+3),'fontweight','bold')
end
 
fig=figure(2);
fig.Units = 'inches';
p = [7 3 6 5];
name = {'a','b','c'};
set(fig,'pos',p)
for n = 1:length(L)
    h = subplot(3,3,n);
    p = get(h,'pos');
    p(3:4) = p(3:4)*1.2;
    p(2)= p(2)*0.9;
    set(h,'pos',p)
    imagesc(log10(c(1,:)*10^6),vfrac(:,1)*100,limmat(:,:,n),[0 4]), hold on
    scatter(log10(C_opt_p(n)*10^6),v_opt(n)*100,'d','markerfacecolor','k','markeredgecolor','k')
    if n==1
        ylabel('Vacuole (%)')
    end
    if n==2
        xlabel('Cell carbon (pg)')
    end
    % xlim([0 100])
    % ylim([0 rmax])
    if n==2 || n== 3
        set(gca, 'YTickLabel', []);
    end
    title(monthstr(2))
    text(-0.8,5,name(n),'fontweight','bold')
end
 
%% Find optimal growth in trait space
 
% non vacuolated cell
carb    = C_opt;                    % Carbon in cytopasm (µg C)
nitro   = carb*1/RCN_g;             % Nitrogen in cytoplasm (µg C)
size    = C_opt./ccb_g;             % Volume of cytoplasm
radi    = (size./(4/3*pi)).^(1/3);  % radius of cytoplasm (µm)
memvol  = 4/3*pi*(radi+m).^3-size;  % plasmalemma volume (µm3)
totsize = size + memvol;            % Volume including plasmalemma membrane
Cmem1   = rhoCtono*memvol;          % Carbon in plasmalemma µg C
Nmem1   = 1/CNtono*Cmem1*12/14;     % Nitrogen in plasmalemma 
CNtot1  = (carb+Cmem1)./(nitro + Nmem1);% total new CN ratio
rNon    = ((totsize)./(4/3*pi)).^(1/3)*10^-6;% raduis of of cytoplasm and plasmalemma in meters (m)
AL1     = cL.*(totsize).^(2/3);     % µg C/day/(W/m2) (r,v)
AN1     = 4*pi*DN*(rNon)*10^3;      % l/day
Jr01    = 0.05*(carb+Cmem1);        % respiratory cost
Jmax1   = jmax*carb;                % maximum synthesis
 
for v = 1:length(carb(1,:))
    JL1(:,v)    = AL1(:,v)*L(v);
    JN1(:,v)    = AN1(:,v)*N(v);
    for rad     = 1:res
        pNN(rad,v)    = max(0,min((JL1(rad,v) - Jr01(rad,v) - JL1(rad,v)*betaL)/(JN1(rad,v)*betaN + JN1(rad,v)*CNtot1(rad,v)),1));
        Jlieb1(rad,v) = max(min([JL1(rad,v)-Jr01(rad,v)-betaL*JL1(rad,v)-pNN(rad,v)*betaN*JN1(rad,v) pNN(rad,v)*JN1(rad,v)*CNtot1(rad,v)]),0);
    end
    mg(:,v)     = m0 + fG*G(v)*(totsize(:,v)).^(-1/4);
end
 
Jtot1   = Jmax1.*Jlieb1./(Jlieb1+Jmax1);
div     = Jtot1./carb;
w1      = 2*grav*(10^-6*(size./(4/3*pi)).^(1/3)).^2/(9*nu).*(rhoB-rhoW);
ms1     = w1./H*60*60*24;   % sinkin mortality (day^-1)
grow    = div-mg-ms1;
 
% Plot comparison between diaton and non-diatom
fig=figure(3);
fig.Units = 'inches';
p = [7 3 6 5];
set(fig,'pos',p)
lpos2 = 0.13;
lpos22 = 0.06;
lpos1 = 0.15;
lofset = 0.19;
for n=1:length(L)
    h = subplot(3,3,n);
    p = get(h,'pos');
    p(3:4) = p(3:4)*1.2;
    p(2)= p(2)*0.9;
    set(h,'pos',p)
    plot(log10(carb(:,n)*10^6),growth(I_r(n),:,n),'k','LineWidth',1), hold on
    plot(log10(carb(:,n)*10^6),grow(:,n),'--k','LineWidth',1)
    lh = legend('Net growth');
    set(lh,'orientation','horizontal')
    legend boxoff
    pos = get(lh,'pos');
    pos(1) = lpos1;
    pos(2) = lpos2;
    set(lh,'pos',pos)
    if n ~=1
        set(gca, 'YTickLabel', []);
    elseif n == 1
        ylabel('Rate (day^{-1})')
        set(gca, 'XTickLabel', []);
    end
    if n == 2 || n== 3
        set(gca, 'XTickLabel', []);
    end
    ylim([0 1])
    xlim([-1 5])
    title(monthstr(2))
    text(-0.8,0.95,name1(n),'fontweight','bold')
    h = subplot(3,3,n+3);
    p = get(h,'pos');
    p(3:4) = p(3:4)*1.2;
    p(2)= p(2)*0.8;
    set(h,'pos',p)
    plot(log10(carb(:,n)*10^6),dgrowth(I_r(n),:,n),'b','LineWidth',1), hold on
    plot(log10(carb(:,n)*10^6),-ms(I_r(n),:),'m','LineWidth',1), hold on
    plot(log10(carb(:,n)*10^6),-mu(I_r(n),:,n),'r','LineWidth',1), hold on
    plot(log10(carb(:,n)*10^6),-mg(:,n),'--r','LineWidth',1), hold on
    plot(log10(carb(:,n)*10^6),div(:,n),'--b','LineWidth',1), hold on
    plot(log10(carb(:,n)*10^6),-ms1(:,n),'--m','LineWidth',1)
    if n+3==5
        xlabel('Cell carbon (pg)')
    end
    ylim([-1.5 1.5])
    xlim([-1 5])
    if n+3 ~=4
        set(gca, 'YTickLabel', []);
    elseif n+3 == 4
        ylabel('Rate (day^{-1})')
    end
    lh = legend('Division rate', 'Sinking loss', 'Grazing loss');
    set(lh,'orientation','horizontal')
    legend boxoff
    pos = get(lh,'pos');
    pos(1) = lpos1+lofset;
    pos(2) = lpos2;
    set(lh,'pos',pos)
    text(-0.8,1.35,name1(n+3),'fontweight','bold')
end

