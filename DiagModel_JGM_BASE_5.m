% =============================================================
%%                          Carbonate Dissolution-Precipitation Diagenesis Model
%
%           JGM_V5.0: Dissolution of carbonate and precipitation of diagenetic mineral
%
%           2015-08-04      Initial Build                      Anne-Sofie C. Ahm (Princeton Univ)
%           2015-11-30      Revised to include Sr          Anne-Sofie C. Ahm (Princeton Univ)
%           2018-03-09      Published                        Ahm et al. (2018) GCA 236 (2018) 140-159
%           2018-11-07      Revised to include Li           Jack Geary Murphy  (Princeton Univ)
%           2020-07-15      Revised to include HMC        Jack Geary Murphy  (Princeton Univ)
%           2021-02-06      Revised to include prim Dolo  Jack Geary Murphy  (Princeton Univ)
%
% =============================================================

%%  ================== CLEAN SLATE ================================

clear all; clc;                                         % Clears all variables and the command log
disp('Model runtime can take a while (e.g. ~5-20 min) depending on # boxes and timespan')
tic                                                     % start timer

global Fluid alpha R M tune                       % These variables are accessable to all .m files

project = 'Bahamas';                                % Name the project for printing on figures, etc

%%  ================== MODEL VARIABLES =============================
 
%----------------------------- Lithium Knobs  ----------------------------------------

tune.A =   0.63;                % Prim Arag [Li]   (ppm)
tune.B =   0.63;                % Prim LMC [Li]    
tune.C =  4.4;                  % Prim HMC [Li]   
tune.D =  1;                    % Prim Dolo [Li]    

tune.G =  20;                   % Prim Arag d7Li  (permil)
tune.E = 29;                     % Prim LMC d7Li   
tune.F = 24;                     % Prim HMC d7Li   
tune.M = 31;                    % Prim Dolo d7Li

% tune.G =  14;                     % Prim Arag d7Li  (permil)
% tune.E = 23;                       % Prim LMC d7Li   
% tune.F = 18;                       % Prim HMC d7Li   
% tune.M = 25.5;                    % Prim Dolo d7Li


tune.H = 2.56e-5;             % Prim fluid [Li]    (mol/kg)
tune.I =  31.5;                 % Prim fluid d7Li   (permil)
% tune.I =  25.5;                 % Prim fluid d7Li   (permil)

tune.L = 0.01;                 % Distribution coeffient (Kd = solid/fluid) for Sr/Ca 
tune.J = 0.008;                % Kd for Li into diagenetic carbonate (Li/Ca)   
                                                  % Note: Marriott (2004) found dependence on T and salinity
                                                  %          inorg. calcite: D =      0.01e(-0.045*T)
                                                  %                                     0.010 @ 5*C
                                                  %                                     0.003 @ 30*C
                                                  %                                     0.0009 @ 10 psu, 25*C
                                                  %                                     0.0030 @ 50 psu, 25*C
                                                  %  Estimated temperature of alteration 12-16*C (Higgins 2018)
tune.K = 1.000;                % alpha Li into diagenetic carbonate      



%----------------------- Composition of Primary Rock -----------------------------------

% Initial mineralogy:
Dolomite = 1e-6;                            % Fraction of initial rock that is primary dolomite (should be > 0)!
LMC     = 1e-6;                             %      "          "           "         LMC (has to be > 0)
HMC = 1e-6;                                  %       "          "           "        HMC (has to be > 0)
Aragonite = 1-(Dolomite+LMC+HMC);   % Fraction of aragonite to make sure total sums to 100%
                                               % Note: no initial diagenetic mineral


% INITIAL ISOTOPES
% Aragonite:                        LMC:                               HMC:
Solid.arag.dCa   = -1.6;          Solid.LMC.dCa   = -1.1;         Solid.HMC.dCa   = -1.1;            
Solid.arag.dMg  = -3.7;           Solid.LMC.dMg  = -3;          Solid.HMC.dMg  = -3;                 
Solid.arag.dC     =  5;           Solid.LMC.dC    =  0;            Solid.HMC.dC =  0;                   
Solid.arag.dO    =  0;            Solid.LMC.dO    = 0;             Solid.HMC.dO    = 2;                  
Solid.arag.dSr   =   0.70900;    Solid.LMC.dSr   =   0.70900;    Solid.HMC.dSr   =   0.70900;      
Solid.arag.dLi   =   tune.G;       Solid.LMC.dLi   =   tune.E;       Solid.HMC.dLi   =   tune.F;         

% Dolomite:                        % Diagenetic Mineral: (need to be defined -- leave as 0)
Solid.dolo.dCa = -1.1;              Solid.DIAG.dCa = 0;      % d44Ca value                
Solid.dolo.dMg = -2.8;              Solid.DIAG.dMg = 0;     % d26Mg value
Solid.dolo.dC = 4;                   Solid.DIAG.dC = 0;       % d13C value
Solid.dolo.dO = 1;                   Solid.DIAG.dO = 0;       % d18O value
Solid.dolo.dSr =  0.70900;         Solid.DIAG.dSr = 0;      % 87Sr/86Sr ratio
Solid.dolo.dLi = tune.M;             Solid.DIAG.dLi = 0;      % d7Li value

%---------------------- Composition of Diagenetic Fluid ---------------------------------

% Concentration
Ca_fluid  = 0.0106;           % Ca in fluid (mol/kg)
Mg_fluid  = 0.0528;            % Mg in fluid (mol/kg)
C_fluid    = 0.0024;           % C in fluid (mol/kg)
Sr_fluid   = 0.00009;          % Sr in fluid (mol/kg)
% Sr_fluid   = 0.000600;          % Sr in fluid (mol/kg)
Li_fluid   = tune.H;             % Li in fluid (mol/kg)

% Isotopic composition 
Fluid.dCa   = 0;                % d44Ca_fluid value
Fluid.dMg  = -0.8;             % d26Mg_fluid value  
Fluid.dC    = -2;               % d13C_fluid value   
Fluid.dO    = -30.5;           % d18O_fluid value (normalized to V-PDB)
Fluid.dSr   = 0.709250;       % 86Sr/87Sr_fluid value
Fluid.dLi    = tune.I;           % d7Li_fluid value


% Fractionation factors of diagenetic mineral 
% alpha.Ca  = 1.000;            % d44Ca_alpha 
alpha.Ca  = 0.9998;            % d44Ca_alpha 
alpha.Mg = 0.998;             % d26Mg_alpha 
alpha.C   = 1.001;             % d13C_alpha
alpha.O   = 1.0345;           % d18O_alpha 
alpha.Sr   = 1;                 % 87Sr/86Sr_alpha (should be 1)   
alpha.Li = tune.K;              % d7Li_alpha 


%------------------------- Diagenetic Length Scale -------------------------------------

u = 0.1;                     % Advection Rate (m/yr)   
% u = 0.005;                 % Advection Rate (m/yr)   

R.arag = 1e-5;             % Reaction Rates (1/yr) 
% R.arag = 1e-2;             % Reaction Rates (1/yr) 
R.HMC = R.arag;
R.LMC = R.arag;
R.dolo = R.arag;
R.DIAG = 0;                 % if you want to re-dissolve diagenetic mineral

    % NOTE! 
    % Depending on the reaction rate, the problem can become too stiff for the ODE solver
    % if there are too many boxes relative to the advection/reaction rate. 
    % You can try a different # boxes / rates, or try a different solver --
    % ode15s is probably he best solver for stiff problems - change @ line 200
    
time = 1e7;                 % Timespan modeled (yrs)
% time = 1e6;                 % Timespan modeled (yrs)
box = 40;                   % Number of boxes (= sediment depth) (m)
% box = 10;                   % Number of boxes (= sediment depth) (m)


%--------------------- Stoichiometry of diagenetic mineral --------------------------------

M =  0.5;                 % Mg/Ca (mol/mol) ratio of diagenetic mineral 
                               %          dolomite        = 0.5 
                               %          low-Mg calcite = 0.035 
                               %          meteoric LMC  = 0.0035

                             
%% =================== CONSTANTS ================================

% Function that calls constants:

[Molar, Box, Solid] = DiagModel_JGM_constants_5(Aragonite, LMC, HMC, Dolomite, Ca_fluid, Mg_fluid, C_fluid, Sr_fluid, Li_fluid, Solid, alpha);

        % Molar = elemental molar mass
        % Box    = box volume and porosity
        % Solid   = Initial composition of bulk rock
        %            (e.g. concentration of elements in primary minerals)
                             
                  
%%  ================== MODEL INPUT VECTOR  ==========================


% Initial conditions vector (to be passed to ODE):

J = [ Solid.LMC.Ca;   Solid.LMC.Mg;    Solid.LMC.C;   Solid.LMC.O;    Solid.LMC.Sr;    Solid.LMC.Li;...        % moles X in LMC
      Solid.HMC.Ca;   Solid.HMC.Mg;    Solid.HMC.C;   Solid.HMC.O;    Solid.HMC.Sr;    Solid.HMC.Li;...     % moles X in HMC
      Solid.arag.Ca;   Solid.arag.Mg;   Solid.arag.C;   Solid.arag.O;   Solid.arag.Sr;   Solid.arag.Li;...         % moles X in aragonite
      Solid.dolo.Ca;   Solid.dolo.Mg;   Solid.dolo.C;   Solid.dolo.O;   Solid.dolo.Sr;   Solid.dolo.Li;...         % moles X in dolomite
      Solid.DIAG.Ca;   Solid.DIAG.Mg;   Solid.DIAG.C;   Solid.DIAG.O;   Solid.DIAG.Sr;   Solid.DIAG.Li;...      % moles X in diagenetic mineral
      Solid.LMC.dCa;  Solid.LMC.dMg;  Solid.LMC.dC;  Solid.LMC.dO; Solid.LMC.dSr;  Solid.LMC.dLi;...        % delta X of LMC
      Solid.HMC.dCa;  Solid.HMC.dMg;  Solid.HMC.dC;  Solid.HMC.dO; Solid.HMC.dSr;  Solid.HMC.dLi;...     % delta X of HMC
      Solid.arag.dCa; Solid.arag.dMg;  Solid.arag.dC; Solid.arag.dO; Solid.arag.dSr; Solid.arag.dLi;...         % delta X of aragonite
      Solid.dolo.dCa; Solid.dolo.dMg;  Solid.dolo.dC; Solid.dolo.dO; Solid.dolo.dSr; Solid.dolo.dLi;...         % delta X of dolomite
      Solid.DIAG.dCa; Solid.DIAG.dMg;  Solid.DIAG.dC; Solid.DIAG.dO; Solid.DIAG.dSr; Solid.DIAG.dLi;...      % delta X of diagenetic mineral
      Solid.dCa;        Solid.dMg;         Solid.dC;         Solid.dO;        Solid.dSr;        Solid.dLi;...        % molesX*dX in all solid 
      Fluid.Ca;         Fluid.Mg;            Fluid.C;          Fluid.O;          Fluid.Sr;          Fluid.Li;...        % moles X in fluid
      Fluid.dCa_init; Fluid.dMg_init;   Fluid.dC_init;   Fluid.dO_init;  Fluid.dSr_init;  Fluid.dLi_init;...          % mol*dX in fluid
      u];                                                                                                                  % Advection Rate (m/yr)


%%  =================== Initialize the Model ===========================

% MODEL DIMENSIONS
    tspan = [0 time];           % time span (yr)
    m = 1:box;                   % vector w/ number of boxes
    v = length(J);                % Number of variables tracked by the model

    
% INITIALIZE ALL BOXES
%   i.e. apply the same initial conditions from J to all boxes, 
%   size( Jt0 ) = [ var , boxes ] 

    for i = 1:length(m) 
        Jt0(:,i) = J(:,1); 
    end  

% ODE SOLVER OPTIONS    
    options = odeset('OutputFcn',@odewbar);     % monitor model progress

%% =================== Evolve the Model =============================
%                                 solve (dJ/dt) using ODE solver

%   Inputs
%       ode15s         identify the ODE solver to use (this one's good for "stiff" problems)
%       @(t,j)odefun   filename of dJdt function that stores the eqns to be solved
%       (t,j,m,v)         input arguments used to ODE functions script
%       tspan           timespan range to integrate over -- solver will determine most efficient timesteps to take                      
%       Jt0               Initial conditions of every cell, Bm.Vv(t0), stored as [v m] 2D array
%       options         option to apply to ODE solver     
% 
%       Note: ODE solver takes in 2d array Jt0 and reshapes it into a row vector, 
%              then evolves each cell by adding rows with each timesteps.
%
%      [ B1.V1(t0)   ...   Bm.V1(t0) ]
%      [ B1.V2(t0)   ...   Bm.V2(t0) ]  ==>  [ B1.V1 (t0)  B1.V2 (t0) ....  B2.V1 (t0)  B2.V2 (t0) .... Bm.Vv (t0) ]
%      [      :          ...        :          ]
%
%   Outputs
%       t:  timesteps taken by the model [ will vary depending on slope of dJ(i)/dt(i) ]
%       J:  [ B1.V1 (t0)   B1.V2 (t0) ....  B2.V1 (t0)   B2.V2 (t0) .... Bm.Vv (t0)  ]
%           [ B1.V1 (t1)   B1.V2 (t1) ....  B2.V1 (t1)   B1.V2 (t1) .... Bm.Vv (t1)  ]
%           [        :               :                    :                  :                  :           ]
%           [ B1.V1 (tf)    B1.V2 (tf) ....  B2.V1 (tf)    B2.V2 (tf)  ....  Bm.Vv (tf) ]

     [t,J] = ode15s(@(t,J)DiagModel_JGM_dJdt_5(t,J,m,v),tspan,Jt0, options);   
     
     
     
%% ====================== Model Output ============================
    
% STORE MODEL RUN

    for i = 1:length(m) 
            j = v*i - v;                                                          
            a = j + 1;                                 % index beginning of variables for new box i                                   
            c = j + v;                                 % index end of variables for box i
            JPrime(:,:,i) = J(1:length(t),a:c)';  % transpose J (time vs tracked variables) for box  i and
                                                           % store variables(row) vs time (col) into box i (file)
    end

% Note: JPrime(var, t, box) stores all output data from model run. 
%          rows= tracked variables     
%          columns = timestep     
%          sheet = box
%
%  save model run:    save(['DiagModelRun_', datestr(now, 'yyyy-mm-dd (HH|MM)')]) 


% Extract variables from JPrime (3d matrix)

for i = 1:length(t)             % timestep
    for j = 1:length(m)       % box

% Elemental composition of LMC (mol):
Solid.LMC.Ca(i,j) = JPrime(1,i,j);          
Solid.LMC.Mg(i,j) = JPrime(2,i,j);          
Solid.LMC.C(i,j) = JPrime(3,i,j);            
Solid.LMC.O(i,j) = JPrime(4,i,j);           
Solid.LMC.Sr(i,j) = JPrime(5,i,j);
Solid.LMC.Li(i,j) = JPrime(6,i,j);

% Elemental composition of HMC (mol):
Solid.HMC.Ca(i,j) = JPrime(7,i,j);          
Solid.HMC.Mg(i,j) = JPrime(8,i,j);          
Solid.HMC.C(i,j) = JPrime(9,i,j);            
Solid.HMC.O(i,j) = JPrime(10,i,j);           
Solid.HMC.Sr(i,j) = JPrime(11,i,j);
Solid.HMC.Li(i,j) = JPrime(12,i,j);

% Elemental composition of aragonite (mol):
Solid.arag.Ca(i,j) = JPrime(13,i,j);
Solid.arag.Mg(i,j) = JPrime(14,i,j);
Solid.arag.C(i,j) = JPrime(15,i,j);
Solid.arag.O(i,j) = JPrime(16,i,j);
Solid.arag.Sr(i,j) = JPrime(17,i,j);
Solid.arag.Li(i,j) = JPrime(18,i,j);

% Elemental composition of the inital dolomite (mol):
Solid.dolo.Ca(i,j) = JPrime(19,i,j);
Solid.dolo.Mg(i,j) = JPrime(20,i,j);
Solid.dolo.C(i,j) = JPrime(21,i,j);
Solid.dolo.O(i,j) = JPrime(22,i,j);
Solid.dolo.Sr(i,j) = JPrime(23,i,j);
Solid.dolo.Li(i,j) = JPrime(24,i,j);

% Elemental composition of the diagenetic mineral (mol):
Solid.DIAG.Ca(i,j) = JPrime(25,i,j);
Solid.DIAG.Mg(i,j) = JPrime(26,i,j);
Solid.DIAG.C(i,j) = JPrime(27,i,j);
Solid.DIAG.O(i,j) = JPrime(28,i,j);
Solid.DIAG.Sr(i,j) = JPrime(29,i,j);
Solid.DIAG.Li(i,j) = JPrime(30,i,j);

% Total elemental composition of the bulk sediment (mol): 
Total_Ca(i,j) = (Solid.LMC.Ca(i,j) + Solid.HMC.Ca(i,j) + Solid.arag.Ca(i,j) + Solid.dolo.Ca(i,j) + Solid.DIAG.Ca(i,j));
Total_Mg(i,j) = (Solid.LMC.Mg(i,j) + Solid.HMC.Mg(i,j) + Solid.arag.Mg(i,j) + Solid.dolo.Mg(i,j) + Solid.DIAG.Mg(i,j));
Total_C(i,j) = (Solid.LMC.C(i,j) + Solid.HMC.C(i,j) + Solid.arag.C(i,j) + Solid.dolo.C(i,j) + Solid.DIAG.C(i,j));
Total_O(i,j) = (Solid.LMC.O(i,j)+ Solid.HMC.O(i,j) + Solid.arag.O(i,j) + Solid.dolo.O(i,j) + Solid.DIAG.O(i,j));
Total_Sr(i,j) = (Solid.LMC.Sr(i,j) + Solid.HMC.Sr(i,j) + Solid.arag.Sr(i,j) + Solid.dolo.Sr(i,j) + Solid.DIAG.Sr(i,j));
Total_Li(i,j) = (Solid.LMC.Li(i,j) + Solid.HMC.Li(i,j) + Solid.arag.Li(i,j) + Solid.dolo.Li(i,j) + Solid.DIAG.Li(i,j));

% Calculate elemetal ratios in the bulk sediment :
Solid.Sr_Ca = Total_Sr./Total_Ca*1e6;    % (umol/mol)
Solid.Mg_Ca = Total_Mg./Total_Ca*1e3;  % (mmol/mol)
Solid.Li_Ca = Total_Li./Total_Ca*1e6;     % (umol/mol)

Solid.Sr_CaMg = Total_Sr./(Total_Ca+Total_Mg)*1e6;     % (umol/mol)
Solid.Mg_CaMg = Total_Mg./(Total_Ca+Total_Mg)*1e3;   % (mmol/mol)
Solid.Li_CaMg = Total_Li./(Total_Ca+Total_Mg)*1e6;      % (umol/mol)


% Total moles of the bulk sediment:
FWs(i,j) = Total_Ca(i,j) + Total_Mg(i,j) + Total_C(i,j) + Total_O(i,j) + Total_Sr(i,j) + Total_Li(i,j);

% Isotopic composition of the diagenetic mineral (permil):
Solid.DIAG.dCa(i,j) = JPrime(55,i,j)/Solid.DIAG.Ca(i,j);
Solid.DIAG.dMg(i,j) = JPrime(56,i,j)/Solid.DIAG.Mg(i,j);
Solid.DIAG.dC(i,j) = JPrime(57,i,j)/Solid.DIAG.C(i,j);
Solid.DIAG.dO(i,j) = JPrime(58,i,j)/Solid.DIAG.O(i,j);
Solid.DIAG.dSr(i,j) = JPrime(59,i,j)/Solid.DIAG.Sr(i,j);
Solid.DIAG.dLi(i,j) = JPrime(60,i,j)/Solid.DIAG.Li(i,j);

% Isotopic composition of the bulk sediment (moles*dX)/(mol) = permil:
Solid.dCa(i,j) = JPrime(61,i,j)/Total_Ca(i,j);               
Solid.dMg(i,j) = JPrime(62,i,j)/Total_Mg(i,j);
Solid.dC(i,j) = JPrime(63,i,j)/Total_C(i,j);
Solid.dO(i,j) = JPrime(64,i,j)/Total_O(i,j);
Solid.dSr(i,j) = JPrime(65,i,j)/Total_Sr(i,j);
Solid.dLi(i,j) = JPrime(66,i,j)/Total_Li(i,j);

% Elemental composition of the diagenetic fluid (mol):
Fluid.Ca(i,j) = JPrime(67,i,j);
Fluid.Mg(i,j) = JPrime(68,i,j);
Fluid.C(i,j) = JPrime(69,i,j);
Fluid.O(i,j) = JPrime(70,i,j);
Fluid.Sr(i,j) = JPrime(71,i,j);
Fluid.Li(i,j) = JPrime(72,i,j);

% Total moles of the diagenetic fluid:
FWf(i,j) = Fluid.Ca(i,j) + Fluid.Mg(i,j) + Fluid.C(i,j) + Fluid.O(i,j) + Fluid.Sr(i,j) + Fluid.Li(i,j);

% Isotopic composition of the diagenetic fluid (moles*dX)/(mol) = permil:
Fluid.dCa(i,j) = JPrime(73,i,j)/Fluid.Ca(i,j);
Fluid.dMg(i,j) = JPrime(74,i,j)/Fluid.Mg(i,j);
Fluid.dC(i,j) = JPrime(75,i,j)/Fluid.C(i,j);
Fluid.dO(i,j) = JPrime(76,i,j)/Fluid.O(i,j);
Fluid.dSr(i,j) = JPrime(77,i,j)/Fluid.Sr(i,j);
Fluid.dLi(i,j) = JPrime(78,i,j)/Fluid.Li(i,j);

% Cummulative fluid-to-rock ratio (unitless):
% [yrs*(box/yr)*(kg of fluid in box/box)+(kg fluid init. in box) ] / (kg solid in box)

WR(i) = (t(i)*u*Box.fluid + Box.fluid)/Box.sed;         

% Percent alteration:
Solid.procent(i,j) = Solid.DIAG.C(i,j)./Total_C(i,j);

    end 
end


%%  ====================== Figures ================================

% ------------------------ FIG 1: DIAGENETIC EVOLUTION --------------------------------
% ------------ Change in bulk rock chemistry with increasing fluid-to-rock ratios ------------ 

figure

suptitle('Progressive Diagenetic Alteration')

% Water/Rock Ratio
subplot(7,1,1)

    semilogx(WR,Solid.procent(:,end)); hold on
    semilogx(WR,1-Solid.procent(:,end));
        xlim([1e1 WR(end)]); ylabel('%')

% Carbon Isotopes       
subplot(7,1,2)

    semilogx(WR,Solid.dC(:,[1,end])); hold on
        xlim([1e1 WR(end)]); ylabel('\delta^{13}C')

% Oxygen Isotopes          
subplot(7,1,3)

    semilogx(WR,Solid.dO(:,[1,end])); hold on
        xlim([1e1 WR(end)]); ylabel('\delta^{18}O')

% Calcium Isotopes 
subplot(7,1,4)

    semilogx(WR,Solid.dCa(:,[1,end])); hold on
        xlim([1e1 WR(end)]); ylabel('\delta^{44}Ca')

% Magnesium Isotopes        
subplot(7,1,5)

    semilogx(WR,Solid.dMg(:,[1,end])); hold on
        xlim([1e1 WR(end)]); ylabel('\delta^{26}Mg')

% Sr/(Ca+Mg)           
subplot(7,1,6)

    semilogx(WR,Solid.Sr_CaMg(:,[1,end])); hold on
        xlim([1e1 WR(end)]); ylabel('Sr/(Ca+Mg) (\mumol/mol)')

% Radiogenic Strontium Isotopes        
subplot(7,1,7)

semilogx(WR,Solid.dSr(:,[1,end])); hold on
        xlim([1e1 WR(end)]); ylabel('^{87}Sr/^{86}Sr')
        xlabel('Cumulative fluid-to-rock ratio')
        
legend('box 1','box n')
        
% % ------------------------ FIG 2: LITHIUM DIAGENETIC EVOLUTION -----------------------

figure     
suptitle('Progressive Diagenetic Alteration')
 
% Li/(Ca+Mg)        
subplot(2,1,1)

    semilogx(WR,Solid.Li_CaMg(:,[1,end])); hold on
        xlim([1e1 WR(end)]); ylabel('Li/(Ca+Mg) (\mumol/mol)')
        
% Lithium Isotopes        
subplot(2,1,2)

semilogx(WR,Solid.dLi(:,[1,end])); hold on
        xlim([1e1 WR(end)]); ylabel('\delta^{7}Li')
        xlabel('Cumulative fluid-to-rock ratio')

legend('box 1','box n')


% % ---------------------------- FIG 3: CROSSPLOTS -------------------------------------
% % ------------ Model phase space defined by covariation between pairs of proxies ----------
   
figure

suptitle('Model Phase Space')

colr = [0.7 0.7 0.7]; 

% Sr/(Ca+Mg) vs d44Ca
subplot(3,2,1)

    plot(Solid.dCa(:,[1,end]),Solid.Sr_CaMg(:,[1,end]),'Color','k'); hold on
    contour(Solid.dCa,Solid.Sr_CaMg,Solid.procent,'-','Color',[0.8 0.8 0.8]);
    plot(Solid.dCa(end,:),Solid.Sr_CaMg(end,:),'Color',[0.1 0.7 1],'LineWidth',1);
        xlabel(['\delta^{44/40}Ca (',char(8240),')']); 
        ylabel('Sr/(Ca+Mg) (\mumol/mol)'); 

% d18O vs d44Ca       
subplot(3,2,2)

    plot(Solid.dCa(:,[1,end]),Solid.dO(:,[1,end]),'Color','k'); hold on
    contour(Solid.dCa,Solid.dO,Solid.procent,'-','Color',[0.8 0.8 0.8]);
    plot(Solid.dCa(end,:),Solid.dO(end,:),'Color',[0.1 0.7 1],'LineWidth',1);

        xlabel(['\delta^{44/40}Ca (',char(8240),')']); 
        ylabel(['\delta^1^8O(',char(8240),')']); 

% d13C vs d44Ca        
subplot(3,2,3)

    plot(Solid.dCa(:,[1,end]),Solid.dC(:,[1,end]),'Color','k'); hold on
    contour(Solid.dCa,Solid.dC,Solid.procent,'-','Color',[0.8 0.8 0.8]);
    plot(Solid.dCa(end,:),Solid.dC(end,:),'Color',[0.1 0.7 1],'LineWidth',1);

        xlabel(['\delta^{44/40}Ca (',char(8240),')']);
        ylabel('\delta^1^3C');

% d24Mg vs d44Ca        
subplot(3,2,4)

    plot(Solid.dCa(:,[1,end]),Solid.dMg(:,[1,end]),'Color','k'); hold on
    contour(Solid.dCa,Solid.dMg,Solid.procent,'-','Color',[0.8 0.8 0.8]);
    plot(Solid.dCa(end,:),Solid.dMg(end,:),'Color',[0.1 0.7 1],'LineWidth',1);

        xlabel(['\delta^{44/40}Ca (',char(8240),')']); 
        ylabel(['\delta^2^6Mg(',char(8240),')']); 

% Li/(Ca+Mg) vs d44Ca
subplot(3,2,5)

    plot(Solid.dCa(:,[1,end]),Solid.Li_CaMg(:,[1,end]),'Color','k'); hold on
    contour(Solid.dCa,Solid.Li_CaMg,Solid.procent,'-','Color',[0.8 0.8 0.8]);
    plot(Solid.dCa(end,:),Solid.Li_CaMg(end,:),'Color',[0.1 0.7 1],'LineWidth',1);

        xlabel(['\delta^{44/40}Ca (',char(8240),')']); 
        ylabel('Li/(Ca+Mg) (\mumol/mol)');    
        
% d7Li vs d44Ca        
subplot(3,2,6)

    plot(Solid.dCa(:,[1,end]),Solid.dLi(:,[1,end]),'Color','k'); hold on
    contour(Solid.dCa,Solid.dLi,Solid.procent,'-','Color',[0.8 0.8 0.8]);
    plot(Solid.dCa(end,:),Solid.dLi(end,:),'Color',[0.1 0.7 1],'LineWidth',1);

        xlabel(['\delta^{44/40}Ca (',char(8240),')']); 
        ylabel(['\delta^7Li(',char(8240),')']); 
              
 figure
 suptitle('Lithium Model Phase Space')
 
 % d7Li vs d44Ca       
subplot(2,1,1)
    plot(Solid.dCa(:,[1,end]),Solid.dLi(:,[1,end]),'Color','k'); hold on
    contour(Solid.dCa,Solid.dLi,Solid.procent,'-','Color',[0.8 0.8 0.8]);
    plot(Solid.dCa(end,:),Solid.dLi(end,:),'Color',[0.1 0.7 1],'LineWidth',1);

        xlabel(['\delta^{44/40}Ca (',char(8240),')']); 
        ylabel(['\delta^7Li(',char(8240),')']); 
        xlim([-2 0]); ylim([15 35])
        
  % d7Li vs Sr/(Ca+Mg)       
subplot(2,1,2)
    plot(Solid.Sr_CaMg(:,[1,end]),Solid.dLi(:,[1,end]),'Color','k'); hold on
    contour(Solid.Sr_CaMg,Solid.dLi,Solid.procent,'-','Color',[0.8 0.8 0.8]);
    plot(Solid.Sr_CaMg(end,:),Solid.dLi(end,:),'Color',[0.1 0.7 1],'LineWidth',1);

        xlabel('Sr/(Ca+Mg) (\mumol/mol)'); 
        ylabel(['\delta^7Li(',char(8240),')']); 
        xlim([0 15000]); ylim([15 35])


   
% savediagmodelparaset_v5
% savediagmodelrun
     
toc                 % Stop timer
disp('              ~~   FIN   ~~ ')
    
% =============== ALL DONE! ============================
