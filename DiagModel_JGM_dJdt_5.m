% =============================================================
%%                                                    MODEL EVOLUTION
%                                                          by ODE SOLVER 
%
%                   Function that reads the value of each tracked variable (v[i]) for each box (n) at  
%                   timestep (t), then calculates the finite difference dv(i)/dt for timestep dt from 
%                   the mass balance of the fluxes [ dv(i)/dt = inputs - outputs ].
%
%                   Calculates next timestep for all variables for each box, then repeats
%                    
% =============================================================

function dy = DiagModel_JGM_dJdt_5(t,J,m,v)

global Fluid

%% ==================== Initialization ===============================

box = length(m);       % number of boxes
dy = zeros(box*v,1);   % initialize dy vector which will store all variables (v) 
                           % for all boxes (m) for a specific timestep (t)
                                 
fluid.Ca = zeros(box+1,1);       fluid.dCa = zeros(box+1,1);    % initialize fluid chemistry vectors
fluid.Mg = zeros(box+1,1);      fluid.dMg = zeros(box+1,1);
fluid.C = zeros(box+1,1);        fluid.dC = zeros(box+1,1);
fluid.O = zeros(box+1,1);        fluid.dO = zeros(box+1,1);
fluid.Sr = zeros(box+1,1);       fluid.dSr = zeros(box+1,1);
fluid.Li = zeros(box+1,1);       fluid.dLi = zeros(box+1,1);

%% ================ Box 1 Boundary Conditions ==========================

% Input to Box 1 is always sea-water / meteoric water

fluid(1).Ca = Fluid.Ca;             % (mol)
fluid(1).Mg = Fluid.Mg;
fluid(1).C = Fluid.C;
fluid(1).O = Fluid.O;
fluid(1).Sr = Fluid.Sr;
fluid(1).Li = Fluid.Li;

fluid(1).dCa = Fluid.dCa;         % (dX) 
fluid(1).dMg = Fluid.dMg;
fluid(1).dC = Fluid.dC;
fluid(1).dO = Fluid.dO;
fluid(1).dSr = Fluid.dSr;
fluid(1).dLi = Fluid.dLi;

%% ================ Read Present State of Each Box  =======================

for n = 1:box
   
    i = v*n - v;                        %  index into a vector with ODE results for all  
                                          %  variables in boxes 1:n at a given timestep.
                                          %  (note: we input a vector rather than a matrix so that the ...    
                                          %          ODE solver is iterating a single vector through time.)
    
% Update Elemental Composition of Solid
% (update values for box n at timestep t)

LMC.Ca = J(i+1);                            % Ca in LMC (mol)
LMC.Mg = J(i+2);                           % Mg in LMC (mol)
LMC.C = J(i+3);                             % C in LMC (mol)
LMC.O = J(i+4);                             % O in LMC (mol)
LMC.Sr = J(i+5);                            % Sr in LMC (mol)
LMC.Li = J(i+6);                            % Li in LMC (mol)

HMC.Ca = J(i+7);                            % Ca in HMC (mol)
HMC.Mg = J(i+8);                           % Mg in HMC (mol)
HMC.C = J(i+9);                             % C in HMC (mol)
HMC.O = J(i+10);                             % O in HMC (mol)
HMC.Sr = J(i+11);                            % Sr in HMC (mol)
HMC.Li = J(i+12);                            % Li in HMC (mol)

arag.Ca = J(i+13);                           % Ca in aragonite (mol)                  
arag.Mg = J(i+14);                          % Mg in aragonite (mol)
arag.C = J(i+15);                            % C in aragonite (mol)
arag.O = J(i+16);                            % O in aragonite (mol)
arag.Sr = J(i+17);                          % Sr in aragonite (mol)
arag.Li = J(i+18);                          % Li in aragonite (mol)

dolo.Ca = J(i+19);                          % Ca in dolomite (mol)
dolo.Mg = J(i+20);                         % Mg in dolomite (mol)
dolo.C = J(i+21);                           % C in dolomite (mol)
dolo.O = J(i+22);                           % O in dolomite (mol)
dolo.Sr = J(i+23);                          % Sr in dolomite (mol)
dolo.Li = J(i+24);                          % Li in dolomite (mol)

DIAG.Ca = J(i+25);                          % Ca in diagenetic mineral (mol)
DIAG.Mg = J(i+26);                         % Mg in diagenetic mineral (mol)
DIAG.C = J(i+27);                           % C in diagenetic mineral (mol)
DIAG.O = J(i+28);                           % O in diagenetic mineral (mol)
DIAG.Sr = J(i+29);                          % Sr in diagenetic mineral (mol)
DIAG.Li = J(i+30);                          % Li in diagenetic mineral (mol)


% Update Isotope Composition of Solid
% (update values for box n at timestep t)

LMC.dCa = J(i+31);                         % d44Ca in LMC (per mil) 
LMC.dMg = J(i+32);                        % d26Mg in LMC (per mil)
LMC.dC = J(i+33);                          % d13C in LMC (per mil)
LMC.dO = J(i+34);                          % d18O in LMC (per mil)
LMC.dSr = J(i+35);                         % 86Sr/87Sr in LMC
LMC.dLi = J(i+36);                         % d7Li  in LMC (per mil)

HMC.dCa = J(i+37);                        % d44Ca in HMC (per mil)
HMC.dMg = J(i+38);                       % d26Mg in HMC (per mil)
HMC.dC = J(i+39);                         % d13C in HMC (per mil)
HMC.dO = J(i+40);                         % d18O in HMC (per mil)
HMC.dSr = J(i+41);                         % 86Sr/87Sr in HMC 
HMC.dLi = J(i+42);                         % d7Li in HMC (per mil)

arag.dCa = J(i+43);                        % d44Ca in aragonite (per mil)
arag.dMg = J(i+44);                       % d26Mg in aragonite (per mil)
arag.dC = J(i+45);                         % d13C in aragonite (per mil)
arag.dO = J(i+46);                         % d18O in aragonite (per mil)
arag.dSr = J(i+47);                         % 86Sr/87Sr in aragonite
arag.dLi = J(i+48);                         % d7Li in aragonite (per mil)

dolo.dCa = J(i+49);                         % d44Ca in dolomite (per mil)
dolo.dMg = J(i+50);                         % d26Mg in dolomite (per mil)
dolo.dC = J(i+51);                           % d13C in dolomite (per mil)
dolo.dO = J(i+52);                           % d18O in dolomite (per mil)
dolo.dSr = J(i+53);                          % 86Sr/87Sr in dolomite
dolo.dLi = J(i+54);                          % d7Li in dolomite (per mil)

DIAG.dCa = J(i+55);                         % d44Ca in diagenetic mineral (per mil)
DIAG.dMg = J(i+56);                         % d26Mg in diagenetic mineral (per mil)
DIAG.dC = J(i+57);                           % d13C in diagenetic mineral (per mil)
DIAG.dO = J(i+58);                           % d18O in diagenetic mineral (per mil)
DIAG.dSr = J(i+59);                          % 86Sr/87Sr in diagenetic mineral
DIAG.dLi = J(i+60);                          % d7Li in diagenetic mineral (per mil)

% Mass of the fluid
% (store values for box n at timesetp t)
fluid(n+1).Ca = J(i+67);                    % Ca in fluid (mol)
fluid(n+1).Mg = J(i+68);                   % Mg in fluid (mol)
fluid(n+1).C = J(i+69);                     % C in fluid (mol)
fluid(n+1).O = J(i+70);                     % O in fluid (mol)
fluid(n+1).Sr = J(i+71);                    % Sr in fluid (mol)
fluid(n+1).Li = J(i+72);                    % Li in fluid (mol)

% Isotope composition of fluid:
% (store values for box n at timesetp t)
fluid(n+1).dCa = J(i+73)/fluid(n+1).Ca;      % d44Ca in fluid (mol*permil/mol)
fluid(n+1).dMg = J(i+74)/fluid(n+1).Mg;    % d26Mg in fluid (mol*permil/mol)
fluid(n+1).dC = J(i+75)/fluid(n+1).C;         % d13C in fluid (mol*permil/mol)
fluid(n+1).dO = J(i+76)/fluid(n+1).O;        % d18O in fluid (mol*permil/mol)
fluid(n+1).dSr = J(i+77)/fluid(n+1).Sr;       % 86Sr/87Sr in fluid (mol*permil/mol)
fluid(n+1).dLi = J(i+78)/fluid(n+1).Li;        % d7Li in fluid (mol*permil/mol)

% Flux
flux = J(i+79);      % Fluid flow rate (m/yr)


%% ==================== Calculate Fluxes =============================
%                          Call function that solves all input and output fluxes to each box

 [in, out] = DiagModel_JGM_fluxes_5(LMC, HMC, arag, dolo, DIAG, fluid, flux, n);

%% ================= Solve ODE for next timestep ========================
%                                       (dM/dt = input fluxes - output fluxes)


dy(i+1) = in.LMC.Ca - out.LMC.Ca;
dy(i+2) = in.LMC.Mg - out.LMC.Mg;
dy(i+3) = in.LMC.C - out.LMC.C;
dy(i+4) = in.LMC.O - out.LMC.O;
dy(i+5) = in.LMC.Sr - out.LMC.Sr;
dy(i+6) = in.LMC.Li - out.LMC.Li;

dy(i+7) = in.HMC.Ca - out.HMC.Ca;
dy(i+8) = in.HMC.Mg - out.HMC.Mg;
dy(i+9) = in.HMC.C - out.HMC.C;
dy(i+10) = in.HMC.O - out.HMC.O;
dy(i+11) = in.HMC.Sr - out.HMC.Sr;
dy(i+12) = in.HMC.Li - out.HMC.Li;

dy(i+13) = in.arag.Ca - out.arag.Ca;
dy(i+14) = in.arag.Mg - out.arag.Mg;
dy(i+15) = in.arag.C - out.arag.C;
dy(i+16) = in.arag.O - out.arag.O;
dy(i+17) = in.arag.Sr - out.arag.Sr;
dy(i+18) = in.arag.Li - out.arag.Li;

dy(i+19) = in.dolo.Ca - out.dolo.Ca;
dy(i+20) = in.dolo.Mg - out.dolo.Mg;
dy(i+21) = in.dolo.C - out.dolo.C;
dy(i+22) = in.dolo.O - out.dolo.O;
dy(i+23) = in.dolo.Sr - out.dolo.Sr;
dy(i+24) = in.dolo.Li - out.dolo.Li;

dy(i+25) = in.DIAG.Ca - out.DIAG.Ca;
dy(i+26) = in.DIAG.Mg - out.DIAG.Mg;
dy(i+27) = in.DIAG.C - out.DIAG.C;
dy(i+28) = in.DIAG.O - out.DIAG.O;
dy(i+29) = in.DIAG.Sr - out.DIAG.Sr;
dy(i+30) = in.DIAG.Li - out.DIAG.Li;

dy(i+31) = 0;     % Constant isotopic end-members of the initial composition of primary LMC ...
dy(i+32) = 0;
dy(i+33) = 0;
dy(i+34) = 0;
dy(i+35) = 0;
dy(i+36) = 0;

dy(i+37) = 0;     % ... and HMC ... 
dy(i+38) = 0;
dy(i+39) = 0;
dy(i+40) = 0;
dy(i+41) = 0;
dy(i+42) = 0;

dy(i+43) = 0;     % ... and arag... 
dy(i+44) = 0;
dy(i+45) = 0;
dy(i+46) = 0;
dy(i+47) = 0;
dy(i+48) = 0;

dy(i+49) = 0;     % ... and initial dolomite. 
dy(i+50) = 0;
dy(i+51) = 0;
dy(i+52) = 0;
dy(i+53) = 0;
dy(i+54) = 0;

dy(i+55) = in.DIAG.dCa - out.DIAG.dCa;
dy(i+56) = in.DIAG.dMg - out.DIAG.dMg;
dy(i+57) = in.DIAG.dC - out.DIAG.dC;
dy(i+58) = in.DIAG.dO - out.DIAG.dO;
dy(i+59) = in.DIAG.dSr - out.DIAG.dSr;
dy(i+60) = in.DIAG.dLi - out.DIAG.dLi;

dy(i+61) = in.solid.dCa - out.solid.dCa;    % note: not divided by mass
dy(i+62) = in.solid.dMg - out.solid.dMg;
dy(i+63) = in.solid.dC - out.solid.dC;
dy(i+64) = in.solid.dO - out.solid.dO;
dy(i+65) = in.solid.dSr - out.solid.dSr;
dy(i+66) = in.solid.dLi - out.solid.dLi;

dy(i+67) = in.fluid.Ca - out.fluid.Ca;
dy(i+68) = in.fluid.Mg - out.fluid.Mg;
dy(i+69) = in.fluid.C - out.fluid.C;
dy(i+70) = in.fluid.O - out.fluid.O;
dy(i+71) = in.fluid.Sr - out.fluid.Sr;
dy(i+72) = in.fluid.Li - out.fluid.Li;

dy(i+73) = in.fluid.dCa - out.fluid.dCa;
dy(i+74) = in.fluid.dMg - out.fluid.dMg;
dy(i+75) = in.fluid.dC - out.fluid.dC;
dy(i+76) = in.fluid.dO - out.fluid.dO;
dy(i+77) = in.fluid.dSr - out.fluid.dSr;
dy(i+78) = in.fluid.dLi - out.fluid.dLi;

dy(i+79) = 0;     % Fluid flow rate kept constant

end

end
