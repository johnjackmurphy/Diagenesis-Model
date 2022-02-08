% =============================================================
%%                                                    CONSTANTS 
%  
%              Function that sets up physical parameters of the model   
%              and calculates constants used to initialize the model   
% =============================================================

function [Molar, Box, Solid] = DiagModel_JGM_constants_5(Aragonite, LMC, HMC, Dolomite, Ca_fluid, Mg_fluid, C_fluid, Sr_fluid, Li_fluid, Solid, alpha)

global Fluid tune

%===================== Molar Masses ===============================
                                                 
Molar.Mg = 24.305;      % (g/mol) 
Molar.Ca  = 40.078;      
Molar.C    = 12.011;
Molar.O   = 15.999;
Molar.Sr   = 87.62;
Molar.Li   = 6.941;
 
% ================== Physical Features of Box ============================

Box.vol = 100^3;                % Box volume (cm3/m3)
Box.Phi = 0.5;                    % Porosity (w. fraction)
Box.RhoS = 1.8;                 % Density of sediment (g/cm3)
Box.RhoF = 1.0125;             % Density of fluid (g/cm3)

Box.sed = (1-Box.Phi)*Box.vol*Box.RhoS/1e3;     % Sediment mass (kg)
Box.fluid = Box.Phi*Box.vol*Box.RhoF/1e3;         % Fluid mass (kg)

% Weight fraction of fluid in the system (cf. Banner and Hanson, 1990)

Box.F = (Box.Phi*Box.RhoF)/(Box.Phi*Box.RhoF + (1 - Box.Phi)*Box.RhoS);


% =============== Initial Elemental Composition of Solid =======================

% LMC (CaCO3)

Solid.LMC.init = LMC*Box.sed;      % Initial mass of calcite (kg) 

ppm.LMC.Ca = 394000;              % ppm Ca in LMC
ppm.LMC.Mg = 5000-tune.B;        % ppm Mg in LMC           
ppm.LMC.C = 120000;               % ppm C in LMC             
ppm.LMC.O = 480000;               % ppm O in LMC            
ppm.LMC.Sr = 1000;                  % ppm Sr in LMC  
ppm.LMC.Li = tune.B;                 % ppm Li in LMC  

Solid.LMC.Ca = ppm.LMC.Ca*Solid.LMC.init/(Molar.Ca*1e3);            % moles Ca in calcite (LMC) 
Solid.LMC.Mg = ppm.LMC.Mg*Solid.LMC.init/(Molar.Mg*1e3);          % moles Mg in calcite (LMC)
Solid.LMC.C = ppm.LMC.C*Solid.LMC.init/(Molar.C*1e3);                 % moles C in calcite (LMC)
Solid.LMC.O = ppm.LMC.O*Solid.LMC.init/(Molar.O*1e3);                % moles O in calcite (LMC) 
Solid.LMC.Sr = ppm.LMC.Sr*Solid.LMC.init/(Molar.Sr*1e3);               % moles Sr in calcite (LMC) 
Solid.LMC.Li = ppm.LMC.Li*Solid.LMC.init/(Molar.Li*1e3);                % moles Li in calcite (LMC) 


% Aragonite (CaCO3)


Solid.arag.init = Aragonite*Box.sed;     % Initial mass of aragonite (kg)    

ppm.arag.Ca = 388500;           % ppm elements in aragonite               
ppm.arag.Mg = 2500-tune.A;                         
ppm.arag.C = 120000;                          
ppm.arag.O = 480000;                         
ppm.arag.Sr = 9000;
ppm.arag.Li = tune.A;
                            
Solid.arag.Ca = ppm.arag.Ca*Solid.arag.init/(Molar.Ca*1e3);     % moles element in aragonite
Solid.arag.Mg = ppm.arag.Mg*Solid.arag.init/(Molar.Mg*1e3);
Solid.arag.C = ppm.arag.C*Solid.arag.init/(Molar.C*1e3);
Solid.arag.O = ppm.arag.O*Solid.arag.init/(Molar.O*1e3);
Solid.arag.Sr = ppm.arag.Sr*Solid.arag.init/(Molar.Sr*1e3);
Solid.arag.Li = ppm.arag.Li*Solid.arag.init/(Molar.Li*1e3);


% HMC (CaCO3 + some Mg in Ca sites)

Solid.HMC.init = HMC*Box.sed;      % Initial mass of HMC (kg) 

ppm.HMC.Ca = 367000;              % ppm Ca in HMC
ppm.HMC.Mg = 30000-tune.C;      % ppm Mg in HMC           
ppm.HMC.C = 120000;               % ppm C in HMC             
ppm.HMC.O = 480000;               % ppm O in HMC            
ppm.HMC.Sr = 3000;                  % ppm Sr in HMC  
ppm.HMC.Li = tune.C;                 % ppm Li in HMC  

Solid.HMC.Ca = ppm.HMC.Ca*Solid.HMC.init/(Molar.Ca*1e3);            % moles Ca in HMC
Solid.HMC.Mg = ppm.HMC.Mg*Solid.HMC.init/(Molar.Mg*1e3);          % moles Mg in HMC
Solid.HMC.C = ppm.HMC.C*Solid.HMC.init/(Molar.C*1e3);                 % moles C in HMC 
Solid.HMC.O = ppm.HMC.O*Solid.HMC.init/(Molar.O*1e3);                % moles O in HMC
Solid.HMC.Sr = ppm.HMC.Sr*Solid.HMC.init/(Molar.Sr*1e3);               % moles Sr in HMC 
Solid.HMC.Li = ppm.HMC.Li*Solid.HMC.init/(Molar.Li*1e3);                % moles Li in HMC



% Dolomite (CaMg[CO3]2)

Solid.dolo.init = Dolomite*Box.sed;          % Initial mass of dolomite (kg)

ppm.dolo.Ca = 217250;       % ppm elements in dolomite 
ppm.dolo.Mg = 131800-tune.D;                         
ppm.dolo.C = 130300;                          
ppm.dolo.O = 520600;                         
ppm.dolo.Sr = 50;
ppm.dolo.Li = tune.D;

Solid.dolo.Ca = ppm.dolo.Ca*Solid.dolo.init/(Molar.Ca*1e3);         % moles element in dolomite
Solid.dolo.Mg = ppm.dolo.Mg*Solid.dolo.init/(Molar.Mg*1e3);
Solid.dolo.C = ppm.dolo.C*Solid.dolo.init/(Molar.C*1e3);
Solid.dolo.O = ppm.dolo.O*Solid.dolo.init/(Molar.O*1e3);
Solid.dolo.Sr = ppm.dolo.Sr*Solid.dolo.init/(Molar.Sr*1e3);
Solid.dolo.Li = ppm.dolo.Li*Solid.dolo.init/(Molar.Li*1e3);


% Diagenetic mineral (CaMg[CO3]2 or CaCO3)

Solid.DIAG.init = 0;          % Initial mass of diagenetic mineral (kg)

ppm.DIAG.Ca = 0;            % ppm elements in diagenetic mineral 
ppm.DIAG.Mg = 0;                         
ppm.DIAG.C = 0;                          
ppm.DIAG.O = 0;                         
ppm.DIAG.Sr = 0;
ppm.DIAG.Li = 0;

Solid.DIAG.Ca = 0;         % moles element in diagenetic mineral
Solid.DIAG.Mg = 0;
Solid.DIAG.C = 0;
Solid.DIAG.O = 0;
Solid.DIAG.Sr = 0;
Solid.DIAG.Li = 0;



% =============== Initial Isotopic Composition of Solid ========================

% Isotopic composition of bulk carbonate 

    %  molesXtotal*dXtotal = molesXLMC*dXLMC + molesXHMC*dXHMC + molesXarag*dXarag + molesXdolo*dXdolo 
    %           note: solving d(M*delta)/dt requires initial conditions to be in units of (mol*permil): 
    %           note: Solid.dX is in (mol*permil)

Solid.dCa = (Solid.LMC.dCa*Solid.LMC.Ca + Solid.HMC.dCa*Solid.HMC.Ca + Solid.arag.dCa*Solid.arag.Ca+ Solid.dolo.dCa*Solid.dolo.Ca);        % molesX*dX in all solid 
Solid.dMg = (Solid.LMC.dMg*Solid.LMC.Mg + Solid.HMC.dMg*Solid.HMC.Mg + Solid.arag.dMg*Solid.arag.Mg+ Solid.dolo.dMg*Solid.dolo.Mg); 
Solid.dC = (Solid.LMC.dC*Solid.LMC.C +Solid.HMC.dC*Solid.HMC.C +  Solid.arag.dC*Solid.arag.C+ Solid.dolo.dC*Solid.dolo.C);
Solid.dO = (Solid.LMC.dO*Solid.LMC.O + Solid.HMC.dO*Solid.HMC.O + Solid.arag.dO*Solid.arag.O+ Solid.dolo.dO*Solid.dolo.O);
Solid.dSr = (Solid.LMC.dSr*Solid.LMC.Sr +Solid.HMC.dSr*Solid.HMC.Sr + Solid.arag.dSr*Solid.arag.Sr+ Solid.dolo.dSr*Solid.dolo.Sr);
Solid.dLi = (Solid.LMC.dLi*Solid.LMC.Li + Solid.HMC.dLi*Solid.HMC.Li + Solid.arag.dLi*Solid.arag.Li+ Solid.dolo.dLi*Solid.dolo.Li);



% ===============  Initial composition of fluid ============================= 

% Moles X in Fluid:
Fluid.Ca = Ca_fluid*Box.fluid;                     % (mol/kg) * (kg) = mol
Fluid.Mg = Mg_fluid*Box.fluid;                   % (mol/kg) * (kg) = mol
Fluid.C = C_fluid*Box.fluid;                        % (mol/kg) * (kg) = mol
Fluid.O = 889e3*Box.fluid/(Molar.O*1e3);   % (mg O/kg H2O) (kg H2O) (mol O/mg O) = mol O 
Fluid.Sr = Sr_fluid*Box.fluid;                      % (mol/kg) * (kg) = mol
Fluid.Li = Li_fluid*Box.fluid;                      % (mol/kg) * (kg) = mol

% Initial moles*permil of fluid to solve ODE:
Fluid.dCa_init = Fluid.dCa*Fluid.Ca;             % mol*dX 
Fluid.dMg_init = Fluid.dMg*Fluid.Mg;           
Fluid.dC_init = Fluid.dC*Fluid.C;                  
Fluid.dO_init = Fluid.dO*Fluid.O;  
Fluid.dSr_init = Fluid.dSr*Fluid.Sr;
Fluid.dLi_init = Fluid.dLi*Fluid.Li;

end
