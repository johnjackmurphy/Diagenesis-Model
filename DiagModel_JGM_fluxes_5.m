% =============================================================
%%                                                    FLUXES
%  
%                Function that solve all input and output fluxes to each box at timestep t  
% =============================================================

function [in, out] = DiagModel_JGM_fluxes_5(LMC, HMC, arag, dolo, DIAG, fluid, flux, n)

global alpha R M tune

%% ==================== Reaction Rates ==============================

% Reaction rates are now defined in the BASE script for each mineral (1/yr):

%% ==================== Stoichiometry ==============================

% Stoichiometry of the diagenetic mineral:
MgC = M;                                % Mg/C for dolomite M = 0.5; for low-Mg calcite M = 0.035
CaC = (1 - MgC);                       % Ca/C is inverse of Mg
Kd.Sr = tune.L;                           % Distribution coeffient (Kd = solid/fluid) for Sr/Ca into diagenetic mineral
Kd.Li = tune.J;                            % Distribution coeffient (Kd = solid/fluid) for Li/Ca into diagenetic mineral



%% ==================== Fluxes into Diagenetic Mineral  ===================

%     To conserve moles of carbon in the sediment (i.e. all C that dissolves is re-precipitated
%     within the same box), the precipitation of diagenetic minerals is scaled directly to the 
%     dissolution rate of carbon in the primary minerals using the stoichiometry of the diagenetic
%     mineral, e.g. dolomite CaMg(CO3)2
%
%     note: 
%     * Elements other than carbon need not be conserved within each box. 
%     * Trace elements are normalized to the Ca flux 


% LMC to Diagenetic Mineral
% i.e. flux of diagenetic material that came from LMC.
dia.LMC.C = R.LMC*LMC.C;       % (mol/yr): all C dissolved is reprecipitated; determines the diag. min. flux
dia.LMC.Mg = MgC*dia.LMC.C;  % (mol/yr): from stochiometry of diag. min. to match C flux. 
dia.LMC.Ca = CaC*dia.LMC.C;   % (mol/yr): from stochiometry of diag. min. to match C flux. 
dia.LMC.O = 3*dia.LMC.C;        % (mol/yr): from stochiometry of diag. min. to match C flux. 
dia.LMC.Sr = R.LMC*LMC.Ca*Kd.Sr*fluid(n+1).Sr/fluid(n+1).Ca;    %  Kd = Solid(Sr/Ca) / Fluid(Sr/Ca)                                                                                             
dia.LMC.Li = R.LMC*LMC.Ca*Kd.Li*fluid(n+1).Li/fluid(n+1).Ca;    %   Kd = Solid(Li/Ca) / Fluid(Li/Ca)
                                                                          %   Both Scaled to the flux of Ca dissolved from LMC 

% HMC to Diagenetic Mineral
% i.e. flux of diagenetic material that came from HMC.
dia.HMC.C = R.HMC*HMC.C;       % (mol/yr): all C dissolved is reprecipitated; determines the diag. min. flux
dia.HMC.Mg = MgC*dia.HMC.C;  % (mol/yr): from stochiometry of diag. min. to match C flux. 
dia.HMC.Ca = CaC*dia.HMC.C;   % (mol/yr): from stochiometry of diag. min. to match C flux. 
dia.HMC.O = 3*dia.HMC.C;        % (mol/yr): from stochiometry of diag. min. to match C flux. 
dia.HMC.Sr = R.HMC*HMC.Ca*Kd.Sr*fluid(n+1).Sr/fluid(n+1).Ca;    %  Kd = Solid(Sr/Ca) / Fluid(Sr/Ca)                                                                                             
dia.HMC.Li = R.HMC*HMC.Ca*Kd.Li*fluid(n+1).Li/fluid(n+1).Ca;    %   Kd = Solid(Li/Ca) / Fluid(Li/Ca)
                                                                          %   Both Scaled to the flux of Ca dissolved from HMC                                                                                          

% Aragonite to Diagenetic Mineral
% i.e. flux of diagenetic material that came from aragonite.
dia.arag.C = R.arag*arag.C;             % Same as above, but from aragonite
dia.arag.Mg = MgC*dia.arag.C;
dia.arag.Ca = CaC*dia.arag.C;
dia.arag.O = 3*dia.arag.C; 
dia.arag.Sr = Kd.Sr*R.arag*arag.Ca*fluid(n+1).Sr/fluid(n+1).Ca;
dia.arag.Li = Kd.Li*R.arag*arag.Ca*fluid(n+1).Li/fluid(n+1).Ca;

% Dolomite to Diagenetic Mineral
% i.e. flux of diagenetic material that came from dolomite.
dia.dolo.C = R.dolo*dolo.C;             % Same as above, but from dolomite
dia.dolo.Mg = MgC*dia.dolo.C;
dia.dolo.Ca = CaC*dia.dolo.C;
dia.dolo.O = 3*dia.dolo.C;     
dia.dolo.Sr = R.dolo*Kd.Sr*dolo.Ca*fluid(n+1).Sr/fluid(n+1).Ca;
dia.dolo.Li = R.dolo*Kd.Li*dolo.Ca*fluid(n+1).Li/fluid(n+1).Ca;

% Diagenetic Mineral to Diagenetic Mineral
% i.e. flux of diagenetic material that came from previous diagenetic mineral.
dia.DIAG.C = R.DIAG*DIAG.C;             % Same as above, but from diagenetic mineral
dia.DIAG.Mg = MgC*dia.DIAG.C;
dia.DIAG.Ca = CaC*dia.DIAG.C;
dia.DIAG.O = 3*dia.DIAG.C;     
dia.DIAG.Sr = R.DIAG*Kd.Sr*DIAG.Ca*fluid(n+1).Sr/fluid(n+1).Ca;
dia.DIAG.Li = R.DIAG*Kd.Li*DIAG.Ca*fluid(n+1).Li/fluid(n+1).Ca;


%% =============  Sediment input and output fluxes (F_in - F_out)  ==============
%                                   Fluxes are set by the rate constant (R) 

% PRIMARY LMC
in.LMC.Ca = 0;                         % Precipitation of primary LMC (mol/yr)
in.LMC.Mg = 0; 
in.LMC.C = 0;
in.LMC.O = 0;
in.LMC.Sr = 0;
in.LMC.Li = 0;

out.LMC.Ca = R.LMC*LMC.Ca;      % Dissolution of primary LMC (mol/yr)
out.LMC.Mg = R.LMC*LMC.Mg;
out.LMC.C = R.LMC*LMC.C;
out.LMC.O = R.LMC*LMC.O;
out.LMC.Sr = R.LMC*LMC.Sr;
out.LMC.Li = R.LMC*LMC.Li;


% PRIMARY HMC
in.HMC.Ca = 0;                         % Precipitation of primary HMC (mol/yr)
in.HMC.Mg = 0; 
in.HMC.C = 0;
in.HMC.O = 0;
in.HMC.Sr = 0;
in.HMC.Li = 0;

out.HMC.Ca = R.HMC*HMC.Ca;       % Dissolution of primary HMC (mol/yr)
out.HMC.Mg = R.HMC*HMC.Mg;
out.HMC.C = R.HMC*HMC.C;
out.HMC.O = R.HMC*HMC.O;
out.HMC.Sr = R.HMC*HMC.Sr;
out.HMC.Li = R.HMC*HMC.Li;

% PRIMARY ARAGONITE
in.arag.Ca  = 0;                         % Precipitation of primary aragonite (mol/yr)
in.arag.Mg = 0;  
in.arag.C   = 0;
in.arag.O   = 0;
in.arag.Sr  = 0;
in.arag.Li  = 0;

out.arag.Ca = R.arag*arag.Ca;       % Dissolution of primary aragonite (mol/yr)                    
out.arag.Mg = R.arag*arag.Mg;
out.arag.C = R.arag*arag.C;
out.arag.O = R.arag*arag.O;
out.arag.Sr = R.arag*arag.Sr;
out.arag.Li = R.arag*arag.Li;

% PRIMARY DOLOMITE
in.dolo.Ca  = 0;                         % Precipitation of primary dolomite (mol/yr)
in.dolo.Mg = 0;  
in.dolo.C   = 0;
in.dolo.O   = 0;
in.dolo.Sr  = 0;
in.dolo.Li  = 0;

out.dolo.Ca = R.dolo*dolo.Ca;       % Dissolution of primary dolomite (mol/yr)                    
out.dolo.Mg = R.dolo*dolo.Mg;
out.dolo.C = R.dolo*dolo.C;
out.dolo.O = R.dolo*dolo.O;
out.dolo.Sr = R.dolo*dolo.Sr;
out.dolo.Li = R.dolo*dolo.Li;



% DIAGENETIC MINERAL
% Precipitation of the diagenetic mineral has to equal the dissolution of
% the solid in order to approximately conserve mass of the rock. 
% Hence, the dissolution flux of primary minerals is equal to the 
% precipitation flux of diagenetic mineral, but scaled stoichimetric (by M).

in.DIAG.Ca = dia.LMC.Ca + dia.HMC.Ca + dia.arag.Ca + dia.dolo.Ca + dia.DIAG.Ca;       % Precipitation of diagenetic mineral (mol/yr)
in.DIAG.Mg = dia.LMC.Mg + dia.HMC.Mg + dia.arag.Mg + dia.dolo.Mg + dia.DIAG.Mg;
in.DIAG.C = dia.LMC.C + dia.HMC.C + dia.arag.C + dia.dolo.C + dia.DIAG.C;
in.DIAG.O = dia.LMC.O + dia.HMC.O + dia.arag.O + dia.dolo.O + dia.DIAG.O;
in.DIAG.Sr = dia.LMC.Sr + dia.HMC.Sr + dia.arag.Sr + dia.dolo.Sr + dia.DIAG.Sr;
in.DIAG.Li = dia.LMC.Li + dia.HMC.Li + dia.arag.Li + dia.dolo.Li + dia.DIAG.Li;

out.DIAG.Ca = R.DIAG*DIAG.Ca;                                       % Dissolution of diagenetic mineral  (mol/yr)
out.DIAG.Mg = R.DIAG*DIAG.Mg;
out.DIAG.C = R.DIAG*DIAG.C;
out.DIAG.O = R.DIAG*DIAG.O;
out.DIAG.Sr = R.DIAG*DIAG.Sr;
out.DIAG.Li = R.DIAG*DIAG.Li;


%% =============  Fluid input and output fluxes (F_in - F_out)  =================
%
%                          The composition of the fluid is controlled by initial seawater flushing 
%                          through the sediment, the continous dissolution of primary minerals, 
%                          and re-precipitation of diagenetic minerals. 
%
%                           fluid.Ca(n)     refers to the composition of the previous box 
%                           fluid.Ca(n+1)  refers to the updated fluid composition.
%
%                           in.fluid.X      =   X in fluid flowing in    +  X dissolved from within the box
%                           out.fluid.X    =   X in fluid flowing out  +  X precipitated into diagenetic mineral

% Calcium in the fluid:
in.fluid.Ca   = flux*fluid(n).Ca + R.LMC*LMC.Ca + R.HMC*HMC.Ca + R.arag*arag.Ca + R.dolo*dolo.Ca + R.DIAG*DIAG.Ca;  
out.fluid.Ca = flux*fluid(n+1).Ca + dia.LMC.Ca +  dia.HMC.Ca + dia.arag.Ca + dia.dolo.Ca + dia.DIAG.Ca;

% Magnesium in the fluid:
in.fluid.Mg = flux*fluid(n).Mg + R.LMC*LMC.Mg + R.HMC*HMC.Mg + R.arag*arag.Mg + R.dolo*dolo.Mg + R.DIAG*DIAG.Mg;
out.fluid.Mg = flux*fluid(n+1).Mg + dia.LMC.Mg + dia.HMC.Mg + dia.arag.Mg + dia.dolo.Mg + dia.DIAG.Mg;

% Carbon in the fluid:
in.fluid.C = flux*fluid(n).C + R.LMC*LMC.C + R.HMC*HMC.C + R.arag*arag.C + R.dolo*dolo.C + R.DIAG*DIAG.C;
out.fluid.C = flux*fluid(n+1).C + dia.LMC.C + dia.HMC.C + dia.arag.C + dia.dolo.C + dia.DIAG.C;

% Oxygen in the fluid:
in.fluid.O = flux*fluid(n).O + R.LMC*LMC.O + R.HMC*HMC.O + R.arag*arag.O + R.dolo*dolo.O + R.DIAG*DIAG.O;
out.fluid.O = flux*fluid(n+1).O + dia.LMC.O + dia.HMC.O + dia.arag.O + dia.dolo.O + dia.DIAG.O;

% Strontium in the fluid:
in.fluid.Sr = flux*fluid(n).Sr + R.LMC*LMC.Sr + R.HMC*HMC.Sr + R.arag*arag.Sr + R.dolo*dolo.Sr + R.DIAG*DIAG.Sr;
out.fluid.Sr = flux*fluid(n+1).Sr + dia.LMC.Sr + dia.HMC.Sr + dia.arag.Sr + dia.dolo.Sr + dia.DIAG.Sr;

% Lithium in the fluid:
in.fluid.Li = flux*fluid(n).Li + R.LMC*LMC.Li + R.HMC*HMC.Li + R.arag*arag.Li + R.dolo*dolo.Li + R.DIAG*DIAG.Li;
out.fluid.Li = flux*fluid(n+1).Li + dia.LMC.Li + dia.HMC.Li + dia.arag.Li + dia.dolo.Li + dia.DIAG.Li;


%% ==================  Isotopic Composition  ==========================
%
%                   Mass [mol] of the individual isotopes in the solid is not constant so
%                   dM*delta/dt has to be solved without isolating delta! This means that 
%                   all delta values that are solved for have units of mol*delta and has to be
%                   divided by mol to get delta. This is accounted for in all equations so that
%                   the units of flux is always mol/yr.


% ISOTOPIC COMPOSITION OF PRECIPITATING DIAGENETIC MINERAL
%   using epsilon = 1000 ln(alpha)

dia.d.Ca = fluid(n+1).dCa + 1e3*log(alpha.Ca);    % dX +  epsilon
dia.d.Mg = fluid(n+1).dMg + 1e3*log(alpha.Mg);
dia.d.C = fluid(n+1).dC + 1e3*log(alpha.C);
dia.d.O = fluid(n+1).dO + 1e3*log(alpha.O);
dia.d.Sr = fluid(n+1).dSr + 1e3*log(alpha.Sr);
dia.d.Li = fluid(n+1).dLi + 1e3*log(alpha.Li);


% ISOTOPE FLUXES OF DIAGENETIC MINERAL 
%   Only changing in the case where the diagenetic mineral is
%   being re-dissolved by setting R.DIAG > 0

in.DIAG.dCa = (dia.LMC.Ca + dia.HMC.Ca + dia.arag.Ca + dia.dolo.Ca + dia.DIAG.Ca)*dia.d.Ca;       % (mol/yr)*dX
out.DIAG.dCa = R.DIAG*DIAG.Ca*DIAG.dCa;                                                                     %  mol*(1/yr)*dX

in.DIAG.dMg = (dia.LMC.Mg + dia.HMC.Mg + dia.arag.Mg + dia.dolo.Mg + dia.DIAG.Mg)*dia.d.Mg;
out.DIAG.dMg = R.DIAG*DIAG.Mg*DIAG.dMg;

in.DIAG.dC = (dia.LMC.C + dia.HMC.C + dia.arag.C + dia.dolo.C  + dia.DIAG.C)*dia.d.C;
out.DIAG.dC = R.DIAG*DIAG.C*DIAG.dC;

in.DIAG.dO = (dia.LMC.O + dia.HMC.O + dia.arag.O + dia.dolo.O + dia.DIAG.O)*dia.d.O;
out.DIAG.dO = R.DIAG*DIAG.O*DIAG.dO;

in.DIAG.dSr = (dia.LMC.Sr + dia.HMC.Sr + dia.arag.Sr + dia.dolo.Sr + dia.DIAG.Sr)*dia.d.Sr;
out.DIAG.dSr = R.DIAG*DIAG.Sr*DIAG.dSr;

in.DIAG.dLi = (dia.LMC.Li + dia.HMC.Li + dia.arag.Li + dia.dolo.Li + dia.DIAG.Li)*dia.d.Li;
out.DIAG.dLi = R.DIAG*DIAG.Li*DIAG.dLi;


% CALCIUM ISOTOPE FLUXES:

in.solid.dCa = (dia.LMC.Ca + dia.HMC.Ca + dia.arag.Ca + dia.dolo.Ca + dia.DIAG.Ca)*dia.d.Ca;       %  (mol/yr)*dX
out.solid.dCa = R.LMC*LMC.Ca*LMC.dCa +  R.HMC*HMC.Ca*HMC.dCa + R.arag*arag.Ca*arag.dCa ... %  (1/yr)*mol*dX
                   + R.dolo*dolo.Ca*dolo.dCa + R.DIAG*DIAG.Ca*DIAG.dCa;
in.fluid.dCa = flux*fluid(n).Ca*fluid(n).dCa + R.LMC*LMC.Ca*LMC.dCa + R.HMC*HMC.Ca*HMC.dCa ...                                                                        % (1/yr)*mol*dX
                 + R.arag*arag.Ca*arag.dCa + R.dolo*dolo.Ca*dolo.dCa + R.DIAG*DIAG.Ca*DIAG.dCa;   % (1/yr)*mol*dX
out.fluid.dCa = flux*fluid(n+1).Ca*fluid(n+1).dCa +...                                                         % (1/yr)*mol*dX
                       (dia.LMC.Ca + dia.HMC.Ca + dia.arag.Ca + dia.dolo.Ca + dia.DIAG.Ca)*dia.d.Ca;                    % (mol/yr)*dX


% MAGNESIUM ISOTOPE FLUXES:

in.solid.dMg = (dia.LMC.Mg + dia.HMC.Mg + dia.arag.Mg + dia.dolo.Mg + dia.DIAG.Mg)*dia.d.Mg;       
out.solid.dMg = R.LMC*LMC.Mg*LMC.dMg +  R.HMC*HMC.Mg*HMC.dMg + R.arag*arag.Mg*arag.dMg ... 
                   + R.dolo*dolo.Mg*dolo.dMg + R.DIAG*DIAG.Mg*DIAG.dMg;
in.fluid.dMg = flux*fluid(n).Mg*fluid(n).dMg + R.LMC*LMC.Mg*LMC.dMg + R.HMC*HMC.Mg*HMC.dMg ...                                                                        
                 + R.arag*arag.Mg*arag.dMg + R.dolo*dolo.Mg*dolo.dMg + R.DIAG*DIAG.Mg*DIAG.dMg;   
out.fluid.dMg = flux*fluid(n+1).Mg*fluid(n+1).dMg +...                                                         
                       (dia.LMC.Mg + dia.HMC.Mg + dia.arag.Mg + dia.dolo.Mg + dia.DIAG.Mg)*dia.d.Mg;    


                    
% CARBON ISOTOPE FLUXES:

in.solid.dC = (dia.LMC.C + dia.HMC.C + dia.arag.C + dia.dolo.C + dia.DIAG.C)*dia.d.C;       
out.solid.dC = R.LMC*LMC.C*LMC.dC +  R.HMC*HMC.C*HMC.dC + R.arag*arag.C*arag.dC ... 
                   + R.dolo*dolo.C*dolo.dC + R.DIAG*DIAG.C*DIAG.dC;
in.fluid.dC = flux*fluid(n).C*fluid(n).dC + R.LMC*LMC.C*LMC.dC + R.HMC*HMC.C*HMC.dC ...                                                                        
                 + R.arag*arag.C*arag.dC + R.dolo*dolo.C*dolo.dC + R.DIAG*DIAG.C*DIAG.dC;   
out.fluid.dC = flux*fluid(n+1).C*fluid(n+1).dC +...                                                         
                       (dia.LMC.C + dia.HMC.C + dia.arag.C + dia.dolo.C + dia.DIAG.C)*dia.d.C;  

                  
% OXYGEN ISOTOPE FLUXES: 

in.solid.dO = (dia.LMC.O + dia.HMC.O + dia.arag.O + dia.dolo.O + dia.DIAG.O)*dia.d.O;       
out.solid.dO = R.LMC*LMC.O*LMC.dO +  R.HMC*HMC.O*HMC.dO + R.arag*arag.O*arag.dO ... 
                   + R.dolo*dolo.O*dolo.dO + R.DIAG*DIAG.O*DIAG.dO;
in.fluid.dO = flux*fluid(n).O*fluid(n).dO + R.LMC*LMC.O*LMC.dO + R.HMC*HMC.O*HMC.dO ...                                                                        
                 + R.arag*arag.O*arag.dO + R.dolo*dolo.O*dolo.dO + R.DIAG*DIAG.O*DIAG.dO;   
out.fluid.dO = flux*fluid(n+1).O*fluid(n+1).dO +...                                                         
                       (dia.LMC.O + dia.HMC.O + dia.arag.O + dia.dolo.O + dia.DIAG.O)*dia.d.O;  


% STRONTIUM ISOTOPE FLUXES:

in.solid.dSr = (dia.LMC.Sr + dia.HMC.Sr + dia.arag.Sr + dia.dolo.Sr + dia.DIAG.Sr)*dia.d.Sr;       
out.solid.dSr = R.LMC*LMC.Sr*LMC.dSr +  R.HMC*HMC.Sr*HMC.dSr + R.arag*arag.Sr*arag.dSr ... 
                   + R.dolo*dolo.Sr*dolo.dSr + R.DIAG*DIAG.Sr*DIAG.dSr;
in.fluid.dSr = flux*fluid(n).Sr*fluid(n).dSr + R.LMC*LMC.Sr*LMC.dSr + R.HMC*HMC.Sr*HMC.dSr ...                                                                        
                 + R.arag*arag.Sr*arag.dSr + R.dolo*dolo.Sr*dolo.dSr + R.DIAG*DIAG.Sr*DIAG.dSr;   
out.fluid.dSr = flux*fluid(n+1).Sr*fluid(n+1).dSr +...                                                         
                       (dia.LMC.Sr + dia.HMC.Sr + dia.arag.Sr + dia.dolo.Sr + dia.DIAG.Sr)*dia.d.Sr;  
                  
% LITHIUM ISOTOPE FLUXES: 

in.solid.dLi = (dia.LMC.Li + dia.HMC.Li + dia.arag.Li + dia.dolo.Li + dia.DIAG.Li)*dia.d.Li;       
out.solid.dLi = R.LMC*LMC.Li*LMC.dLi +  R.HMC*HMC.Li*HMC.dLi + R.arag*arag.Li*arag.dLi ... 
                   + R.dolo*dolo.Li*dolo.dLi + R.DIAG*DIAG.Li*DIAG.dLi;
in.fluid.dLi = flux*fluid(n).Li*fluid(n).dLi + R.LMC*LMC.Li*LMC.dLi + R.HMC*HMC.Li*HMC.dLi ...                                                                        
                 + R.arag*arag.Li*arag.dLi + R.dolo*dolo.Li*dolo.dLi + R.DIAG*DIAG.Li*DIAG.dLi;   
out.fluid.dLi = flux*fluid(n+1).Li*fluid(n+1).dLi +...                                                         
                       (dia.LMC.Li + dia.HMC.Li + dia.arag.Li + dia.dolo.Li + dia.DIAG.Li)*dia.d.Li;  


end
