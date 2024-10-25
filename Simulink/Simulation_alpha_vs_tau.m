% LATAM simulation: Alpha vs Tau verification
% Script to test simulation in a 6 Zones DC MG for article


%General Configurations
%cvx_setup
%cvx_solver mosek
%%
%                                  SEcond level
%
clc 
clear

V_ref   = 1050;
tau = 0.066;
h = 0.99; 
mu = h; 


time_1 = 0.5; % Time delay active + P_min
time_2 = 5; % P_max / V_min
time_3 = 2;
time_a_up = 100; %positive step on a
time_4 = 15; % first ptot perturbation + 
time_a_down = 100;
time_5 = 100;  % change in 'a' 
time_6 = 13;  % Tau drop change  
time_tau_down = 25;
time_7 = 100;  % drop 90% Ptot back to stablity
time_8 = 100 ; % change in PI controllers gain
time_9 = 100;
time_10 = 35; % second ptot perturbation
time_R_up = 5;
       
g12    = 1.2           ;
g13    = ((0.383+1.756) / (0.383*1.756))^(-1);
g14    = 0             ;
g23    = 0             ;
g24    = 1.55          ;
g34    = 0.75          ;

% Interconnection Graph    
        aij = 1;
        
% Control Gains
        % Voltage I controller   
            k_iv  =   1               ;

        % Current PI controller   
            k_ii  =   60              ;
            k_pi  =   25           ;  
            I_rate  =   [100,50,33,33];
            
            k_ii2  =   10              ;
            k_pi2  =   0.5             ;  
            I_rate  =   [100,50,33,33];
            

            
% Defining the initial alpha value
alpha = 0.85;
alpha2 = 0.85;


P_tot = 225e3;
%dist = [0.1,0.3,0.45,0.15]; %dist nº1
%dist = [0.45,0.3,0.1,0.15]; %dist nº2
dist = [0.3,0.45,0.15,0.1]; %dist nº3
P_MG = P_tot * dist;





Perturbation_gain_1 = 1;
Perturbation_gain_2 = 1 ;
Perturbation_gain_3 = 1;
Perturbation_gain_4 = 1;
        
%%
%                                   CLUSTER Nï¿½1
%
% Parameter Definition

    % General Parameters
        
        n       = 3 ; %N of buses or Zones
        dist    = [0.4,0.5,0.1];
        
    % Passive filter components
        Ctot        = 1.7e-3;
        C_mg_1    	= Ctot .* dist;
        
        Ltot        = 0.75e-3;
        L_mg_1      = ones(n) .* n .* Ltot;
        
        rtot        = 0.25;
        rd_mg_1     = n^2 * rtot .* dist;  

    % Interconnection Impedances
        r_jk    = 0.1             ;
       
        g12_1     = 1/r_jk        ;
        g13_1     = 1/r_jk        ;
        g23_1     = 1/r_jk        ;
        
    % Interconnection Graph    
        aij_1 = 10;

    % Control Gains
        % Voltage I controller   
            k_iv_1  =   30               ;

        % Current PI controller   
            k_ii_1  =   1000            ;
            
            
            k_pi_1  =   10           ;
    
    % Droop Controller
        I_rate_1    =   [5, 2.5, 2.5];
        
        r_d0    =   1.3                 ;
        kd_1      =   [5, 2.5, 2.5]        ;
        rd1    =   r_d0/kd_1(1)          ;
        rd2     =  r_d0/kd_1(2)          ;
        rd3    =   r_d0/kd_1(3)          ;

    % Time Delay
        %tau     =   0.01        ;
        %mu      =   0.1         ; h= mu;  
        status2008 = 'Solved';


% ---------------------------------------------
%   MATRIX POLYTOPIC GENERATION
% ---------------------------------------------

    vertex      = 2^n;
    P_CPL_max   = (1-alpha) * P_MG(1);
    P_CPL_min   = (1-alpha) * 10e3;
    P_DER_max   = 1e3;
    P_DER_min   = 1e3;
    P_Rload_max = (alpha) * P_MG(1);
    P_Rload_min = (alpha) * 10e3;
    
    %dist = [1, 0.6, 0.25 , 0.3];
    %dist = [0.6,1,0.3,0.25]
    P_CPL_sup   = P_CPL_max*dist; %ESS charging and collapse bus 1 and 2 because they are RES and share bus connections
    P_CPL_inf   = P_CPL_min*dist; % Also collapse bus 4 and 5 as they have similar waveforms for being residential loads
    P_DER_sup   = P_DER_max*[1      , 0   ,  0  ];
    P_DER_inf   = P_DER_min*[0      , 0   ,  0  ];
    P_Rload_sup   = P_Rload_max*dist;
    P_Rload_inf   = P_Rload_min*dist;
    
    

    P_load_sup  = P_CPL_sup - P_DER_inf ;
    P_load_inf  = P_CPL_inf - P_DER_sup ;

    P_aux           =   [P_load_sup;P_load_inf];
    G_load_aux      =   [P_Rload_sup;P_Rload_inf]./(V_ref^2);
    P_load_vertex_1   =   [];
    G_load_vertex_1   =   [];

    for i1=1:2
            for i2=1:2
                for i3=1:2
                            P_load_vertex_1 = [   P_load_vertex_1;...
                                                P_aux(i1,1) , P_aux(i2,2) , P_aux(i3,3) ];
                            G_load_vertex_1 = [G_load_vertex_1;...
                                                G_load_aux(i1,1) , G_load_aux(i2,2) , G_load_aux(i3,3)];
                
                end
            end
    end
%%
%                                   CLUSTER Nï¿½2
%
% Parameter Definition

    % General Parameters
        
        n       = 2 ; %N of buses or Zones
        dist    = [0.4,0.6];
        
    % Passive filter components
        Ctot        = 1e-3;
        C_mg_2    	= Ctot .* dist;
        
        Ltot        = 0.9e-3;
        L_mg_2      = ones(n) .* n .* Ltot;
        
        rtot        = 0.5;
        rd_mg_2     = n^2 * rtot .* dist;  

    % Interconnection Impedances
        r_jk    = 0.1             ;
       
        g12_2     = 1/r_jk        ;
        
    % Interconnection Graph    
        aij_2 = 10;

    % Control Gains
        % Voltage I controller   
            k_iv_2  =   30               ;

        % Current PI controller   
            k_ii_2  =   1000            ;
            
            
            k_pi_2  =   10           ;
    
    % Droop Controller
        I_rate_2    =   [2.5, 5];
        
        r_d0    =   1.3                 ;
        kd_2      =   [2.5, 50]        ;
        rd1    =   r_d0/kd_2(1)          ;
        rd2     =  r_d0/kd_2(2)          ;

    % Time Delay
        %tau     =   0.01        ;
        %mu      =   0.1         ; h= mu;  
        status2008 = 'Solved';


% ---------------------------------------------
%   MATRIX POLYTOPIC GENERATION
% ---------------------------------------------

    vertex      = 2^n;
    P_CPL_max   = (1 - alpha) * P_MG(2);
    P_CPL_min   = (1 - alpha) * 10e3;
    P_DER_max   = 1e3;
    P_DER_min   = 1e3;  
    P_Rload_max = (alpha) * P_MG(2);
    P_Rload_min = (alpha) * 10e3;
    
    %dist = [1, 0.6, 0.25 , 0.3];
    %dist = [0.6,1,0.3,0.25]
    P_CPL_sup   = P_CPL_max*dist; %ESS charging and collapse bus 1 and 2 because they are RES and share bus connections
    P_CPL_inf   = P_CPL_min*dist; % Also collapse bus 4 and 5 as they have similar waveforms for being residential loads
    P_DER_sup   = P_DER_max*[1      , 0   ];
    P_DER_inf   = P_DER_min*[0      , 0   ];
    P_Rload_sup   = P_Rload_max*dist;
    P_Rload_inf   = P_Rload_min*dist;

    P_load_sup  = P_CPL_sup - P_DER_inf ;
    P_load_inf  = P_CPL_inf - P_DER_sup ;

    P_aux           =   [P_load_sup;P_load_inf];
    G_load_aux      =   [P_Rload_sup;P_Rload_inf]./(V_ref^2);
    P_load_vertex_2   =   [];
    G_load_vertex_2   =   [];


    for i1=1:2
            for i2=1:2
                            P_load_vertex_2 = [   P_load_vertex_2;...
                                                P_aux(i1,1) , P_aux(i2,2) ]; 
                            G_load_vertex_2 = [G_load_vertex_2;...
                                                G_load_aux(i1,1) , G_load_aux(i2,2) ];
            end
    end
%%
%                                   CLUSTER Nï¿½3
%
% Parameter Definition

    % General Parameters
        
        n       = 4 ; %N of buses or Zones
        dist    = [0.25,0.4,0.15,0.2];
        
    % Passive filter components
        Ctot        = 1e-3;
        C_mg_3    	= Ctot .* dist;
        
        Ltot        = 1.25e-3;
        L_mg_3      = ones(n) .* n .* Ltot;
        
        rtot        = 0.25;
        rd_mg_3     = n^2 * rtot .* dist;  

    % Interconnection Impedances
        r_jk    = 0.1             ;
       
        g12_3     = 1/r_jk           ;
        g13_3     = 1/r_jk        ;
        g14_3     = 0             ;
        g23_3     = 0             ;
        g24_3     = 1/r_jk          ;
        g34_3     = 1/r_jk          ;
        
    % Interconnection Graph    
        aij_3 =  10;
    % Control Gains
        % Voltage I controller   
            k_iv_3  =   30               ;

        % Current PI controller   
            k_ii_3  =   1000            ;
            
            
            k_pi_3  =   10           ;
    
    % Droop Controller
        I_rate_3    =   [4, 3, 2, 1];
        
        r_d0    =   1.3                 ;
        kd_3      =   [4, 3, 2, 1]        ;
        rd1    =   r_d0/kd_3(1)          ;
        rd2     =  r_d0/kd_3(2)          ;
        rd3    =   r_d0/kd_3(3)          ;
        rd4     =  r_d0/kd_3(4)          ;

    % Time Delay
        %tau     =   0.01        ;
        %mu      =   0.1         ; h= mu;  
        status2008 = 'Solved';


% ---------------------------------------------
%   MATRIX POLYTOPIC GENERATION
% ---------------------------------------------

    vertex      = 2^n;
    P_CPL_max   = (1- alpha) * P_MG(3);
    P_CPL_min   = (1-alpha) * 10e3;
    P_DER_max   = 1e3;
    P_DER_min   = 1e3;
    P_Rload_max = (alpha) * P_MG(3);
    P_Rload_min = (alpha) * 10e3;
    
    %dist = [1, 0.6, 0.25 , 0.3];
    %dist = [0.6,1,0.3,0.25]
    P_CPL_sup   = P_CPL_max*dist; %ESS charging and collapse bus 1 and 2 because they are RES and share bus connections
    P_CPL_inf   = P_CPL_min*dist; % Also collapse bus 4 and 5 as they have similar waveforms for being residential loads
    P_DER_sup   = P_DER_max*[1      , 0   ,  0  , 0   ];
    P_DER_inf   = P_DER_min*[0      , 0   ,  0  , 0   ];
    P_Rload_sup   = P_Rload_max*dist;
    P_Rload_inf   = P_Rload_min*dist;

    P_load_sup  = P_CPL_sup - P_DER_inf ;
    P_load_inf  = P_CPL_inf - P_DER_sup ;

    P_aux           =   [P_load_sup;P_load_inf];
    G_load_aux      =   [P_Rload_sup;P_Rload_inf]./(V_ref^2);
    P_load_vertex_3   =   [];
    G_load_vertex_3   =   [];

    for i1=1:2
            for i2=1:2
                for i3=1:2
                    for i4=1:2
                            P_load_vertex_3 = [   P_load_vertex_3;...
                                                P_aux(i1,1) , P_aux(i2,2) , P_aux(i3,3), P_aux(i4,4) ];
                                            
                            G_load_vertex_3 = [G_load_vertex_3;...
                                                G_load_aux(i1,1) , G_load_aux(i2,2) , G_load_aux(i3,3), G_load_aux(i4,4)];
                
                
                    end
                end
            end
    end

    
    
%% 
%                           CLUSTER 4
%

% Parameter Definition

    % General Parameters
        
        n       = 4 ; %N of buses or Zones
        dist    = [0.4,0.2,0.3,0.1];
        
    % Passive filter components
        Ctot        = 1e-3;
        C_mg_4    	= Ctot .* dist;
        
        Ltot        = 0.8e-3;
        L_mg_4      = ones(n) .* n .* Ltot;
        
        rtot        = 0.5;
        rd_mg_4     = n^2 * rtot .* dist;  

    % Interconnection Impedances
        r_jk    = 0.1             ;
       
        g12_4     = 1/r_jk           ;
        g13_4     = 1/r_jk        ;
        g14_4     = 0             ;
        g23_4     = 1/r_jk             ;
        g24_4     = 1/r_jk          ;
        g34_4     = 1/r_jk          ;
        
    % Interconnection Graph    
        aij_4 =  10;

    % Control Gains
        % Voltage I controller   
            k_iv_4  =   30               ;

        % Current PI controller   
            k_ii_4  =   1000            ;
            
            
            k_pi_4  =   10           ;
    
    % Droop Controller
        I_rate_4    =   [4, 3, 2, 1];
        
        r_d0    =   1.3                 ;
        kd_4      =   [4, 3, 2, 1]        ;
        rd1    =   r_d0/kd_4(1)          ;
        rd2     =  r_d0/kd_4(2)          ;
        rd3    =   r_d0/kd_4(3)          ;
        rd4     =  r_d0/kd_4(4)          ;

    % Time Delay
        %tau     =   0.01        ;
        %mu      =   0.1         ; h= mu;  
        status2008 = 'Solved';


% ---------------------------------------------
%   MATRIX POLYTOPIC GENERATION
% ---------------------------------------------

    vertex      = 2^n;
    P_CPL_max   = (1-alpha) * P_MG(4);
    P_CPL_min   = (1-alpha) * 10e3;
    P_DER_max   = 1e3;
    P_DER_min   = 1e3;    
    P_Rload_max = (alpha) * P_MG(4);
    P_Rload_min = (alpha) * 10e3;
    
    %dist = [1, 0.6, 0.25 , 0.3];
    %dist = [0.6,1,0.3,0.25]
    P_CPL_sup   = P_CPL_max*dist; %ESS charging and collapse bus 1 and 2 because they are RES and share bus connections
    P_CPL_inf   = P_CPL_min*dist; % Also collapse bus 4 and 5 as they have similar waveforms for being residential loads
    P_DER_sup   = P_DER_max*[0      , 0   ,  0  , 0   ];
    P_DER_inf   = P_DER_min*[0      , 0   ,  0  , 0   ];
    P_Rload_sup   = P_Rload_max*dist;
    P_Rload_inf   = P_Rload_min*dist;
    
    P_load_sup  = P_CPL_sup - P_DER_inf ;
    P_load_inf  = P_CPL_inf - P_DER_sup ;

    P_aux           =   [P_load_sup;P_load_inf];
    G_load_aux      =   [P_Rload_sup;P_Rload_inf]./(V_ref^2);
    P_load_vertex_4   =   [];
    G_load_vertex_4   =   [];

    for i1=1:2
            for i2=1:2
                for i3=1:2
                    for i4=1:2
                            P_load_vertex_4 = [   P_load_vertex_4;...
                                                P_aux(i1,1) , P_aux(i2,2) , P_aux(i3,3), P_aux(i4,4) ];
                            G_load_vertex_4 = [G_load_vertex_4;...
                                                G_load_aux(i1,1) , G_load_aux(i2,2) , G_load_aux(i3,3), G_load_aux(i4,4)];
                
                
                    end
                end
            end
    end