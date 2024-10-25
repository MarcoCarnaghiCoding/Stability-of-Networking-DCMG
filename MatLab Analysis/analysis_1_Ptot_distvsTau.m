% This Script is to generate the simulations for the plot of: Tau_max vs Peq

% Full version:
%               Diferent plot are going to be defined for different 
%               Power distribution among the clusters.

%README:
%       The idea is to make 2 nested loops:
%           * Outer loop: takes a value of P_eq = P_eq_i.
%           * Inner loop: 
%                        * Takes values of Tau for a constant d(tau)/dt and P_eq_i
%                        * Defines the maximum Tau tolerated to the given conditions using LMI
                    

%-----------------------------------------------------------------------------
%                   CONFIGURATION AND BASIC DEFINITIONS
%-----------------------------------------------------------------------------

cvx_setup
cvx_solver mosek
clc 
clear

%General
V_ref = 700;
n = 3;

%Passive filter
    C1 = 0.33e-3    ;
    C2 = 0.5e-3     ;
    C3 = 1e-3       ;
    
    L1 = 1e-3       ;
    L2 = 0.85e-3    ;
    L3 = 1.15e-3    ;
    
    rs1 = 0.45      ;
    rs2 = 0.55      ;
    rs3 = 0.5       ;

%Interconnection
    r_jk = 2            ;
    g12 = 1/r_jk        ;
    g13 = 1/r_jk        ;
    g23 = 1/r_jk        ;
    
    a12 = 1;
    a13 = 1;
    a23 = 1;
% 
% %load
%     P_load      = 100e3;
%     p_fraction  = [0.25 1 0.5];
%     
% 
%     g_cpl_1     = -(V_ref^2)/(P_load*p_fraction(1)); 
%     g_cpl_2     = -(V_ref^2)/(P_load*p_fraction(2));
%     g_cpl_3     = -(V_ref^2)/(P_load*p_fraction(3));
%     g_load_1    = 0;
%     g_load_2    = 0;
%     g_load_3    = 0;


% control

    k_iv_1  =   1              ;
    k_iv_2  =   k_iv_1          ;
    k_iv_3  =   k_iv_1          ;
    k_ii_1  =   10              ;
    k_ii_2  =   k_ii_1          ;
    k_ii_3  =   k_ii_1          ; 
    k_pi_1  =   0.5             ;
    k_pi_2  =   k_pi_1          ;
    k_pi_3  =   k_pi_1          ;
    I_rate  =   [10, 20, 10]    ;
    k_1     =   1/I_rate(1)     ;
    k_2     =   1/I_rate(2)     ;
    k_3     =   1/I_rate(3)     ;
    
    r_d0    =   1.3             ;
    kd      =   [1, 2, 1]       ;
    rd1     =   r_d0/kd(1)      ;
    rd2     =   r_d0/kd(2)      ;
    rd3     =   r_d0/kd(3)      ;
%Delay

    tau = 0.45;
    tau_max = tau;
    h   =   tau_max;
    mu = 0.9;
    
    %initilize register index
    aux = 10e3:5e3:300e3;
    iteration = 0;
    register = zeros(length(aux),3);
for P=10e3:5e3:300e3 
    iteration = iteration +1;
    tau = 0.00;
    mu  = 0.5 ;
    status2008 = 'Solved';
    while status2008 == 'Solved'
        
        %-----------------------------------------------------------------
        % UPDATE REGISTER       
        %-----------------------------------------------------------------
        register(iteration,:) = [P,tau,mu];
        tau = tau + 0.01;
        h   = tau; tau_max = tau;     
        
        %-----------------------------------------------------------------
        % MATRIX POLYTOPIC GENERATION        
        %-----------------------------------------------------------------
            vertex = 2^n;
            %increase P1 and P2 and keeping P3 constant
            p_fraction  =  [0.25 1 0.5];p_fraction    = p_fraction(1:n);
            P_load_sup  = P * p_fraction ;
            P_load_inf  = 0.1 * P * p_fraction ;
            
            P_load_sup(3) = 10e3;
            P_load_inf(3) = 1e3;


             p_aux = [P_load_sup;P_load_inf];
             p_load_vertex=[];

             for i=1:2
                for j=1:2
                    for k=1:2
                       p_load_vertex = [p_load_vertex; p_aux(i,1) , p_aux(j,2) , p_aux(k,3)]; 
                    end
                end
             end

            D = zeros(n*4,vertex);
            for i=1:vertex
                %params
                g_cpl_1     = -(p_load_vertex(i,1))/(V_ref^2); 
                g_cpl_2     = -(p_load_vertex(i,2))/(V_ref^2); 
                g_cpl_3     = -(p_load_vertex(i,3))/(V_ref^2); 
                g_load_1    = 0;
                g_load_2    = 0;
                g_load_3    = 0;    
                % coef
                d1= -(g12*k_1*k_ii_1*k_ii_2*k_iv_3 - g12*k_1*k_ii_1*k_ii_3*k_iv_2 + g12*k_2*k_ii_1*k_ii_2*k_iv_3 - g12*k_2*k_ii_2*k_ii_3*k_iv_1 - g13*k_1*k_ii_1*k_ii_2*k_iv_3 + g13*k_1*k_ii_1*k_ii_3*k_iv_2 + g13*k_3*k_ii_1*k_ii_3*k_iv_2 - g13*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g23*k_2*k_ii_1*k_ii_2*k_iv_3 - 2*g23*k_2*k_ii_2*k_ii_3*k_iv_1 + 2*g23*k_3*k_ii_1*k_ii_3*k_iv_2 - 2*g23*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_2*k_2*k_ii_1*k_ii_2*k_iv_3 - g_cpl_2*k_2*k_ii_2*k_ii_3*k_iv_1 + g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2 - g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1 + g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3 - g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1 + g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2 - g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g13*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_3*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_load_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n)/(3*g12*g13*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g13*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g12*g13*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g12*g23*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g23*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g13*g23*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g23*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g13*g23*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g13*g23*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g12*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g13*g_cpl_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g12*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g13*g_cpl_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g13*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g23*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g23*g_cpl_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_cpl_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g23*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g12*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g13*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g12*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g13*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g13*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g23*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g23*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g23*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_2*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_2*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_3*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_2*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_3*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_load_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_load_1*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_load_2*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g13*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_2*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_3*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_load_1*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n);
                d2= -(g12*k_1*k_ii_1*k_ii_2*k_iv_3 - g12*k_1*k_ii_1*k_ii_3*k_iv_2 + g12*k_2*k_ii_1*k_ii_2*k_iv_3 - g12*k_2*k_ii_2*k_ii_3*k_iv_1 + 2*g13*k_1*k_ii_1*k_ii_2*k_iv_3 - 2*g13*k_1*k_ii_1*k_ii_3*k_iv_2 - 2*g13*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*k_3*k_ii_2*k_ii_3*k_iv_1 - g23*k_2*k_ii_1*k_ii_2*k_iv_3 + g23*k_2*k_ii_2*k_ii_3*k_iv_1 - g23*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*k_1*k_ii_1*k_ii_2*k_iv_3 - g_cpl_1*k_1*k_ii_1*k_ii_3*k_iv_2 - g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1 + g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3 - g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2 - g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2 + g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g13*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_3*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_load_1*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n)/(3*g12*g13*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g13*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g12*g13*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g12*g23*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g23*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g13*g23*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g23*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g13*g23*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g13*g23*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g12*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g13*g_cpl_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g12*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g13*g_cpl_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g13*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g23*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g23*g_cpl_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_cpl_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g23*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g12*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g13*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g12*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g13*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g13*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g23*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g23*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g23*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_2*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_2*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_3*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_2*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_3*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_load_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_load_1*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_load_2*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g13*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_2*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_3*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_load_1*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n);
                d3= -(2*g12*k_1*k_ii_1*k_ii_3*k_iv_2 - 2*g12*k_1*k_ii_1*k_ii_2*k_iv_3 - 2*g12*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*k_2*k_ii_2*k_ii_3*k_iv_1 - g13*k_1*k_ii_1*k_ii_2*k_iv_3 + g13*k_1*k_ii_1*k_ii_3*k_iv_2 + g13*k_3*k_ii_1*k_ii_3*k_iv_2 - g13*k_3*k_ii_2*k_ii_3*k_iv_1 - g23*k_2*k_ii_1*k_ii_2*k_iv_3 + g23*k_2*k_ii_2*k_ii_3*k_iv_1 - g23*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*k_3*k_ii_2*k_ii_3*k_iv_1 - g_cpl_1*k_1*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*k_1*k_ii_1*k_ii_3*k_iv_2 - g_cpl_2*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_2*k_2*k_ii_2*k_ii_3*k_iv_1 - g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3 + g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2 - g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3 + g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1 + g12*g13*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g_load_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n)/(3*g12*g13*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g13*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g12*g13*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g12*g23*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g23*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g13*g23*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g23*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g13*g23*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g13*g23*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g12*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g13*g_cpl_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g12*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g13*g_cpl_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g13*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g23*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g23*g_cpl_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_cpl_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g23*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g12*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g13*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g12*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g13*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g13*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g23*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g23*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g23*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_2*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_2*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_3*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_2*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_3*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_load_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_load_1*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_load_2*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g13*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_2*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_3*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_load_1*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n);
                d4= -(g12*k_1*k_ii_1*k_ii_2*k_iv_3 - g12*k_1*k_ii_1*k_ii_3*k_iv_2 + g12*k_2*k_ii_1*k_ii_2*k_iv_3 - g12*k_2*k_ii_2*k_ii_3*k_iv_1 - g13*k_1*k_ii_1*k_ii_2*k_iv_3 + g13*k_1*k_ii_1*k_ii_3*k_iv_2 + g13*k_3*k_ii_1*k_ii_3*k_iv_2 - g13*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g23*k_2*k_ii_1*k_ii_2*k_iv_3 - 2*g23*k_2*k_ii_2*k_ii_3*k_iv_1 + 2*g23*k_3*k_ii_1*k_ii_3*k_iv_2 - 2*g23*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_2*k_2*k_ii_1*k_ii_2*k_iv_3 - g_cpl_2*k_2*k_ii_2*k_ii_3*k_iv_1 + g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2 - g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1 + g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3 - g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1 + g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2 - g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g12*g13*k_2*k_ii_1*k_ii_2*k_iv_3*rd1 - 3*g12*g13*k_2*k_ii_2*k_ii_3*k_iv_1*rd1 + 3*g12*g13*k_3*k_ii_1*k_ii_3*k_iv_2*rd1 - 3*g12*g13*k_3*k_ii_2*k_ii_3*k_iv_1*rd1 + 3*g12*g23*k_2*k_ii_1*k_ii_2*k_iv_3*rd1 - 3*g12*g23*k_2*k_ii_2*k_ii_3*k_iv_1*rd1 + 3*g12*g23*k_3*k_ii_1*k_ii_3*k_iv_2*rd1 - 3*g12*g23*k_3*k_ii_2*k_ii_3*k_iv_1*rd1 + 3*g13*g23*k_2*k_ii_1*k_ii_2*k_iv_3*rd1 - 3*g13*g23*k_2*k_ii_2*k_ii_3*k_iv_1*rd1 + 3*g13*g23*k_3*k_ii_1*k_ii_3*k_iv_2*rd1 - 3*g13*g23*k_3*k_ii_2*k_ii_3*k_iv_1*rd1 + 3*g12*g13*k_2*k_ii_1*k_ii_2*k_iv_3*rs1 - 3*g12*g13*k_2*k_ii_2*k_ii_3*k_iv_1*rs1 + 3*g12*g13*k_3*k_ii_1*k_ii_3*k_iv_2*rs1 - 3*g12*g13*k_3*k_ii_2*k_ii_3*k_iv_1*rs1 + 3*g12*g23*k_2*k_ii_1*k_ii_2*k_iv_3*rs1 - 3*g12*g23*k_2*k_ii_2*k_ii_3*k_iv_1*rs1 + 3*g12*g23*k_3*k_ii_1*k_ii_3*k_iv_2*rs1 - 3*g12*g23*k_3*k_ii_2*k_ii_3*k_iv_1*rs1 + 3*g13*g23*k_2*k_ii_1*k_ii_2*k_iv_3*rs1 - 3*g13*g23*k_2*k_ii_2*k_ii_3*k_iv_1*rs1 + 3*g13*g23*k_3*k_ii_1*k_ii_3*k_iv_2*rs1 - 3*g13*g23*k_3*k_ii_2*k_ii_3*k_iv_1*rs1 + g12*g_cpl_1*k_2*k_ii_1*k_ii_2*k_iv_3*rd1 - g12*g_cpl_1*k_2*k_ii_2*k_ii_3*k_iv_1*rd1 + g12*g_cpl_2*k_2*k_ii_1*k_ii_2*k_iv_3*rd1 - g12*g_cpl_2*k_2*k_ii_2*k_ii_3*k_iv_1*rd1 + g13*g_cpl_1*k_3*k_ii_1*k_ii_3*k_iv_2*rd1 - g13*g_cpl_1*k_3*k_ii_2*k_ii_3*k_iv_1*rd1 + 2*g13*g_cpl_2*k_2*k_ii_1*k_ii_2*k_iv_3*rd1 - 2*g13*g_cpl_2*k_2*k_ii_2*k_ii_3*k_iv_1*rd1 + 2*g12*g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2*rd1 - 2*g12*g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1*rd1 + g13*g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2*rd1 - g13*g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1*rd1 + 2*g23*g_cpl_1*k_2*k_ii_1*k_ii_2*k_iv_3*rd1 - 2*g23*g_cpl_1*k_2*k_ii_2*k_ii_3*k_iv_1*rd1 + 2*g23*g_cpl_1*k_3*k_ii_1*k_ii_3*k_iv_2*rd1 - 2*g23*g_cpl_1*k_3*k_ii_2*k_ii_3*k_iv_1*rd1 + g12*g_load_1*k_2*k_ii_1*k_ii_2*k_iv_3*rd1 - g12*g_load_1*k_2*k_ii_2*k_ii_3*k_iv_1*rd1 + g12*g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3*rd1 - g12*g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1*rd1 + g13*g_load_1*k_3*k_ii_1*k_ii_3*k_iv_2*rd1 - g13*g_load_1*k_3*k_ii_2*k_ii_3*k_iv_1*rd1 + 2*g13*g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3*rd1 - 2*g13*g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1*rd1 + 2*g12*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2*rd1 - 2*g12*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1*rd1 + g13*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2*rd1 - g13*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1*rd1 + 2*g23*g_load_1*k_2*k_ii_1*k_ii_2*k_iv_3*rd1 - 2*g23*g_load_1*k_2*k_ii_2*k_ii_3*k_iv_1*rd1 + 2*g23*g_load_1*k_3*k_ii_1*k_ii_3*k_iv_2*rd1 - 2*g23*g_load_1*k_3*k_ii_2*k_ii_3*k_iv_1*rd1 + g12*g_cpl_1*k_2*k_ii_1*k_ii_2*k_iv_3*rs1 - g12*g_cpl_1*k_2*k_ii_2*k_ii_3*k_iv_1*rs1 + g12*g_cpl_2*k_2*k_ii_1*k_ii_2*k_iv_3*rs1 - g12*g_cpl_2*k_2*k_ii_2*k_ii_3*k_iv_1*rs1 + g13*g_cpl_1*k_3*k_ii_1*k_ii_3*k_iv_2*rs1 - g13*g_cpl_1*k_3*k_ii_2*k_ii_3*k_iv_1*rs1 + 2*g13*g_cpl_2*k_2*k_ii_1*k_ii_2*k_iv_3*rs1 - 2*g13*g_cpl_2*k_2*k_ii_2*k_ii_3*k_iv_1*rs1 + 2*g12*g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2*rs1 - 2*g12*g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1*rs1 + g13*g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2*rs1 - g13*g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1*rs1 + 2*g23*g_cpl_1*k_2*k_ii_1*k_ii_2*k_iv_3*rs1 - 2*g23*g_cpl_1*k_2*k_ii_2*k_ii_3*k_iv_1*rs1 + 2*g23*g_cpl_1*k_3*k_ii_1*k_ii_3*k_iv_2*rs1 - 2*g23*g_cpl_1*k_3*k_ii_2*k_ii_3*k_iv_1*rs1 + g12*g_load_1*k_2*k_ii_1*k_ii_2*k_iv_3*rs1 - g12*g_load_1*k_2*k_ii_2*k_ii_3*k_iv_1*rs1 + g12*g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3*rs1 - g12*g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1*rs1 + g13*g_load_1*k_3*k_ii_1*k_ii_3*k_iv_2*rs1 - g13*g_load_1*k_3*k_ii_2*k_ii_3*k_iv_1*rs1 + 2*g13*g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3*rs1 - 2*g13*g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1*rs1 + 2*g12*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2*rs1 - 2*g12*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1*rs1 + g13*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2*rs1 - g13*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1*rs1 + 2*g23*g_load_1*k_2*k_ii_1*k_ii_2*k_iv_3*rs1 - 2*g23*g_load_1*k_2*k_ii_2*k_ii_3*k_iv_1*rs1 + 2*g23*g_load_1*k_3*k_ii_1*k_ii_3*k_iv_2*rs1 - 2*g23*g_load_1*k_3*k_ii_2*k_ii_3*k_iv_1*rs1 + g_cpl_1*g_cpl_2*k_2*k_ii_1*k_ii_2*k_iv_3*rd1 - g_cpl_1*g_cpl_2*k_2*k_ii_2*k_ii_3*k_iv_1*rd1 + g_cpl_1*g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2*rd1 - g_cpl_1*g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1*rd1 + g_cpl_1*g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3*rd1 - g_cpl_1*g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1*rd1 + g_cpl_2*g_load_1*k_2*k_ii_1*k_ii_2*k_iv_3*rd1 - g_cpl_2*g_load_1*k_2*k_ii_2*k_ii_3*k_iv_1*rd1 + g_cpl_1*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2*rd1 - g_cpl_1*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1*rd1 + g_cpl_3*g_load_1*k_3*k_ii_1*k_ii_3*k_iv_2*rd1 - g_cpl_3*g_load_1*k_3*k_ii_2*k_ii_3*k_iv_1*rd1 + g_load_1*g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3*rd1 - g_load_1*g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1*rd1 + g_load_1*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2*rd1 - g_load_1*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1*rd1 + g_cpl_1*g_cpl_2*k_2*k_ii_1*k_ii_2*k_iv_3*rs1 - g_cpl_1*g_cpl_2*k_2*k_ii_2*k_ii_3*k_iv_1*rs1 + g_cpl_1*g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2*rs1 - g_cpl_1*g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1*rs1 + g_cpl_1*g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3*rs1 - g_cpl_1*g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1*rs1 + g_cpl_2*g_load_1*k_2*k_ii_1*k_ii_2*k_iv_3*rs1 - g_cpl_2*g_load_1*k_2*k_ii_2*k_ii_3*k_iv_1*rs1 + g_cpl_1*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2*rs1 - g_cpl_1*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1*rs1 + g_cpl_3*g_load_1*k_3*k_ii_1*k_ii_3*k_iv_2*rs1 - g_cpl_3*g_load_1*k_3*k_ii_2*k_ii_3*k_iv_1*rs1 + g_load_1*g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3*rs1 - g_load_1*g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1*rs1 + g_load_1*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2*rs1 - g_load_1*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1*rs1 + g12*g13*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_3*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_load_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g12*g13*g_cpl_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g12*g13*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g12*g23*g_cpl_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g12*g23*g_cpl_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g13*g23*g_cpl_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g12*g23*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g13*g23*g_cpl_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g13*g23*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g12*g13*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g12*g13*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g12*g13*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g12*g23*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g12*g23*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g13*g23*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g12*g23*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g13*g23*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g13*g23*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g12*g13*g_cpl_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g12*g13*g_cpl_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g12*g13*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g12*g23*g_cpl_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g12*g23*g_cpl_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g13*g23*g_cpl_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g12*g23*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g13*g23*g_cpl_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g13*g23*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g12*g13*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g12*g13*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g12*g13*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g12*g23*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g12*g23*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g13*g23*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g12*g23*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g13*g23*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g13*g23*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g12*g_cpl_1*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g13*g_cpl_1*g_cpl_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g12*g_cpl_2*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g13*g_cpl_2*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g23*g_cpl_1*g_cpl_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g23*g_cpl_1*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g12*g_cpl_1*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g12*g_cpl_3*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g13*g_cpl_1*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g13*g_cpl_2*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g12*g_cpl_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g12*g_cpl_3*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g13*g_cpl_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g13*g_cpl_3*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g23*g_cpl_1*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g23*g_cpl_2*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g23*g_cpl_1*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g23*g_cpl_3*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g12*g_load_1*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g13*g_load_1*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g12*g_load_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g13*g_load_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g23*g_load_1*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g23*g_load_1*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g12*g_cpl_1*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g13*g_cpl_1*g_cpl_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g12*g_cpl_2*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g13*g_cpl_2*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g23*g_cpl_1*g_cpl_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g23*g_cpl_1*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g12*g_cpl_1*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g12*g_cpl_3*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g13*g_cpl_1*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g13*g_cpl_2*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g12*g_cpl_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g12*g_cpl_3*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g13*g_cpl_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g13*g_cpl_3*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g23*g_cpl_1*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g23*g_cpl_2*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g23*g_cpl_1*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g23*g_cpl_3*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g12*g_load_1*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g13*g_load_1*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g12*g_load_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g13*g_load_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g23*g_load_1*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g23*g_load_1*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g_cpl_1*g_cpl_2*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g_cpl_1*g_cpl_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g_cpl_1*g_cpl_3*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g_cpl_2*g_cpl_3*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g_cpl_1*g_load_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g_cpl_2*g_load_1*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g_cpl_3*g_load_1*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g_load_1*g_load_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd1 + g_cpl_1*g_cpl_2*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g_cpl_1*g_cpl_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g_cpl_1*g_cpl_3*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g_cpl_2*g_cpl_3*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g_cpl_1*g_load_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g_cpl_2*g_load_1*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g_cpl_3*g_load_1*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1 + g_load_1*g_load_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs1)/(3*g12*g13*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g13*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g12*g13*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g12*g23*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g23*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g13*g23*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g23*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g13*g23*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g13*g23*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g12*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g13*g_cpl_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g12*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g13*g_cpl_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g13*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g23*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g23*g_cpl_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_cpl_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g23*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g12*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g13*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g12*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g13*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g13*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g23*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g23*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g23*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_2*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_2*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_3*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_2*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_3*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_load_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_load_1*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_load_2*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g13*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_2*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_3*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_load_1*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n);
                d5= -(g12*k_1*k_ii_1*k_ii_2*k_iv_3 - g12*k_1*k_ii_1*k_ii_3*k_iv_2 + g12*k_2*k_ii_1*k_ii_2*k_iv_3 - g12*k_2*k_ii_2*k_ii_3*k_iv_1 + 2*g13*k_1*k_ii_1*k_ii_2*k_iv_3 - 2*g13*k_1*k_ii_1*k_ii_3*k_iv_2 - 2*g13*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*k_3*k_ii_2*k_ii_3*k_iv_1 - g23*k_2*k_ii_1*k_ii_2*k_iv_3 + g23*k_2*k_ii_2*k_ii_3*k_iv_1 - g23*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*k_1*k_ii_1*k_ii_2*k_iv_3 - g_cpl_1*k_1*k_ii_1*k_ii_3*k_iv_2 - g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1 + g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3 - g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2 - g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2 + g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g12*g13*k_1*k_ii_1*k_ii_2*k_iv_3*rd2 - 3*g12*g13*k_1*k_ii_1*k_ii_3*k_iv_2*rd2 - 3*g12*g13*k_3*k_ii_1*k_ii_3*k_iv_2*rd2 + 3*g12*g13*k_3*k_ii_2*k_ii_3*k_iv_1*rd2 + 3*g12*g23*k_1*k_ii_1*k_ii_2*k_iv_3*rd2 - 3*g12*g23*k_1*k_ii_1*k_ii_3*k_iv_2*rd2 + 3*g13*g23*k_1*k_ii_1*k_ii_2*k_iv_3*rd2 - 3*g13*g23*k_1*k_ii_1*k_ii_3*k_iv_2*rd2 - 3*g12*g23*k_3*k_ii_1*k_ii_3*k_iv_2*rd2 + 3*g12*g23*k_3*k_ii_2*k_ii_3*k_iv_1*rd2 - 3*g13*g23*k_3*k_ii_1*k_ii_3*k_iv_2*rd2 + 3*g13*g23*k_3*k_ii_2*k_ii_3*k_iv_1*rd2 + 3*g12*g13*k_1*k_ii_1*k_ii_2*k_iv_3*rs2 - 3*g12*g13*k_1*k_ii_1*k_ii_3*k_iv_2*rs2 - 3*g12*g13*k_3*k_ii_1*k_ii_3*k_iv_2*rs2 + 3*g12*g13*k_3*k_ii_2*k_ii_3*k_iv_1*rs2 + 3*g12*g23*k_1*k_ii_1*k_ii_2*k_iv_3*rs2 - 3*g12*g23*k_1*k_ii_1*k_ii_3*k_iv_2*rs2 + 3*g13*g23*k_1*k_ii_1*k_ii_2*k_iv_3*rs2 - 3*g13*g23*k_1*k_ii_1*k_ii_3*k_iv_2*rs2 - 3*g12*g23*k_3*k_ii_1*k_ii_3*k_iv_2*rs2 + 3*g12*g23*k_3*k_ii_2*k_ii_3*k_iv_1*rs2 - 3*g13*g23*k_3*k_ii_1*k_ii_3*k_iv_2*rs2 + 3*g13*g23*k_3*k_ii_2*k_ii_3*k_iv_1*rs2 + g12*g_cpl_1*k_1*k_ii_1*k_ii_2*k_iv_3*rd2 - g12*g_cpl_1*k_1*k_ii_1*k_ii_3*k_iv_2*rd2 + g12*g_cpl_2*k_1*k_ii_1*k_ii_2*k_iv_3*rd2 - g12*g_cpl_2*k_1*k_ii_1*k_ii_3*k_iv_2*rd2 + 2*g13*g_cpl_2*k_1*k_ii_1*k_ii_2*k_iv_3*rd2 - 2*g13*g_cpl_2*k_1*k_ii_1*k_ii_3*k_iv_2*rd2 - 2*g12*g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2*rd2 + 2*g12*g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1*rd2 - 2*g13*g_cpl_2*k_3*k_ii_1*k_ii_3*k_iv_2*rd2 + 2*g13*g_cpl_2*k_3*k_ii_2*k_ii_3*k_iv_1*rd2 + 2*g23*g_cpl_1*k_1*k_ii_1*k_ii_2*k_iv_3*rd2 - 2*g23*g_cpl_1*k_1*k_ii_1*k_ii_3*k_iv_2*rd2 - g23*g_cpl_2*k_3*k_ii_1*k_ii_3*k_iv_2*rd2 + g23*g_cpl_2*k_3*k_ii_2*k_ii_3*k_iv_1*rd2 - g23*g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2*rd2 + g23*g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1*rd2 + g12*g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3*rd2 - g12*g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2*rd2 + g12*g_load_2*k_1*k_ii_1*k_ii_2*k_iv_3*rd2 - g12*g_load_2*k_1*k_ii_1*k_ii_3*k_iv_2*rd2 + 2*g13*g_load_2*k_1*k_ii_1*k_ii_2*k_iv_3*rd2 - 2*g13*g_load_2*k_1*k_ii_1*k_ii_3*k_iv_2*rd2 - 2*g12*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2*rd2 + 2*g12*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1*rd2 - 2*g13*g_load_2*k_3*k_ii_1*k_ii_3*k_iv_2*rd2 + 2*g13*g_load_2*k_3*k_ii_2*k_ii_3*k_iv_1*rd2 + 2*g23*g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3*rd2 - 2*g23*g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2*rd2 - g23*g_load_2*k_3*k_ii_1*k_ii_3*k_iv_2*rd2 + g23*g_load_2*k_3*k_ii_2*k_ii_3*k_iv_1*rd2 - g23*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2*rd2 + g23*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1*rd2 + g12*g_cpl_1*k_1*k_ii_1*k_ii_2*k_iv_3*rs2 - g12*g_cpl_1*k_1*k_ii_1*k_ii_3*k_iv_2*rs2 + g12*g_cpl_2*k_1*k_ii_1*k_ii_2*k_iv_3*rs2 - g12*g_cpl_2*k_1*k_ii_1*k_ii_3*k_iv_2*rs2 + 2*g13*g_cpl_2*k_1*k_ii_1*k_ii_2*k_iv_3*rs2 - 2*g13*g_cpl_2*k_1*k_ii_1*k_ii_3*k_iv_2*rs2 - 2*g12*g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2*rs2 + 2*g12*g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1*rs2 - 2*g13*g_cpl_2*k_3*k_ii_1*k_ii_3*k_iv_2*rs2 + 2*g13*g_cpl_2*k_3*k_ii_2*k_ii_3*k_iv_1*rs2 + 2*g23*g_cpl_1*k_1*k_ii_1*k_ii_2*k_iv_3*rs2 - 2*g23*g_cpl_1*k_1*k_ii_1*k_ii_3*k_iv_2*rs2 - g23*g_cpl_2*k_3*k_ii_1*k_ii_3*k_iv_2*rs2 + g23*g_cpl_2*k_3*k_ii_2*k_ii_3*k_iv_1*rs2 - g23*g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2*rs2 + g23*g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1*rs2 + g12*g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3*rs2 - g12*g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2*rs2 + g12*g_load_2*k_1*k_ii_1*k_ii_2*k_iv_3*rs2 - g12*g_load_2*k_1*k_ii_1*k_ii_3*k_iv_2*rs2 + 2*g13*g_load_2*k_1*k_ii_1*k_ii_2*k_iv_3*rs2 - 2*g13*g_load_2*k_1*k_ii_1*k_ii_3*k_iv_2*rs2 - 2*g12*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2*rs2 + 2*g12*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1*rs2 - 2*g13*g_load_2*k_3*k_ii_1*k_ii_3*k_iv_2*rs2 + 2*g13*g_load_2*k_3*k_ii_2*k_ii_3*k_iv_1*rs2 + 2*g23*g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3*rs2 - 2*g23*g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2*rs2 - g23*g_load_2*k_3*k_ii_1*k_ii_3*k_iv_2*rs2 + g23*g_load_2*k_3*k_ii_2*k_ii_3*k_iv_1*rs2 - g23*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2*rs2 + g23*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1*rs2 + g_cpl_1*g_cpl_2*k_1*k_ii_1*k_ii_2*k_iv_3*rd2 - g_cpl_1*g_cpl_2*k_1*k_ii_1*k_ii_3*k_iv_2*rd2 - g_cpl_2*g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2*rd2 + g_cpl_2*g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1*rd2 + g_cpl_1*g_load_2*k_1*k_ii_1*k_ii_2*k_iv_3*rd2 - g_cpl_1*g_load_2*k_1*k_ii_1*k_ii_3*k_iv_2*rd2 + g_cpl_2*g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3*rd2 - g_cpl_2*g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2*rd2 - g_cpl_2*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2*rd2 + g_cpl_2*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1*rd2 - g_cpl_3*g_load_2*k_3*k_ii_1*k_ii_3*k_iv_2*rd2 + g_cpl_3*g_load_2*k_3*k_ii_2*k_ii_3*k_iv_1*rd2 + g_load_1*g_load_2*k_1*k_ii_1*k_ii_2*k_iv_3*rd2 - g_load_1*g_load_2*k_1*k_ii_1*k_ii_3*k_iv_2*rd2 - g_load_2*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2*rd2 + g_load_2*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1*rd2 + g_cpl_1*g_cpl_2*k_1*k_ii_1*k_ii_2*k_iv_3*rs2 - g_cpl_1*g_cpl_2*k_1*k_ii_1*k_ii_3*k_iv_2*rs2 - g_cpl_2*g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2*rs2 + g_cpl_2*g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1*rs2 + g_cpl_1*g_load_2*k_1*k_ii_1*k_ii_2*k_iv_3*rs2 - g_cpl_1*g_load_2*k_1*k_ii_1*k_ii_3*k_iv_2*rs2 + g_cpl_2*g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3*rs2 - g_cpl_2*g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2*rs2 - g_cpl_2*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2*rs2 + g_cpl_2*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1*rs2 - g_cpl_3*g_load_2*k_3*k_ii_1*k_ii_3*k_iv_2*rs2 + g_cpl_3*g_load_2*k_3*k_ii_2*k_ii_3*k_iv_1*rs2 + g_load_1*g_load_2*k_1*k_ii_1*k_ii_2*k_iv_3*rs2 - g_load_1*g_load_2*k_1*k_ii_1*k_ii_3*k_iv_2*rs2 - g_load_2*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2*rs2 + g_load_2*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1*rs2 + g12*g13*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_3*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_load_1*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g12*g13*g_cpl_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g12*g13*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g12*g23*g_cpl_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g12*g23*g_cpl_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g13*g23*g_cpl_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g12*g23*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g13*g23*g_cpl_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g13*g23*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g12*g13*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g12*g13*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g12*g13*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g12*g23*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g12*g23*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g13*g23*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g12*g23*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g13*g23*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g13*g23*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g12*g13*g_cpl_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g12*g13*g_cpl_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g12*g13*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g12*g23*g_cpl_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g12*g23*g_cpl_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g13*g23*g_cpl_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g12*g23*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g13*g23*g_cpl_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g13*g23*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g12*g13*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g12*g13*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g12*g13*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g12*g23*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g12*g23*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g13*g23*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g12*g23*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g13*g23*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g13*g23*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g12*g_cpl_1*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g13*g_cpl_1*g_cpl_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g12*g_cpl_2*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g13*g_cpl_2*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g23*g_cpl_1*g_cpl_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g23*g_cpl_1*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g12*g_cpl_1*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g12*g_cpl_3*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g13*g_cpl_1*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g13*g_cpl_2*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g12*g_cpl_2*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g12*g_cpl_3*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g13*g_cpl_2*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g13*g_cpl_3*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g23*g_cpl_1*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g23*g_cpl_2*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g23*g_cpl_1*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g23*g_cpl_3*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g12*g_load_1*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g13*g_load_1*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g12*g_load_2*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g13*g_load_2*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g23*g_load_1*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g23*g_load_1*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g12*g_cpl_1*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g13*g_cpl_1*g_cpl_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g12*g_cpl_2*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g13*g_cpl_2*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g23*g_cpl_1*g_cpl_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g23*g_cpl_1*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g12*g_cpl_1*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g12*g_cpl_3*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g13*g_cpl_1*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g13*g_cpl_2*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g12*g_cpl_2*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g12*g_cpl_3*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g13*g_cpl_2*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g13*g_cpl_3*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g23*g_cpl_1*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g23*g_cpl_2*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g23*g_cpl_1*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g23*g_cpl_3*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g12*g_load_1*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g13*g_load_1*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g12*g_load_2*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g13*g_load_2*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g23*g_load_1*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g23*g_load_1*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g_cpl_1*g_cpl_2*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g_cpl_1*g_cpl_2*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g_cpl_1*g_cpl_3*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g_cpl_2*g_cpl_3*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g_cpl_1*g_load_2*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g_cpl_2*g_load_1*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g_cpl_3*g_load_1*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g_load_1*g_load_2*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rd2 + g_cpl_1*g_cpl_2*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g_cpl_1*g_cpl_2*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g_cpl_1*g_cpl_3*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g_cpl_2*g_cpl_3*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g_cpl_1*g_load_2*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g_cpl_2*g_load_1*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g_cpl_3*g_load_1*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2 + g_load_1*g_load_2*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n*rs2)/(3*g12*g13*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g13*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g12*g13*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g12*g23*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g23*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g13*g23*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g23*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g13*g23*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g13*g23*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g12*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g13*g_cpl_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g12*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g13*g_cpl_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g13*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g23*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g23*g_cpl_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_cpl_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g23*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g12*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g13*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g12*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g13*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g13*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g23*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g23*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g23*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_2*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_2*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_3*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_2*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_3*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_load_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_load_1*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_load_2*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g13*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_2*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_3*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_load_1*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n);
                d6= -(2*g12*k_1*k_ii_1*k_ii_3*k_iv_2 - 2*g12*k_1*k_ii_1*k_ii_2*k_iv_3 - 2*g12*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*k_2*k_ii_2*k_ii_3*k_iv_1 - g13*k_1*k_ii_1*k_ii_2*k_iv_3 + g13*k_1*k_ii_1*k_ii_3*k_iv_2 + g13*k_3*k_ii_1*k_ii_3*k_iv_2 - g13*k_3*k_ii_2*k_ii_3*k_iv_1 - g23*k_2*k_ii_1*k_ii_2*k_iv_3 + g23*k_2*k_ii_2*k_ii_3*k_iv_1 - g23*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*k_3*k_ii_2*k_ii_3*k_iv_1 - g_cpl_1*k_1*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*k_1*k_ii_1*k_ii_3*k_iv_2 - g_cpl_2*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_2*k_2*k_ii_2*k_ii_3*k_iv_1 - g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3 + g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2 - g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3 + g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1 - 3*g12*g13*k_1*k_ii_1*k_ii_2*k_iv_3*rd3 + 3*g12*g13*k_1*k_ii_1*k_ii_3*k_iv_2*rd3 - 3*g12*g13*k_2*k_ii_1*k_ii_2*k_iv_3*rd3 + 3*g12*g13*k_2*k_ii_2*k_ii_3*k_iv_1*rd3 - 3*g12*g23*k_1*k_ii_1*k_ii_2*k_iv_3*rd3 + 3*g12*g23*k_1*k_ii_1*k_ii_3*k_iv_2*rd3 - 3*g12*g23*k_2*k_ii_1*k_ii_2*k_iv_3*rd3 + 3*g12*g23*k_2*k_ii_2*k_ii_3*k_iv_1*rd3 - 3*g13*g23*k_1*k_ii_1*k_ii_2*k_iv_3*rd3 + 3*g13*g23*k_1*k_ii_1*k_ii_3*k_iv_2*rd3 - 3*g13*g23*k_2*k_ii_1*k_ii_2*k_iv_3*rd3 + 3*g13*g23*k_2*k_ii_2*k_ii_3*k_iv_1*rd3 - 3*g12*g13*k_1*k_ii_1*k_ii_2*k_iv_3*rs3 + 3*g12*g13*k_1*k_ii_1*k_ii_3*k_iv_2*rs3 - 3*g12*g13*k_2*k_ii_1*k_ii_2*k_iv_3*rs3 + 3*g12*g13*k_2*k_ii_2*k_ii_3*k_iv_1*rs3 - 3*g12*g23*k_1*k_ii_1*k_ii_2*k_iv_3*rs3 + 3*g12*g23*k_1*k_ii_1*k_ii_3*k_iv_2*rs3 - 3*g12*g23*k_2*k_ii_1*k_ii_2*k_iv_3*rs3 + 3*g12*g23*k_2*k_ii_2*k_ii_3*k_iv_1*rs3 - 3*g13*g23*k_1*k_ii_1*k_ii_2*k_iv_3*rs3 + 3*g13*g23*k_1*k_ii_1*k_ii_3*k_iv_2*rs3 - 3*g13*g23*k_2*k_ii_1*k_ii_2*k_iv_3*rs3 + 3*g13*g23*k_2*k_ii_2*k_ii_3*k_iv_1*rs3 - g13*g_cpl_1*k_1*k_ii_1*k_ii_2*k_iv_3*rd3 + g13*g_cpl_1*k_1*k_ii_1*k_ii_3*k_iv_2*rd3 - 2*g12*g_cpl_3*k_1*k_ii_1*k_ii_2*k_iv_3*rd3 + 2*g12*g_cpl_3*k_1*k_ii_1*k_ii_3*k_iv_2*rd3 - 2*g12*g_cpl_3*k_2*k_ii_1*k_ii_2*k_iv_3*rd3 + 2*g12*g_cpl_3*k_2*k_ii_2*k_ii_3*k_iv_1*rd3 - 2*g13*g_cpl_2*k_2*k_ii_1*k_ii_2*k_iv_3*rd3 + 2*g13*g_cpl_2*k_2*k_ii_2*k_ii_3*k_iv_1*rd3 - g13*g_cpl_3*k_1*k_ii_1*k_ii_2*k_iv_3*rd3 + g13*g_cpl_3*k_1*k_ii_1*k_ii_3*k_iv_2*rd3 - 2*g23*g_cpl_1*k_1*k_ii_1*k_ii_2*k_iv_3*rd3 + 2*g23*g_cpl_1*k_1*k_ii_1*k_ii_3*k_iv_2*rd3 - g23*g_cpl_2*k_2*k_ii_1*k_ii_2*k_iv_3*rd3 + g23*g_cpl_2*k_2*k_ii_2*k_ii_3*k_iv_1*rd3 - g23*g_cpl_3*k_2*k_ii_1*k_ii_2*k_iv_3*rd3 + g23*g_cpl_3*k_2*k_ii_2*k_ii_3*k_iv_1*rd3 - g13*g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3*rd3 + g13*g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2*rd3 - 2*g12*g_load_3*k_1*k_ii_1*k_ii_2*k_iv_3*rd3 + 2*g12*g_load_3*k_1*k_ii_1*k_ii_3*k_iv_2*rd3 - 2*g12*g_load_3*k_2*k_ii_1*k_ii_2*k_iv_3*rd3 + 2*g12*g_load_3*k_2*k_ii_2*k_ii_3*k_iv_1*rd3 - 2*g13*g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3*rd3 + 2*g13*g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1*rd3 - g13*g_load_3*k_1*k_ii_1*k_ii_2*k_iv_3*rd3 + g13*g_load_3*k_1*k_ii_1*k_ii_3*k_iv_2*rd3 - 2*g23*g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3*rd3 + 2*g23*g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2*rd3 - g23*g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3*rd3 + g23*g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1*rd3 - g23*g_load_3*k_2*k_ii_1*k_ii_2*k_iv_3*rd3 + g23*g_load_3*k_2*k_ii_2*k_ii_3*k_iv_1*rd3 - g13*g_cpl_1*k_1*k_ii_1*k_ii_2*k_iv_3*rs3 + g13*g_cpl_1*k_1*k_ii_1*k_ii_3*k_iv_2*rs3 - 2*g12*g_cpl_3*k_1*k_ii_1*k_ii_2*k_iv_3*rs3 + 2*g12*g_cpl_3*k_1*k_ii_1*k_ii_3*k_iv_2*rs3 - 2*g12*g_cpl_3*k_2*k_ii_1*k_ii_2*k_iv_3*rs3 + 2*g12*g_cpl_3*k_2*k_ii_2*k_ii_3*k_iv_1*rs3 - 2*g13*g_cpl_2*k_2*k_ii_1*k_ii_2*k_iv_3*rs3 + 2*g13*g_cpl_2*k_2*k_ii_2*k_ii_3*k_iv_1*rs3 - g13*g_cpl_3*k_1*k_ii_1*k_ii_2*k_iv_3*rs3 + g13*g_cpl_3*k_1*k_ii_1*k_ii_3*k_iv_2*rs3 - 2*g23*g_cpl_1*k_1*k_ii_1*k_ii_2*k_iv_3*rs3 + 2*g23*g_cpl_1*k_1*k_ii_1*k_ii_3*k_iv_2*rs3 - g23*g_cpl_2*k_2*k_ii_1*k_ii_2*k_iv_3*rs3 + g23*g_cpl_2*k_2*k_ii_2*k_ii_3*k_iv_1*rs3 - g23*g_cpl_3*k_2*k_ii_1*k_ii_2*k_iv_3*rs3 + g23*g_cpl_3*k_2*k_ii_2*k_ii_3*k_iv_1*rs3 - g13*g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3*rs3 + g13*g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2*rs3 - 2*g12*g_load_3*k_1*k_ii_1*k_ii_2*k_iv_3*rs3 + 2*g12*g_load_3*k_1*k_ii_1*k_ii_3*k_iv_2*rs3 - 2*g12*g_load_3*k_2*k_ii_1*k_ii_2*k_iv_3*rs3 + 2*g12*g_load_3*k_2*k_ii_2*k_ii_3*k_iv_1*rs3 - 2*g13*g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3*rs3 + 2*g13*g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1*rs3 - g13*g_load_3*k_1*k_ii_1*k_ii_2*k_iv_3*rs3 + g13*g_load_3*k_1*k_ii_1*k_ii_3*k_iv_2*rs3 - 2*g23*g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3*rs3 + 2*g23*g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2*rs3 - g23*g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3*rs3 + g23*g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1*rs3 - g23*g_load_3*k_2*k_ii_1*k_ii_2*k_iv_3*rs3 + g23*g_load_3*k_2*k_ii_2*k_ii_3*k_iv_1*rs3 - g_cpl_1*g_cpl_3*k_1*k_ii_1*k_ii_2*k_iv_3*rd3 + g_cpl_1*g_cpl_3*k_1*k_ii_1*k_ii_3*k_iv_2*rd3 - g_cpl_2*g_cpl_3*k_2*k_ii_1*k_ii_2*k_iv_3*rd3 + g_cpl_2*g_cpl_3*k_2*k_ii_2*k_ii_3*k_iv_1*rd3 - g_cpl_1*g_load_3*k_1*k_ii_1*k_ii_2*k_iv_3*rd3 + g_cpl_1*g_load_3*k_1*k_ii_1*k_ii_3*k_iv_2*rd3 - g_cpl_3*g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3*rd3 + g_cpl_3*g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2*rd3 - g_cpl_2*g_load_3*k_2*k_ii_1*k_ii_2*k_iv_3*rd3 + g_cpl_2*g_load_3*k_2*k_ii_2*k_ii_3*k_iv_1*rd3 - g_cpl_3*g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3*rd3 + g_cpl_3*g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1*rd3 - g_load_1*g_load_3*k_1*k_ii_1*k_ii_2*k_iv_3*rd3 + g_load_1*g_load_3*k_1*k_ii_1*k_ii_3*k_iv_2*rd3 - g_load_2*g_load_3*k_2*k_ii_1*k_ii_2*k_iv_3*rd3 + g_load_2*g_load_3*k_2*k_ii_2*k_ii_3*k_iv_1*rd3 - g_cpl_1*g_cpl_3*k_1*k_ii_1*k_ii_2*k_iv_3*rs3 + g_cpl_1*g_cpl_3*k_1*k_ii_1*k_ii_3*k_iv_2*rs3 - g_cpl_2*g_cpl_3*k_2*k_ii_1*k_ii_2*k_iv_3*rs3 + g_cpl_2*g_cpl_3*k_2*k_ii_2*k_ii_3*k_iv_1*rs3 - g_cpl_1*g_load_3*k_1*k_ii_1*k_ii_2*k_iv_3*rs3 + g_cpl_1*g_load_3*k_1*k_ii_1*k_ii_3*k_iv_2*rs3 - g_cpl_3*g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3*rs3 + g_cpl_3*g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2*rs3 - g_cpl_2*g_load_3*k_2*k_ii_1*k_ii_2*k_iv_3*rs3 + g_cpl_2*g_load_3*k_2*k_ii_2*k_ii_3*k_iv_1*rs3 - g_cpl_3*g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3*rs3 + g_cpl_3*g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1*rs3 - g_load_1*g_load_3*k_1*k_ii_1*k_ii_2*k_iv_3*rs3 + g_load_1*g_load_3*k_1*k_ii_1*k_ii_3*k_iv_2*rs3 - g_load_2*g_load_3*k_2*k_ii_1*k_ii_2*k_iv_3*rs3 + g_load_2*g_load_3*k_2*k_ii_2*k_ii_3*k_iv_1*rs3 + g12*g13*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g_load_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g12*g13*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g12*g13*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g12*g23*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g12*g23*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g13*g23*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g12*g23*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g13*g23*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g13*g23*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g12*g13*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g12*g13*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g12*g13*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g12*g23*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g12*g23*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g13*g23*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g12*g23*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g13*g23*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g13*g23*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g12*g13*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g12*g13*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g12*g13*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g12*g23*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g12*g23*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g13*g23*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g12*g23*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g13*g23*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g13*g23*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g12*g13*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g12*g13*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g12*g13*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g12*g23*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g12*g23*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g13*g23*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g12*g23*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g13*g23*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g13*g23*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g12*g_cpl_1*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g13*g_cpl_1*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g12*g_cpl_2*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g13*g_cpl_2*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g23*g_cpl_1*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g23*g_cpl_1*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g12*g_cpl_1*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g12*g_cpl_3*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g13*g_cpl_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g13*g_cpl_2*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g12*g_cpl_2*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g12*g_cpl_3*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g13*g_cpl_2*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g13*g_cpl_3*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g23*g_cpl_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g23*g_cpl_2*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g23*g_cpl_1*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g23*g_cpl_3*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g12*g_load_1*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g13*g_load_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g12*g_load_2*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g13*g_load_2*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g23*g_load_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g23*g_load_1*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g12*g_cpl_1*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g13*g_cpl_1*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g12*g_cpl_2*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g13*g_cpl_2*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g23*g_cpl_1*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g23*g_cpl_1*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g12*g_cpl_1*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g12*g_cpl_3*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g13*g_cpl_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g13*g_cpl_2*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g12*g_cpl_2*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g12*g_cpl_3*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g13*g_cpl_2*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g13*g_cpl_3*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g23*g_cpl_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g23*g_cpl_2*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g23*g_cpl_1*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g23*g_cpl_3*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g12*g_load_1*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g13*g_load_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g12*g_load_2*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g13*g_load_2*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g23*g_load_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g23*g_load_1*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g_cpl_1*g_cpl_2*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g_cpl_1*g_cpl_2*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g_cpl_1*g_cpl_3*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g_cpl_2*g_cpl_3*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g_cpl_1*g_load_2*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g_cpl_2*g_load_1*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g_cpl_3*g_load_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g_load_1*g_load_2*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rd3 + g_cpl_1*g_cpl_2*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g_cpl_1*g_cpl_2*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g_cpl_1*g_cpl_3*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g_cpl_2*g_cpl_3*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g_cpl_1*g_load_2*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g_cpl_2*g_load_1*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g_cpl_3*g_load_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3 + g_load_1*g_load_2*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n*rs3)/(3*g12*g13*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g13*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g12*g13*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g12*g23*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g23*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g13*g23*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g23*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g13*g23*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g13*g23*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g12*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g13*g_cpl_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g12*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g13*g_cpl_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g13*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g23*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g23*g_cpl_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_cpl_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g23*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g12*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g13*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g12*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g13*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g13*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g23*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g23*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g23*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_2*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_2*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_3*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_2*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_3*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_load_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_load_1*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_load_2*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g13*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_2*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_3*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_load_1*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n);
                d7= -(3*g12*g13*k_2*k_ii_1*k_ii_2*k_iv_3 - 3*g12*g13*k_2*k_ii_2*k_ii_3*k_iv_1 + 3*g12*g13*k_3*k_ii_1*k_ii_3*k_iv_2 - 3*g12*g13*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g12*g23*k_2*k_ii_1*k_ii_2*k_iv_3 - 3*g12*g23*k_2*k_ii_2*k_ii_3*k_iv_1 + 3*g12*g23*k_3*k_ii_1*k_ii_3*k_iv_2 - 3*g12*g23*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g13*g23*k_2*k_ii_1*k_ii_2*k_iv_3 - 3*g13*g23*k_2*k_ii_2*k_ii_3*k_iv_1 + 3*g13*g23*k_3*k_ii_1*k_ii_3*k_iv_2 - 3*g13*g23*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_cpl_1*k_2*k_ii_1*k_ii_2*k_iv_3 - g12*g_cpl_1*k_2*k_ii_2*k_ii_3*k_iv_1 + g12*g_cpl_2*k_2*k_ii_1*k_ii_2*k_iv_3 - g12*g_cpl_2*k_2*k_ii_2*k_ii_3*k_iv_1 + g13*g_cpl_1*k_3*k_ii_1*k_ii_3*k_iv_2 - g13*g_cpl_1*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g13*g_cpl_2*k_2*k_ii_1*k_ii_2*k_iv_3 - 2*g13*g_cpl_2*k_2*k_ii_2*k_ii_3*k_iv_1 + 2*g12*g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2 - 2*g12*g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1 + g13*g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2 - g13*g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g23*g_cpl_1*k_2*k_ii_1*k_ii_2*k_iv_3 - 2*g23*g_cpl_1*k_2*k_ii_2*k_ii_3*k_iv_1 + 2*g23*g_cpl_1*k_3*k_ii_1*k_ii_3*k_iv_2 - 2*g23*g_cpl_1*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_load_1*k_2*k_ii_1*k_ii_2*k_iv_3 - g12*g_load_1*k_2*k_ii_2*k_ii_3*k_iv_1 + g12*g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3 - g12*g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1 + g13*g_load_1*k_3*k_ii_1*k_ii_3*k_iv_2 - g13*g_load_1*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g13*g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3 - 2*g13*g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1 + 2*g12*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2 - 2*g12*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1 + g13*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2 - g13*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g23*g_load_1*k_2*k_ii_1*k_ii_2*k_iv_3 - 2*g23*g_load_1*k_2*k_ii_2*k_ii_3*k_iv_1 + 2*g23*g_load_1*k_3*k_ii_1*k_ii_3*k_iv_2 - 2*g23*g_load_1*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_cpl_2*k_2*k_ii_1*k_ii_2*k_iv_3 - g_cpl_1*g_cpl_2*k_2*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2 - g_cpl_1*g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3 - g_cpl_1*g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1 + g_cpl_2*g_load_1*k_2*k_ii_1*k_ii_2*k_iv_3 - g_cpl_2*g_load_1*k_2*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2 - g_cpl_1*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_3*g_load_1*k_3*k_ii_1*k_ii_3*k_iv_2 - g_cpl_3*g_load_1*k_3*k_ii_2*k_ii_3*k_iv_1 + g_load_1*g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3 - g_load_1*g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1 + g_load_1*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2 - g_load_1*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g13*g_cpl_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_cpl_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_3*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_2*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_3*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_1*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_1*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_cpl_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_3*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_cpl_3*g_load_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_load_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_load_1*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_3*g_load_1*g_load_2*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_load_1*g_load_2*g_load_3*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n)/(3*g12*g13*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g13*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g12*g13*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g12*g23*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g23*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g13*g23*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g23*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g13*g23*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g13*g23*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g12*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g13*g_cpl_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g12*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g13*g_cpl_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g13*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g23*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g23*g_cpl_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_cpl_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g23*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g12*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g13*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g12*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g13*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g13*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g23*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g23*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g23*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_2*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_2*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_3*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_2*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_3*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_load_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_load_1*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_load_2*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g13*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_2*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_3*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_load_1*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n);
                d8= -(3*g12*g13*k_1*k_ii_1*k_ii_2*k_iv_3 - 3*g12*g13*k_1*k_ii_1*k_ii_3*k_iv_2 - 3*g12*g13*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g12*g13*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g12*g23*k_1*k_ii_1*k_ii_2*k_iv_3 - 3*g12*g23*k_1*k_ii_1*k_ii_3*k_iv_2 + 3*g13*g23*k_1*k_ii_1*k_ii_2*k_iv_3 - 3*g13*g23*k_1*k_ii_1*k_ii_3*k_iv_2 - 3*g12*g23*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g12*g23*k_3*k_ii_2*k_ii_3*k_iv_1 - 3*g13*g23*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g13*g23*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_cpl_1*k_1*k_ii_1*k_ii_2*k_iv_3 - g12*g_cpl_1*k_1*k_ii_1*k_ii_3*k_iv_2 + g12*g_cpl_2*k_1*k_ii_1*k_ii_2*k_iv_3 - g12*g_cpl_2*k_1*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_cpl_2*k_1*k_ii_1*k_ii_2*k_iv_3 - 2*g13*g_cpl_2*k_1*k_ii_1*k_ii_3*k_iv_2 - 2*g12*g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g12*g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1 - 2*g13*g_cpl_2*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_cpl_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g23*g_cpl_1*k_1*k_ii_1*k_ii_2*k_iv_3 - 2*g23*g_cpl_1*k_1*k_ii_1*k_ii_3*k_iv_2 - g23*g_cpl_2*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_cpl_2*k_3*k_ii_2*k_ii_3*k_iv_1 - g23*g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3 - g12*g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2 + g12*g_load_2*k_1*k_ii_1*k_ii_2*k_iv_3 - g12*g_load_2*k_1*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_load_2*k_1*k_ii_1*k_ii_2*k_iv_3 - 2*g13*g_load_2*k_1*k_ii_1*k_ii_3*k_iv_2 - 2*g12*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g12*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1 - 2*g13*g_load_2*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_load_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g23*g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3 - 2*g23*g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2 - g23*g_load_2*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_load_2*k_3*k_ii_2*k_ii_3*k_iv_1 - g23*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_cpl_2*k_1*k_ii_1*k_ii_2*k_iv_3 - g_cpl_1*g_cpl_2*k_1*k_ii_1*k_ii_3*k_iv_2 - g_cpl_2*g_cpl_3*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_2*g_cpl_3*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_load_2*k_1*k_ii_1*k_ii_2*k_iv_3 - g_cpl_1*g_load_2*k_1*k_ii_1*k_ii_3*k_iv_2 + g_cpl_2*g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3 - g_cpl_2*g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2 - g_cpl_2*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_2*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1 - g_cpl_3*g_load_2*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_3*g_load_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_load_1*g_load_2*k_1*k_ii_1*k_ii_2*k_iv_3 - g_load_1*g_load_2*k_1*k_ii_1*k_ii_3*k_iv_2 - g_load_2*g_load_3*k_3*k_ii_1*k_ii_3*k_iv_2 + g_load_2*g_load_3*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g13*g_cpl_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_cpl_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_3*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_2*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_3*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_1*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_1*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_2*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_2*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_cpl_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_3*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_cpl_3*g_load_1*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_load_2*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_load_1*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_3*g_load_1*g_load_2*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_load_1*g_load_2*g_load_3*k_1*k_3*k_ii_1*k_ii_2*k_ii_3*n)/(3*g12*g13*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g13*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g12*g13*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g12*g23*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g23*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g13*g23*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g23*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g13*g23*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g13*g23*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g12*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g13*g_cpl_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g12*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g13*g_cpl_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g13*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g23*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g23*g_cpl_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_cpl_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g23*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g12*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g13*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g12*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g13*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g13*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g23*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g23*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g23*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_2*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_2*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_3*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_2*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_3*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_load_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_load_1*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_load_2*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g13*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_2*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_3*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_load_1*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n);
                d9= -(3*g12*g13*k_1*k_ii_1*k_ii_3*k_iv_2 - 3*g12*g13*k_1*k_ii_1*k_ii_2*k_iv_3 - 3*g12*g13*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g13*k_2*k_ii_2*k_ii_3*k_iv_1 - 3*g12*g23*k_1*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g23*k_1*k_ii_1*k_ii_3*k_iv_2 - 3*g12*g23*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g23*k_2*k_ii_2*k_ii_3*k_iv_1 - 3*g13*g23*k_1*k_ii_1*k_ii_2*k_iv_3 + 3*g13*g23*k_1*k_ii_1*k_ii_3*k_iv_2 - 3*g13*g23*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g13*g23*k_2*k_ii_2*k_ii_3*k_iv_1 - g13*g_cpl_1*k_1*k_ii_1*k_ii_2*k_iv_3 + g13*g_cpl_1*k_1*k_ii_1*k_ii_3*k_iv_2 - 2*g12*g_cpl_3*k_1*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_cpl_3*k_1*k_ii_1*k_ii_3*k_iv_2 - 2*g12*g_cpl_3*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_cpl_3*k_2*k_ii_2*k_ii_3*k_iv_1 - 2*g13*g_cpl_2*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g13*g_cpl_2*k_2*k_ii_2*k_ii_3*k_iv_1 - g13*g_cpl_3*k_1*k_ii_1*k_ii_2*k_iv_3 + g13*g_cpl_3*k_1*k_ii_1*k_ii_3*k_iv_2 - 2*g23*g_cpl_1*k_1*k_ii_1*k_ii_2*k_iv_3 + 2*g23*g_cpl_1*k_1*k_ii_1*k_ii_3*k_iv_2 - g23*g_cpl_2*k_2*k_ii_1*k_ii_2*k_iv_3 + g23*g_cpl_2*k_2*k_ii_2*k_ii_3*k_iv_1 - g23*g_cpl_3*k_2*k_ii_1*k_ii_2*k_iv_3 + g23*g_cpl_3*k_2*k_ii_2*k_ii_3*k_iv_1 - g13*g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3 + g13*g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2 - 2*g12*g_load_3*k_1*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_load_3*k_1*k_ii_1*k_ii_3*k_iv_2 - 2*g12*g_load_3*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_load_3*k_2*k_ii_2*k_ii_3*k_iv_1 - 2*g13*g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g13*g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1 - g13*g_load_3*k_1*k_ii_1*k_ii_2*k_iv_3 + g13*g_load_3*k_1*k_ii_1*k_ii_3*k_iv_2 - 2*g23*g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3 + 2*g23*g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2 - g23*g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3 + g23*g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1 - g23*g_load_3*k_2*k_ii_1*k_ii_2*k_iv_3 + g23*g_load_3*k_2*k_ii_2*k_ii_3*k_iv_1 - g_cpl_1*g_cpl_3*k_1*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*g_cpl_3*k_1*k_ii_1*k_ii_3*k_iv_2 - g_cpl_2*g_cpl_3*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_2*g_cpl_3*k_2*k_ii_2*k_ii_3*k_iv_1 - g_cpl_1*g_load_3*k_1*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*g_load_3*k_1*k_ii_1*k_ii_3*k_iv_2 - g_cpl_3*g_load_1*k_1*k_ii_1*k_ii_2*k_iv_3 + g_cpl_3*g_load_1*k_1*k_ii_1*k_ii_3*k_iv_2 - g_cpl_2*g_load_3*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_2*g_load_3*k_2*k_ii_2*k_ii_3*k_iv_1 - g_cpl_3*g_load_2*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_3*g_load_2*k_2*k_ii_2*k_ii_3*k_iv_1 - g_load_1*g_load_3*k_1*k_ii_1*k_ii_2*k_iv_3 + g_load_1*g_load_3*k_1*k_ii_1*k_ii_3*k_iv_2 - g_load_2*g_load_3*k_2*k_ii_1*k_ii_2*k_iv_3 + g_load_2*g_load_3*k_2*k_ii_2*k_ii_3*k_iv_1 + g12*g13*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_3*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_2*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_3*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_1*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_2*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_2*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_cpl_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_3*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_cpl_3*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_load_2*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_load_1*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_3*g_load_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n + g_load_1*g_load_2*g_load_3*k_1*k_2*k_ii_1*k_ii_2*k_ii_3*n)/(3*g12*g13*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g13*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g12*g13*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g12*g23*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g23*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g13*g23*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 3*g12*g23*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 3*g13*g23*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 3*g13*g23*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g12*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g13*g_cpl_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g12*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g13*g_cpl_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g13*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g23*g_cpl_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g23*g_cpl_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_cpl_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g23*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g12*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g13*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g13*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g12*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g12*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + 2*g13*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g13*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + 2*g23*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + 2*g23*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g23*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g23*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_cpl_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*g_cpl_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_2*g_cpl_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_2*g_load_1*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_cpl_1*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_3*g_load_1*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_cpl_2*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_cpl_3*g_load_2*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g_load_1*g_load_2*k_1*k_2*k_ii_1*k_ii_2*k_iv_3 + g_load_1*g_load_3*k_1*k_3*k_ii_1*k_ii_3*k_iv_2 + g_load_2*g_load_3*k_2*k_3*k_ii_2*k_ii_3*k_iv_1 + g12*g13*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g13*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g23*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g23*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_2*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g12*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g13*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g23*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_cpl_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_cpl_3*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_cpl_3*g_load_1*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_1*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_2*g_load_1*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_cpl_3*g_load_1*g_load_2*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n + g_load_1*g_load_2*g_load_3*k_1*k_2*k_3*k_ii_1*k_ii_2*k_ii_3*n);
                d10= -1;
                d11= -1;
                d12= 0;

                aux = [d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12]';
                D(:,i) = aux;
            end
        
        % Passive Filter
            r_sk    = 0.5; 
            R_SK    = r_sk * [0.9 1.05 1];% I consider this value more realistic than 0.05
            R_SK    = R_SK(1:n);

            L_sk    = 1e-3;
            IND     = L_sk * [1 0.85 1.15];IND    = IND(1:n);

            C_bk    = 1e-3; 
            CAP     = C_bk * [0.33 0.5 1];CAP    = CAP(1:n);

        %Interconnection   
            topology =  [0 1 1; ...
                         1 0 1; ...
                         1 1 0];   topology = topology(1:n,1:n);      
            r_jk    = 2; 
            R_JK    = r_jk * ones(1,n);
            RJK     = topology * r_jk;
            for i=1:n
                for j=1:n
                    if RJK(i,j) ~= 0
                    RJK(i,j)     = -1*RJK(i,j)^(-1);
                    end
                end
            end
            RJK     = RJK + diag(-1*sum(RJK)); 


            AG_i = [  0 1 0 ; ...
                      1 0 1 ; ...
                      0 1 0 ]; %Communication matrices
            AG_i = AG_i(1:n,1:n);

        %Load
            R_passive   = [1 1 1]; R_passive    = R_passive(1:n);
            R_load      = 6e8 * R_passive;

        %Adaptive droop
            %ki = 1;
            kiv_i = k_iv_1;  k_iv = ones(1,n) * kiv_i;

            kii_i = k_ii_1;     k_ii = ones(1,n) * kii_i;
            kpi_i = k_pi_1; k_pi = ones(1,n) * kpi_i; 

        %Current sharing
            r_d0 = 1.3; 
            R_D0 = r_d0* ones(1,n); 
            RD0 = r_d0* diag(ones(1,n));

            k_i = [1, 2, 1]; k_i = k_i(1:n);
            K = diag(k_i.^(-1)); 
            RDI = RD0 * K;
            I_rate =[10, 20, 10]; I_rate =I_rate(1:n);
            I_RATE = diag(I_rate.^(-1));

        %MATRIX GENERATION    
            A_glob = zeros(vertex,4*n-1,4*n-1); 
            B_glob = zeros(vertex,4*n-1,4*n-1);  
            eigenvalues = [];

            for i=1:vertex

                 R_CPL = -(V_ref^2)./p_load_vertex(i,:) ;

                 m1_L    = diag(IND.^(-1));
                 m1_C      = diag(CAP.^(-1));
                 m1_C2     = m1_C * ones(n,n);

                 mK_II   = diag(k_ii);
                 mK_IV   = diag(k_iv);
                 mK_PI   = diag(k_pi);


                 G_CPL     = R_CPL.^(-1); G_CPL = diag(G_CPL);
                 G_ast     = [g_load_1, g_load_2 , g_load_3]; G_ast = diag(G_ast);
                 mG        = G_CPL + G_ast;


                 mLE       = RJK;
                 mLA = diag(sum(AG_i)) - AG_i;

                 cero = zeros(n,n);

                 A = [-1*RDI*m1_L   , ((RDI*m1_L) - (1/n.*mK_IV))   , ((RDI*diag(R_SK)*m1_L)-(mK_II*I_RATE))    , (mK_II - mK_PI*diag(diag(mLA)))   ;...
                     cero           , -1*(mLE + mG).*m1_C2          , m1_C                                      ,   cero                            ;...
                     m1_L           , -1.*m1_L                      , -1.*(diag(R_SK)*m1_L)                     , cero                              ;...
                     I_RATE*m1_L    , -1.*( I_RATE*m1_L)            , -1.*(diag(R_SK)*I_RATE*m1_L)              , -1.*(diag(diag(mLA)))             ];

                 aux = ones(n,n) - eye(n);
                 mK_IV2  = mK_IV * ones(n,n);

                 B = [cero          , (-1/n) .* (aux .*mK_IV2)      , cero                                      ,   -1.*mK_PI*(mLA-diag(diag(mLA))) ;...
                      cero          , cero                          , cero                                      ,   cero                            ;...
                      cero          , cero                          , cero                                      ,   cero                            ;...
                      cero          , cero                          , cero                                      , -1.*(mLA-diag(diag(mLA)))         ];   


                  %con u = [VREF, VREFc,I_load_tot]
                 H = [cero          , mK_IV                         , cero      ;...
                     cero           , cero                          , -1.*(m1_C);...
                     cero           , cero                          , cero      ;...
                     cero           , cero                          , cero       ]; 



                 B_sim = [B, H];
                 %C = diag(ones(1,4*n));
                 %D = zeros(size(B_sim));

                 A_correction = D(:,i) * A(end,:) ;
                 B_correction = D(:,i) * B(end,:) ;
                 A_prim = A + A_correction;
                 B_prim = B + B_correction;
                 A_prim = A_prim(1:(end-1),1:(end-1));
                 B_prim = B_prim(1:(end-1),1:(end-1));

                 eigenvalues = [eigenvalues , eig(A_prim+B_prim)];

                 A_glob(i,:,:) = A_prim;
                 B_glob(i,:,:) = B_prim;

            end % for vertex
            
            
        %-----------------------------------------------------------------
        % LMI
        %-----------------------------------------------------------------
            [n_lmi,m_lmi] = size(A_prim);
                                                                                            
              cvx_clear
              %cvx_solver mosek

              cero_lmi = zeros(n_lmi,n_lmi);
              cero_lmi_vector = [cero_lmi, cero_lmi, cero_lmi,cero_lmi];

              cvx_begin sdp
                % P
                variable P1(n_lmi,n_lmi) symmetric
                variable P2(n_lmi,n_lmi) symmetric
                variable P3(n_lmi,n_lmi) symmetric
                variable P4(n_lmi,n_lmi) symmetric
                variable P5(n_lmi,n_lmi) symmetric
                variable P6(n_lmi,n_lmi) symmetric
                variable P7(n_lmi,n_lmi) symmetric
                variable P8(n_lmi,n_lmi) symmetric

                % Q 
                variable Q1(n_lmi,n_lmi) symmetric
                variable Q2(n_lmi,n_lmi) symmetric
                variable Q3(n_lmi,n_lmi) symmetric
                variable Q4(n_lmi,n_lmi) symmetric
                variable Q5(n_lmi,n_lmi) symmetric
                variable Q6(n_lmi,n_lmi) symmetric
                variable Q7(n_lmi,n_lmi) symmetric
                variable Q8(n_lmi,n_lmi) symmetric
                % R 
                variable R1(n_lmi,n_lmi) symmetric
                variable R2(n_lmi,n_lmi) symmetric
                variable R3(n_lmi,n_lmi) symmetric
                variable R4(n_lmi,n_lmi) symmetric
                variable R5(n_lmi,n_lmi) symmetric
                variable R6(n_lmi,n_lmi) symmetric
                variable R7(n_lmi,n_lmi) symmetric
                variable R8(n_lmi,n_lmi) symmetric
                % Z
                variable Z1_1(n_lmi,n_lmi) symmetric
                variable Z1_2(n_lmi,n_lmi) symmetric
                variable Z1_3(n_lmi,n_lmi) symmetric
                variable Z1_4(n_lmi,n_lmi) symmetric
                variable Z1_5(n_lmi,n_lmi) symmetric
                variable Z1_6(n_lmi,n_lmi) symmetric
                variable Z1_7(n_lmi,n_lmi) symmetric
                variable Z1_8(n_lmi,n_lmi) symmetric

                variable Z2_1(n_lmi,n_lmi) symmetric
                variable Z2_2(n_lmi,n_lmi) symmetric
                variable Z2_3(n_lmi,n_lmi) symmetric
                variable Z2_4(n_lmi,n_lmi) symmetric
                variable Z2_5(n_lmi,n_lmi) symmetric
                variable Z2_6(n_lmi,n_lmi) symmetric
                variable Z2_7(n_lmi,n_lmi) symmetric
                variable Z2_8(n_lmi,n_lmi) symmetric

                % N
                %variable N1_1(n_lmi,n_lmi)
                variable N1_1(n_lmi,n_lmi) 
                variable N1_2(n_lmi,n_lmi) 
                variable N1_3(n_lmi,n_lmi) 
                variable N1_4(n_lmi,n_lmi) 
                variable N1_5(n_lmi,n_lmi) 
                variable N1_6(n_lmi,n_lmi) 
                variable N1_7(n_lmi,n_lmi) 
                variable N1_8(n_lmi,n_lmi) 


                %variable N2(n_lmi,n_lmi)
                %variable N2_1(n_lmi,n_lmi)
                variable N2_1(n_lmi,n_lmi) 
                variable N2_2(n_lmi,n_lmi) 
                variable N2_3(n_lmi,n_lmi) 
                variable N2_4(n_lmi,n_lmi) 
                variable N2_5(n_lmi,n_lmi) 
                variable N2_6(n_lmi,n_lmi) 
                variable N2_7(n_lmi,n_lmi) 
                variable N2_8(n_lmi,n_lmi) 

                %variable N3(n_lmi,n_lmi)
                variable N3_1(n_lmi,n_lmi) 
                variable N3_2(n_lmi,n_lmi) 
                variable N3_3(n_lmi,n_lmi) 
                variable N3_4(n_lmi,n_lmi) 
                variable N3_5(n_lmi,n_lmi) 
                variable N3_6(n_lmi,n_lmi) 
                variable N3_7(n_lmi,n_lmi) 
                variable N3_8(n_lmi,n_lmi) 

                %variable N4(n_lmi,n_lmi)
                variable N4_1(n_lmi,n_lmi) 
                variable N4_2(n_lmi,n_lmi) 
                variable N4_3(n_lmi,n_lmi) 
                variable N4_4(n_lmi,n_lmi) 
                variable N4_5(n_lmi,n_lmi) 
                variable N4_6(n_lmi,n_lmi) 
                variable N4_7(n_lmi,n_lmi) 
                variable N4_8(n_lmi,n_lmi) 

                % S   
                %variable S1(n_lmi,n_lmi) 
                variable S1_1(n_lmi,n_lmi) 
                variable S1_2(n_lmi,n_lmi) 
                variable S1_3(n_lmi,n_lmi) 
                variable S1_4(n_lmi,n_lmi) 
                variable S1_5(n_lmi,n_lmi) 
                variable S1_6(n_lmi,n_lmi) 
                variable S1_7(n_lmi,n_lmi) 
                variable S1_8(n_lmi,n_lmi) 


                %variable S2(n_lmi,n_lmi)
                variable S2_1(n_lmi,n_lmi) 
                variable S2_2(n_lmi,n_lmi) 
                variable S2_3(n_lmi,n_lmi) 
                variable S2_4(n_lmi,n_lmi) 
                variable S2_5(n_lmi,n_lmi) 
                variable S2_6(n_lmi,n_lmi) 
                variable S2_7(n_lmi,n_lmi) 
                variable S2_8(n_lmi,n_lmi) 

                %variable S3(n_lmi,n_lmi)
                variable S3_1(n_lmi,n_lmi) 
                variable S3_2(n_lmi,n_lmi) 
                variable S3_3(n_lmi,n_lmi) 
                variable S3_4(n_lmi,n_lmi) 
                variable S3_5(n_lmi,n_lmi) 
                variable S3_6(n_lmi,n_lmi) 
                variable S3_7(n_lmi,n_lmi) 
                variable S3_8(n_lmi,n_lmi) 

                %variable S4(n_lmi,n_lmi)
                variable S4_1(n_lmi,n_lmi) 
                variable S4_2(n_lmi,n_lmi) 
                variable S4_3(n_lmi,n_lmi) 
                variable S4_4(n_lmi,n_lmi) 
                variable S4_5(n_lmi,n_lmi) 
                variable S4_6(n_lmi,n_lmi) 
                variable S4_7(n_lmi,n_lmi) 
                variable S4_8(n_lmi,n_lmi) 


                % M

                %variable M1(n_lmi,n_lmi) 
            %     variable M1_1(n_lmi,n_lmi) 
            %     variable M1_2(n_lmi,n_lmi) 
            %     variable M1_3(n_lmi,n_lmi) 
            %     variable M1_4(n_lmi,n_lmi) 
            %     variable M1_5(n_lmi,n_lmi) 
            %     variable M1_6(n_lmi,n_lmi) 
            %     variable M1_7(n_lmi,n_lmi) 
            %     variable M1_8(n_lmi,n_lmi) 
            %     
            %     
            %     %variable M2(n_lmi,n_lmi) 
            %     variable M2_1(n_lmi,n_lmi) 
            %     variable M2_2(n_lmi,n_lmi) 
            %     variable M2_3(n_lmi,n_lmi) 
            %     variable M2_4(n_lmi,n_lmi) 
            %     variable M2_5(n_lmi,n_lmi) 
            %     variable M2_6(n_lmi,n_lmi) 
            %     variable M2_7(n_lmi,n_lmi) 
            %     variable M2_8(n_lmi,n_lmi) 
            %     
            %     %variable M3(n_lmi,n_lmi)
            %     variable M3_1(n_lmi,n_lmi) 
            %     variable M3_2(n_lmi,n_lmi) 
            %     variable M3_3(n_lmi,n_lmi) 
            %     variable M3_4(n_lmi,n_lmi) 
            %     variable M3_5(n_lmi,n_lmi) 
            %     variable M3_6(n_lmi,n_lmi) 
            %     variable M3_7(n_lmi,n_lmi) 
            %     variable M3_8(n_lmi,n_lmi) 
            %     
            %     % variable M4(n_lmi,n_lmi) 
            %     variable M4_1(n_lmi,n_lmi) 
            %     variable M4_2(n_lmi,n_lmi) 
            %     variable M4_3(n_lmi,n_lmi) 
            %     variable M4_4(n_lmi,n_lmi) 
            %     variable M4_5(n_lmi,n_lmi) 
            %     variable M4_6(n_lmi,n_lmi) 
            %     variable M4_7(n_lmi,n_lmi) 
            %     variable M4_8(n_lmi,n_lmi) 

                % T
                variable T1(n_lmi,n_lmi) 
                variable T2(n_lmi,n_lmi) 
                variable T3(n_lmi,n_lmi)
                variable T4(n_lmi,n_lmi)


                %cvx_precision high

                minimize 1

                subject to:
                    T = [T1 T2 T3 T4]';
            %1
                        N_1 = [N1_1 N2_1 N3_1 N4_1]';

                        aux_1 = reshape(A_glob(1,:,:),[n_lmi,n_lmi]);
                        bux_1 = reshape(B_glob(1,:,:),[n_lmi,n_lmi]);
                        A_lmi_1 = [-aux_1, -bux_1, cero_lmi, eye(n_lmi)];

                        O1_11_1 = Q1 + R1;
                        O1_14_1 = P1;
                        O1_22_1 = -1*(1-mu)*Q1;
                        O1_33_1 = -1*R1;
                        O1_44_1 = h*(Z1_1+Z2_1);

                        O1_1 = [O1_11_1     , cero_lmi  , cero_lmi  , O1_14_1       ;...
                              cero_lmi      , O1_22_1   , cero_lmi  , cero_lmi      ;...
                              cero_lmi      , cero_lmi  , O1_33_1   , cero_lmi      ;...
                              O1_14_1'      , cero_lmi  , cero_lmi  , O1_44_1       ];

                        O2_1 = [N_1, -N_1, cero_lmi_vector', cero_lmi_vector'] + (T*A_lmi_1) + ((A_lmi_1')*(T'));

                        O_1 = O1_1 + O2_1 + O2_1';

                        %Definitions Condition 1

                        Condition_1 =    [O_1       , h*N_1         ,cero_lmi_vector'        , cero_lmi_vector'         ;...
                                         h*N_1'     , -h*Z1_1       , cero_lmi      , cero_lmi      ;...    
                                         cero_lmi_vector    , cero_lmi      , -h*Z1_1       , cero_lmi      ;...
                                         cero_lmi_vector     , cero_lmi      , cero_lmi      , -h*Z2_1       ];

                        Condition_1 <= -eye(7*n_lmi);
                        %Condition1 <= 0;
                        P1      >= eye(n_lmi);
                        Q1      >= 0;
                        R1      >= 0;
                        Z1_1    >= eye(n_lmi);
                        Z2_1    >= eye(n_lmi);

            %2
                        N_2 = [N1_2 N2_2 N3_2 N4_2]';



                        aux_2 = reshape(A_glob(2,:,:),[n_lmi,n_lmi]);
                        bux_2 = reshape(B_glob(2,:,:),[n_lmi,n_lmi]);
                        A_lmi_2 = [-aux_2, -bux_2, cero_lmi, eye(n_lmi)];

                        O1_11_2 = Q2 +R2;
                        O1_14_2 = P2;
                        O1_22_2 = -1*(1-mu)*Q2;
                        O1_33_2 = -1*R2;
                        O1_44_2 = h*(Z1_2+Z2_2);

                        O1_2 = [O1_11_2     , cero_lmi  , cero_lmi  , O1_14_2       ;...
                              cero_lmi      , O1_22_2   , cero_lmi  , cero_lmi      ;...
                              cero_lmi      , cero_lmi  , O1_33_2   , cero_lmi      ;...
                              O1_14_2'      , cero_lmi  , cero_lmi  , O1_44_2       ];

                        O2_2 = [N_2, -N_2, cero_lmi_vector', cero_lmi_vector'] + (T*A_lmi_2) + ((A_lmi_2')*(T'));

                        O_2 = O1_2 + O2_2 + O2_2';

                        %Definitions Condition 1

                        Condition_2 =    [O_2       , h*N_2         , cero_lmi_vector'         , cero_lmi_vector'         ;...
                                         h*N_2'     , -h*Z1_2       , cero_lmi      , cero_lmi      ;...    
                                         cero_lmi_vector   , cero_lmi      , -h*Z1_2       , cero_lmi      ;...
                                         cero_lmi_vector     , cero_lmi      , cero_lmi      , -h*Z2_2       ];

                        Condition_2 <= -eye(7*n_lmi);
                        %Condition1 <= 0;
                        P2      >= eye(n_lmi);
                        Q2      >= 0;
                        R2      >= 0;
                        Z1_2    >= eye(n_lmi);
                        Z2_2    >= eye(n_lmi);
            %3
                        N_3 = [N1_3 N2_3 N3_3 N4_3]';



                        aux_3 = reshape(A_glob(3,:,:),[n_lmi,n_lmi]);
                        bux_3 = reshape(B_glob(3,:,:),[n_lmi,n_lmi]);
                        A_lmi_3 = [-aux_3, -bux_3, cero_lmi, eye(n_lmi)];

                        O1_11_3 = Q3 + R3;
                        O1_14_3 = P3;
                        O1_22_3 = -1*(1-mu)*Q3;
                        O1_33_3 = -1*R3;
                        O1_44_3 = h*(Z1_3+Z2_3);

                        O1_3 = [O1_11_3     , cero_lmi  , cero_lmi  , O1_14_3       ;...
                              cero_lmi      , O1_22_3   , cero_lmi  , cero_lmi      ;...
                              cero_lmi      , cero_lmi  , O1_33_3   , cero_lmi      ;...
                              O1_14_3'      , cero_lmi  , cero_lmi  , O1_44_3       ];

                        O2_3 = [N_3, -N_3, cero_lmi_vector', cero_lmi_vector'] + (T*A_lmi_3) + ((A_lmi_3')*(T'));

                        O_3 = O1_3 + O2_3 + O2_3';

                        %Definitions Condition 1

                        Condition_3 =    [O_3       , h*N_3         , cero_lmi_vector'         , cero_lmi_vector'         ;...
                                         h*N_3'     , -h*Z1_3       , cero_lmi      , cero_lmi      ;...    
                                         cero_lmi_vector    , cero_lmi      , -h*Z1_3       , cero_lmi      ;...
                                         cero_lmi_vector     , cero_lmi      , cero_lmi      , -h*Z2_3       ];

                        Condition_3 <= -eye(7*n_lmi);
                        %Condition1 <= 0;
                        P3      >= eye(n_lmi);
                        Q3      >= 0;
                        R3      >= 0;
                        Z1_3    >= eye(n_lmi);
                        Z2_3    >= eye(n_lmi);
            %4
                        N_4 = [N1_4 N2_4 N3_4 N4_4]';



                        aux_4 = reshape(A_glob(4,:,:),[n_lmi,n_lmi]);
                        bux_4 = reshape(B_glob(4,:,:),[n_lmi,n_lmi]);
                        A_lmi_4 = [-aux_4, -bux_4, cero_lmi, eye(n_lmi)];

                        O1_11_4 = Q4 + R4;
                        O1_14_4 = P4;
                        O1_22_4 = -1*(1-mu)*Q4;
                        O1_33_4 = -1*R4;
                        O1_44_4 = h*(Z1_4+Z2_4);

                        O1_4 = [O1_11_4     , cero_lmi  , cero_lmi  , O1_14_4       ;...
                              cero_lmi      , O1_22_4  , cero_lmi  , cero_lmi      ;...
                              cero_lmi      , cero_lmi  , O1_33_4   , cero_lmi      ;...
                              O1_14_4'      , cero_lmi  , cero_lmi  , O1_44_4       ];

                        O2_4 = [N_4, -N_4, cero_lmi_vector', cero_lmi_vector'] + (T*A_lmi_4) + ((A_lmi_4')*(T'));

                        O_4 = O1_4 + O2_4 + O2_4';

                        %Definitions Condition 1

                        Condition_4 =    [O_4       , h*N_4         , cero_lmi_vector'         , cero_lmi_vector'         ;...
                                         h*N_4'     , -h*Z1_4       , cero_lmi      , cero_lmi      ;...    
                                         cero_lmi_vector     , cero_lmi      , -h*Z1_4       , cero_lmi      ;...
                                         cero_lmi_vector     , cero_lmi      , cero_lmi      , -h*Z2_4       ];

                        Condition_4 <= -eye(7*n_lmi);
                        %Condition1 <= 0;
                        P4      >= eye(n_lmi);
                        Q4      >= 0;
                        R4      >= 0;
                        Z1_4    >= eye(n_lmi);
                        Z2_4    >= eye(n_lmi);
            %5
                        N_5 = [N1_5 N2_5 N3_5 N4_5]';



                        aux_5 = reshape(A_glob(5,:,:),[n_lmi,n_lmi]);
                        bux_5 = reshape(B_glob(5,:,:),[n_lmi,n_lmi]);
                        A_lmi_5 = [-aux_5, -bux_5, cero_lmi, eye(n_lmi)];

                        O1_11_5 = Q5 + R5 ;
                        O1_14_5 = P5;
                        O1_22_5 = -1*(1-mu)*Q5;
                        O1_33_5 = -1*R5;
                        O1_44_5 = h*(Z1_5+Z2_5);

                        O1_5 = [O1_11_5     , cero_lmi  , cero_lmi  , O1_14_5       ;...
                              cero_lmi      , O1_22_5   , cero_lmi  , cero_lmi      ;...
                              cero_lmi      , cero_lmi  , O1_33_5   , cero_lmi      ;...
                              O1_14_5'      , cero_lmi  , cero_lmi  , O1_44_5       ];

                        O2_5 = [N_5, -N_5, cero_lmi_vector', cero_lmi_vector'] + (T*A_lmi_5) + ((A_lmi_5')*(T'));

                        O_5 = O1_5 + O2_5 + O2_5';

                        %Definitions Condition 1

                        Condition_5 =    [O_5       , h*N_5       , cero_lmi_vector'         , cero_lmi_vector'         ;...
                                         h*N_5'     , -h*Z1_5      , cero_lmi      , cero_lmi      ;...    
                                         cero_lmi_vector     , cero_lmi      , -h*Z1_5      , cero_lmi      ;...
                                         cero_lmi_vector     , cero_lmi      , cero_lmi      , -h*Z2_5       ];

                        Condition_5 <= -eye(7*n_lmi);
                        %Condition1 <= 0;
                        P5      >= eye(n_lmi);
                        Q5      >= 0;
                        R5      >= 0;
                        Z1_5    >= eye(n_lmi);
                        Z2_5    >= eye(n_lmi);
            %6
                        N_6 = [N1_6 N2_6 N3_6 N4_6]';



                        aux_6 = reshape(A_glob(6,:,:),[n_lmi,n_lmi]);
                        bux_6 = reshape(B_glob(6,:,:),[n_lmi,n_lmi]);
                        A_lmi_6 = [-aux_6, -bux_6, cero_lmi, eye(n_lmi)];

                        O1_11_6 = Q6 + R6;
                        O1_14_6 = P6;
                        O1_22_6 = -1*(1-mu)*Q6;
                        O1_33_6 = -1*R6;
                        O1_44_6 = h*(Z1_6+Z2_6);

                        O1_6 = [O1_11_6     , cero_lmi  , cero_lmi  , O1_14_6       ;...
                              cero_lmi      , O1_22_6   , cero_lmi  , cero_lmi      ;...
                              cero_lmi      , cero_lmi  , O1_33_6   , cero_lmi      ;...
                              O1_14_6'      , cero_lmi  , cero_lmi  , O1_44_6       ];

                        O2_6 = [N_6, -N_6, cero_lmi_vector', cero_lmi_vector'] + (T*A_lmi_6) + ((A_lmi_6')*(T'));

                        O_6 = O1_6 + O2_6 + O2_6';

                        %Definitions Condition 1

                        Condition_6 =    [O_6       , h*N_6         , cero_lmi_vector'        , cero_lmi_vector'         ;...
                                         h*N_6'     , -h*Z1_6       , cero_lmi      , cero_lmi      ;...    
                                         cero_lmi_vector     , cero_lmi      , -h*Z1_6       , cero_lmi      ;...
                                        cero_lmi_vector     , cero_lmi      , cero_lmi      , -h*Z2_6       ];

                        Condition_6 <= -eye(7*n_lmi);
                        %Condition1 <= 0;
                        P6      >= eye(n_lmi);
                        Q6      >= 0;
                        R6      >= 0;
                        Z1_6    >= eye(n_lmi);
                        Z2_6    >= eye(n_lmi);
            %7
                        N_7 = [N1_7 N2_7 N3_7 N4_7]';



                        aux_7 = reshape(A_glob(7,:,:),[n_lmi,n_lmi]);
                        bux_7 = reshape(B_glob(7,:,:),[n_lmi,n_lmi]);
                        A_lmi_7 = [-aux_7, -bux_7, cero_lmi, eye(n_lmi)];

                        O1_11_7 = Q7 + R7;
                        O1_14_7 = P7;
                        O1_22_7 = -1*(1-mu)*Q7;
                        O1_33_7 = -1*R7;
                        O1_44_7 = h*(Z1_7+Z2_7);

                        O1_7 = [O1_11_7     , cero_lmi  , cero_lmi  , O1_14_7       ;...
                              cero_lmi      , O1_22_7   , cero_lmi  , cero_lmi      ;...
                              cero_lmi      , cero_lmi  , O1_33_7   , cero_lmi      ;...
                              O1_14_7'      , cero_lmi  , cero_lmi  , O1_44_7       ];

                        O2_7 = [N_7, -N_7, cero_lmi_vector', cero_lmi_vector'] + (T*A_lmi_7) + ((A_lmi_7')*(T'));

                        O_7 = O1_7 + O2_7 + O2_7';

                        %Definitions Condition 1

                        Condition_7 =    [O_7       , h*N_7         , cero_lmi_vector'         , cero_lmi_vector'        ;...
                                         h*N_7'     , -h*Z1_7       , cero_lmi      , cero_lmi      ;...    
                                         cero_lmi_vector    , cero_lmi      , -h*Z1_7       , cero_lmi      ;...
                                         cero_lmi_vector     , cero_lmi      , cero_lmi      , -h*Z2_7       ];

                        Condition_7 <= -eye(7*n_lmi);
                        %Condition1 <= 0;
                        P7     >= eye(n_lmi);
                        Q7      >= 0;
                        R7      >= 0;
                        Z1_7    >= eye(n_lmi);
                        Z2_7    >= eye(n_lmi);
            %8
                        N_8 = [N1_8 N2_8 N3_8 N4_8]';



                        aux_8 = reshape(A_glob(8,:,:),[n_lmi,n_lmi]);
                        bux_8 = reshape(B_glob(8,:,:),[n_lmi,n_lmi]);
                        A_lmi_8 = [-aux_8, -bux_8, cero_lmi, eye(n_lmi)];

                        O1_11_8 = Q8 + R8;
                        O1_14_8 = P8;
                        O1_22_8 = -1*(1-mu)*Q8;
                        O1_33_8 = -1*R8;
                        O1_44_8 = h*(Z1_8+Z2_8);

                        O1_8 = [O1_11_8     , cero_lmi  , cero_lmi  , O1_14_8       ;...
                              cero_lmi      , O1_22_8   , cero_lmi  , cero_lmi      ;...
                              cero_lmi      , cero_lmi  , O1_33_8   , cero_lmi      ;...
                              O1_14_8'      , cero_lmi  , cero_lmi  , O1_44_8       ];

                        O2_8 = [N_8, -N_8, cero_lmi_vector', cero_lmi_vector'] + (T*A_lmi_8) + ((A_lmi_8')*(T'));

                        O_8 = O1_8 + O2_8 + O2_8';

                        %Definitions Condition 1

                        Condition_8 =    [O_8       , h*N_8         , cero_lmi_vector'        ,cero_lmi_vector'         ;...
                                         h*N_8'     , -h*Z1_8       , cero_lmi      , cero_lmi      ;...    
                                         cero_lmi_vector     , cero_lmi      , -h*Z1_8       , cero_lmi      ;...
                                         cero_lmi_vector     , cero_lmi      , cero_lmi      , -h*Z2_8       ];

                        Condition_8 <= -eye(7*n_lmi);
                        %Condition1 <= 0;
                        P8      >= eye(n_lmi);
                        Q8      >= 0;
                        R8      >= 0;
                        Z1_8    >= eye(n_lmi);
                        Z2_8    >= eye(n_lmi);
            cvx_end

            status2008 = cvx_status;
    end % while feasible
end % for P=10e3 ...
 

