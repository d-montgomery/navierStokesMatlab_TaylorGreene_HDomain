function [Err, S, grd, N] = SolveH_NS(D,t,prms,IC,u_coeffs,v_coeffs,BCu,BCv,f,flags,Soln)
% Solve rho[ u_t + u u_x + v u_y] = -p_x + mu ( u_xx + u_yy) + fx
%       rho[ v_t + u v_x + v v_y] = -p_y + mu ( v_xx + v_yy) + fy
% in the H domain given by parameters Hlow, Hmid, Hup, WL, Wmid, WR.
% Inputs:
% D    : Domain parameters D = [Nx,Ny,Wleft,Wmid,Wright,Hlow,Hmid,Hup]
% t    : is a row vector for time 0 < t <= Tf
% prms : a vector with parameters prms = [mu, rho, w, beta, h]
% IC   : a 3x1 cell array where IC{1} = ICu = @(x,y) ...;  for u - eqt
%                               IC{2} = ICv = @(x,y) ...;  for v - eqt.
%                               IC{3} = ICp = @(x,y) ...;  for pressure.
% u_coeffs : matrix containing BC coefficients for u-momentum eqt, 
%            associated with \partial \Omega_i, u_coeffs = [a5 a6; b5 b6];   
% v_coeffs : matrix containing BC coefficients for v-momentum eqt, 
%                                   v_coeffs = [c1 c2 c3 c4; d1 d2 d3 d4];
% BCu, BCv : cell arrays containing RHS functions for each BC.
%                                   BCu = {g1,  g2,  g3,  g4,  g5,  g6;
%                                         wu1, wu2, wu3, wu4, wu5, wu6};      
%                                   BCv = {h1,  h2,  h3,  h4,  h5,  h6;
%                                         wv1, wv2, wv3, wv4, wv5, wv6};
% f     : a 2X1 cell array with body forcing functions 
%                                  f{1} = fx = @(x,y,t) ...;  for u - eqt
%                                  f{2} = fy = @(x,y,t) ...;  for v - eqt
% flags : a cell array with flags = {PLOT, GRID, incr};
% Soln  : a cell array with the manufactured solution Soln = {u, v, p};

% unpack flags
PLOT = flags{1};
GRID = flags{2};
incr = flags{3};

% unpack physical parameters and char. scales
mu = prms(1);
rho = prms(2);
Re = prms(3);
Uchar = prms(4);
Xchar = prms(5);
Tchar = prms(6);
Pchar = prms(7);

% Unpack spatial parameters 
Nx = D(1);     
Ny = D(2);
D(3:8) = D(3:8)/Xchar; % Nondimensionalize lengths
WL = D(3); Wmid = D(4); WR = D(5); 
Hlow  = D(6); Hmid = D(7); Hup = D(8);

% Nondimensionalize time and gather time params Nt, dt, Tf
t = t / Tchar;
Nt = length(t) - 1;
dt = t(2)-t(1);
Tf = t(Nt+1);

% Unpack Initial Conds, Boundary Conds Coeffs, and Forcing Terms fx, fy
ICu = IC{1};
ICv = IC{2};
ICp = IC{3};

fx = f{1};
fy = f{2};

% Get grid for FVM 
[XU,YU,XV,YV,xp,yp,N] = NSGridHnoGhst(Nx,Ny,WL,Wmid,WR,Hlow,Hmid,Hup,GRID,PLOT);

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% Create Mesh Grid for Plotting
[XP, YP] = meshgrid(xp(2:Nx+1),yp(2:Ny+1)); %Ignore ghost cells

% Get Matrices Au, Av, Ap
Au = BuildAuH(N,XU,YU,XV,YV,u_coeffs,Re,dt);
Av = BuildAvH(N,XU,YU,XV,YV,v_coeffs,Re,dt);
Ap = BuildApH_Ghst(N,xp,yp,XU,YV,v_coeffs(2,2),v_coeffs(1,2));

% Store factorizations for efficient solving
Au = factorize(Au);
Av = factorize(Av);
Ap = factorize(Ap); 

% Initial Conditions (Includes Ghost Nodes)
U0 = ICu(XU,YU);
V0 = ICv(XV,YV);
P0 = ICp(XP,YP);

% Force zeros in middle to form H (doesn't zero ghost nodes)
[U0,V0,P0] = zeroH(U0,V0,P0,NxL,NyL,NxM,NyM,Ny);

% --- One step of Forward Euler ---
[NLU,NLV] = NonLin(U0,V0,N,XU,YU,XV,YV);
NLU = 2/3*NLU;
NLV = 2/3*NLV;
NLU0 = 0*NLU;
NLV0 = 0*NLV;

% Forcing Terms at t = t0
Quo = fx(XU,YU,t(1)); 
Qvo = fy(XV,YV,t(1)); 

if PLOT > 1
    figure(1)
end

% Step in time
tic

for k = 2:Nt+1
    
    % Forcing Terms
    Qu = fx(XU,YU,t(k)); 
    Qv = fy(XV,YV,t(k)); 
    
    % Exact Solution at current time step = t + dt
    Ue = Soln{1}( Xchar*XU,Xchar*YU,Tchar*t(k) );
    Ve = Soln{2}( Xchar*XV,Xchar*YV,Tchar*t(k) );
    Pe = Soln{3}( Xchar*XP,Xchar*YP,Tchar*t(k) );
    
    % Prediction - Solve u-momentum equation
    Fu = uRHS_H(U0,NLU,NLU0,P0,Qu,Quo,N,dt,Re,XU,YU,XV,YV,u_coeffs,BCu,t(k));
    TMPU = Au \ Fu; 

    % Prediction - Solve v- momentum equation
    Fv = vRHS_H(V0,NLV,NLV0,P0,Qv,Qvo,N,dt,Re,XU,YU,XV,YV,v_coeffs,BCv,t(k));
    TMPV = Av \ Fv;
    
    % Reshape TMPU and TMPV into Matrices
    [U,V] = reshapeUV_H(NxL,NxM,Nx,NyL,NyM,Ny,TMPU,TMPV);
     
    % Correction - Solve PHI Equation
    Fp = pRHS_H_Ghst(N,U,V,dt,XU,YV);
    TMPP = Ap\Fp;
    
    % Reshape TMPP into Matrix
    PHI = reshapePHI_Ghst(NxL,NxM,Nx,NyL,NyM,Ny,TMPP);

    % Correction - Get Current U, V and P
    [U,V,P] = Correction_Ghst(U,V,P0,PHI,N,dt,xp,yp);
    
    % Update Velocities
    U0 = U;
    V0 = V;
    P0 = P;
    Quo = Qu;
    Qvo = Qv;
    NLU0 = NLU;
    NLV0 = NLV;
    [NLU,NLV] = NonLin(U0,V0,N,XU,YU,XV,YV);

    % plot solution at various values of time
    if PLOT == 2 && ~mod(k,incr)
        DIV = mean( mean( abs(div(N,U,V,XU,YV) ) ) ) % check div(u) = 0
        
        % Force zeros in middle to form H
        [U,V,P] = zeroH(U,V,P,NxL,NyL,NxM,NyM,Ny);
        [Ue,Ve,Pe] = zeroH(Ue,Ve,Pe,NxL,NyL,NxM,NyM,Ny);
        makePlots(Uchar*U,Uchar*V,Pchar*P,Ue,Ve,Pe,Xchar*XU,Xchar*YU,...
                  Xchar*XV,Xchar*YV,Xchar*XP,Xchar*YP,Tchar*t(k))
        pause(0.1)
    end
end
toc

% Force zeros in middle to form H
% [U,V,P] = zeroH(U,V,P,NxL,NyL,NxM,NyM,Ny);
[Ue,Ve,Pe] = zeroH(Ue,Ve,Pe,NxL,NyL,NxM,NyM,Ny);
    
% Create and save Surf Plot at Tf
if PLOT == 3
    makePlots(Uchar*U,Uchar*V,Pchar*P,Ue,Ve,Pe,Xchar*XU,Xchar*YU,...
                  Xchar*XV,Xchar*YV,Xchar*XP,Xchar*YP,Tchar*t(k))
end


% Final Error in u-equation
xAbsErr = abs(Ue - Uchar*U);
uErr = max( max(xAbsErr) )/ max( max(abs(Ue)) ); % Store error at Tf

% Final Error in v-equation
yAbsErr = abs(Ve - Uchar*V);
vErr = max( max(yAbsErr) )/ max( max(abs(Ve)) ); % Store error at Tf

% Final Error in P
pAbsErr = abs(Pe - Pchar*P);
pErr = max( max(pAbsErr) )/ max( max(abs(Pe)) ); % Store error at Tf

Err = [uErr, vErr, pErr];

S = {Uchar*U, Uchar*V, Pchar*P};
grd = {Xchar*XU,Xchar*YU,Xchar*XV,Xchar*YV,Xchar*XP,Xchar*YP};



