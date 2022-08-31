%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Taylor-Greene NS Solver with Ghost Nodes for Phi
% Solve rho[ u_t + u u_x + v u_y] =  - p_x + mu ( u_xx + u_yy) + fu
%       rho[ v_t + u v_x + v v_y] =  - p_y + mu ( v_xx + v_yy) + fv
% for x in [0, Lx], y in [0, Ly] with general BC's
%          a_i u + b_i Du = g_i
%          c_i v + d_i Dv = h_i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all;

% plotFlag:  1 -> FV Grid, 
%            2 -> Surf (includes error surf) during time stepping
%            3 -> Surf at Tf and Convergence Plot
%            4 -> Plots Stream Lines
% testFlag: 'temp'     ->  Temporal Test
%            'spatial'  -> Spatial Test with Steady State
% gridFlag:  0 -> uniform grid
%            1 -> nonuniform grid
% incr:      # of iterations between plots       
% BCflag:    0 -> Dirichlet Conditions
%            1 -> Neumann on bottom boundary
    
plotFlag = 3; 
testFlag = 'spatial';
gridFlag = 1;
incr = 10;
BCflag = 1;

% Store flags in cell array
flags = {plotFlag, gridFlag, incr}; 

% Domain H within [0, Wleft +Wmid +Wright] x [0, Hlow +Hmid +Hup]
Wleft = 2*pi/3; Wmid = 2*pi/3; Wright = 2*pi/3; 
Hlow  = 2*pi/3 ; Hmid = 2*pi/3; Hup    = 2*pi/3;

% Wleft = 1; Wmid = 1; Wright = 1; 
% Hlow  = 2; Hmid = 1; Hup    = 1;

mu = 2*pi;
rho = 1;
Uchar = 1; % m/s
Xchar = 1; % m
Pchar = rho * Uchar^2;
Tchar = Xchar/Uchar; % s
Fchar = rho*Uchar^2/Xchar;
Re = rho * Uchar * Xchar / mu

% Odds-n-ends from input parmeters
H = Hlow + Hmid + Hup;
W = Wleft + Wmid + Wright;
prms = [mu, rho, Re, Uchar, Xchar, Tchar, Pchar, Fchar]; % Store Parameters in vector


% Set Nx, Ny and dt for given test --------------------------------------
if strcmpi(testFlag,'temp') % Temporal Test Parameters
    % Set Range for Nt Values
    j = 3;
    Tf = 1; 
    dt = [0.01 0.005 0.0025 0.00125]; % Dirichlet Nx = 600

    %     dt = [0.001 0.0005 0.00025 0.000125]; % Dirichlet Nx = 600
    Nx = 600*ones(1,j+1); % 256 for quick results, 600 for best results
    Ny = Nx;
    j = 0;
    dt = 0.0001;
    Nx = 128; Ny = 128;

    
else % Spatial Tests
    % Set Range for N Values
    js = 1; % 1
    jf = 4; % 4 or 5
    j = jf - js;
    Nx = 8*2.^(js:jf);
    Ny = Nx;
%     Nx = 12;
%     Ny = 10;
    
    % Time
    dt = 0.001;
    dt = dt*ones(1,j+1);
    Tf = 8;

end

% Set Parameters h, beta, w
beta = 1;
h = 0;
if strcmpi(testFlag,'temp')  % set w depending on test
    w = 4*pi; % test temporal convergence
else 
    w = 0; % test spatial convergence
end


% Exact Solution and Derivatives -----------------------------------------
if BCflag == 0
    u = @(x,y,t) -cos(beta*x-h).*sin(beta*y-h)*cos(w*t);
    u_x = @(x,y,t) beta*sin(beta*x-h).*sin(beta*y-h)*cos(w*t);
    u_xx = @(x,y,t) beta^2*cos(beta*x-h).*sin(beta*y-h)*cos(w*t);
    u_y = @(x,y,t) -beta*cos(beta*x-h).*cos(beta*y-h)*cos(w*t);
    u_yy = @(x,y,t) beta^2*cos(beta*x-h).*sin(beta*y-h)*cos(w*t);

    v = @(x,y,t) cos(beta*y-h).*sin(beta*x-h)*cos(w*t);
    v_x = @(x,y,t) beta*cos(beta*y-h).*cos(beta*x-h)*cos(w*t);
    v_xx = @(x,y,t) -beta^2*cos(beta*y-h).*sin(beta*x-h)*cos(w*t);
    v_y = @(x,y,t) -beta*sin(beta*y-h).*sin(beta*x-h)*cos(w*t);
    v_yy = @(x,y,t) -beta^2*cos(beta*y-h).*sin(beta*x-h)*cos(w*t);

    u_t = @(x,y,t) w*cos(beta*x-h).*sin(beta*y-h)*sin(w*t);
    v_t = @(x,y,t) -w*sin(beta*x-h).*cos(beta*y-h)*sin(w*t);

    p = @(x,y,t) sin(beta*x-h).*sin(beta*y-h).*cos(w*t);
    p_x = @(x,y,t) beta*cos(beta*x-h).*sin(beta*y-h).*cos(w*t);
    p_y = @(x,y,t) beta*sin(beta*x-h).*cos(beta*y-h).*cos(w*t);
else
    u = @(x,y,t) cos(beta*y-h).*sin(beta*x-h)*cos(w*t);
    u_x = @(x,y,t) beta*cos(beta*y-h).*cos(beta*x-h)*cos(w*t);
    u_xx = @(x,y,t) -beta^2*cos(beta*y-h).*sin(beta*x-h)*cos(w*t);
    u_y = @(x,y,t) -beta*sin(beta*y-h).*sin(beta*x-h)*cos(w*t);
    u_yy = @(x,y,t) -beta^2*cos(beta*y-h).*sin(beta*x-h)*cos(w*t);
    
    v = @(x,y,t) -cos(beta*x-h).*sin(beta*y-h)*cos(w*t);
    v_x = @(x,y,t) beta*sin(beta*x-h).*sin(beta*y-h)*cos(w*t);
    v_xx = @(x,y,t) beta^2*cos(beta*x-h).*sin(beta*y-h)*cos(w*t);
    v_y = @(x,y,t) -beta*cos(beta*x-h).*cos(beta*y-h)*cos(w*t);
    v_yy = @(x,y,t) beta^2*cos(beta*x-h).*sin(beta*y-h)*cos(w*t);

    u_t = @(x,y,t) -w*sin(beta*x-h).*cos(beta*y-h)*sin(w*t);
    v_t = @(x,y,t) w*cos(beta*x-h).*sin(beta*y-h)*sin(w*t);

    p = @(x,y,t) sin(beta*x-h).*sin(beta*y-h).*cos(w*t);
    p_x = @(x,y,t) beta*cos(beta*x-h).*sin(beta*y-h).*cos(w*t);
    p_y = @(x,y,t) beta*sin(beta*x-h).*cos(beta*y-h).*cos(w*t);
end

% Set Initial Conditions -------------------------------------------------
if strcmpi(testFlag,'temp') % Test for Temporal Convergence
    ICu = @(x,y) u(Xchar*x,Xchar*y,0)/Uchar;
    ICv = @(x,y) v(Xchar*x,Xchar*y,0)/Uchar;
    ICp = @(x,y) p(Xchar*x,Xchar*y,0)/Pchar;
else % Test for spatial convergence (Steady State)
    ICu = @(x,y) 0*x.*y;
    ICv = @(x,y) 0*x.*y;
    ICp = @(x,y) 0*p(x,y,0);
end

% Store IC's in cell array for passing to solver
IC = {ICu ; ICv; ICp};

% Exact Solution
Soln = {u, v, p};

% Set Boundary Conditions for u-equation ------------------------------
% On \partial \Omega_i: a_i Up + b_i Up_y = g_i, for i = 2, 4
if BCflag == 0  
    a1 = 1; a2 = 1; a3 = 1; a4 = 1; 
    b1 = 0; b2 = 0; b3 = 0; b4 = 0;
else
    a1 = 1; a2 = 0; a3 = 1; a4 = 0; 
    b1 = 0; b2 = 1; b3 = 0; b4 = 1;
end


% Set Boundary Conditions for v-equation ------------------------------ 
% On \partial \Omega_i: c_i * Vp + d_i * Vp_y = h_i, for i = 2, 4
if BCflag == 0
    c1 = 1; c2 = 1; c3 = 1; c4 = 1; 
    d1 = 0; d2 = 0; d3 = 0; d4 = 0;
else
    c1 = 1; c2 = 0; c3 = 1; c4 = 0; 
    d1 = 0; d2 = 1; d3 = 0; d4 = 1;
end
    
    
% RHS BC's for u-equation -------------
g1 = @(x,t) u(Xchar*x,H,Tchar*t)/Uchar; % Top Left
g2 = @(x,t) a2*u(Xchar*x,0,Tchar*t)/Uchar ...
            + b2*Xchar*u_y(Xchar*x,0,Tchar*t)/Uchar; % Bottom Left
g3 = @(x,t) u(Xchar*x,H,Tchar*t)/Uchar; % Top Right
g4 = @(x,t) a4*u(Xchar*x,0,Tchar*t)/Uchar ...
            + b4*Xchar*u_y(Xchar*x,0,Tchar*t)/Uchar; % Bottom Right
g5 = @(y,t) u(0,Xchar*y,Tchar*t)/Uchar; % Left
g6 = @(y,t) u(W,Xchar*y,Tchar*t)/Uchar; % Right

wu1 = @(y,t) u(Wleft,Xchar*y,Tchar*t)/Uchar; % Wall 1
wu2 = @(y,t) u(Wleft,Xchar*y,Tchar*t)/Uchar; % Wall 2
wu3 = @(y,t) u(Wleft + Wmid,Xchar*y,Tchar*t)/Uchar; % Wall 3
wu4 = @(y,t) u(Wleft + Wmid,Xchar*y,Tchar*t)/Uchar; % Wall 4
wu5 = @(x,t) u(Xchar*x,Hlow + Hmid,Tchar*t)/Uchar; % Wall 5
wu6 = @(x,t) u(Xchar*x,Hlow,Tchar*t)/Uchar; % Wall 6

% RHS BC's for v-equation --------------
h1 = @(x,t) v(Xchar*x,H,Tchar*t)/Uchar; % Top Left
h2 = @(x,t) c2*v(Xchar*x,0,Tchar*t)/Uchar ...
            + d2*Xchar*v_y(Xchar*x,0,Tchar*t)/Uchar; % Bottom Left
h3 = @(x,t) v(Xchar*x,H,Tchar*t)/Uchar; % Top Right
h4 = @(x,t) c4*v(Xchar*x,0,Tchar*t)/Uchar ...
            + d4*Xchar*v_y(Xchar*x,0,Tchar*t)/Uchar; % Bottom Right
h5 = @(y,t) v(0,Xchar*y,Tchar*t)/Uchar; % Left
h6 = @(y,t) v(W,Xchar*y,Tchar*t)/Uchar; % Right

wv1 = @(y,t) v(Wleft,Xchar*y,Tchar*t)/Uchar; % Wall 1
wv2 = @(y,t) v(Wleft,Xchar*y,Tchar*t)/Uchar; % Wall 2
wv3 = @(y,t) v(Wleft + Wmid,Xchar*y,Tchar*t)/Uchar; % Wall 3
wv4 = @(y,t) v(Wleft + Wmid,Xchar*y,Tchar*t)/Uchar; % Wall 4
wv5 = @(x,t) v(Xchar*x,Hlow + Hmid,Tchar*t)/Uchar; % Wall 5
wv6 = @(x,t) v(Xchar*x,Hlow,Tchar*t)/Uchar; % Wall 6

% Create cell and vector for Boundary Conditions 
BCu = {g1,  g2,  g3,  g4,  g5,  g6;
       wu1, wu2, wu3, wu4, wu5, wu6};
      
BCv = {h1,  h2,  h3,  h4,  h5,  h6;
       wv1, wv2, wv3, wv4, wv5, wv6};

% Create Vectors containing coefficients for boundary conditions
u_coeffs = [a1 a2 a3 a4;
            b1 b2 b3 b4];
v_coeffs = [c1 c2 c3 c4; 
            d1 d2 d3 d4];

% Set Body Force Terms fx and fy ------------------------------------------
          
fx = @(x,y,t)  1/Fchar * (p_x(Xchar*x,Xchar*y,Tchar*t)... 
               + rho*(u_t(Xchar*x,Xchar*y,Tchar*t) ...
            + u(Xchar*x,Xchar*y,Tchar*t).*u_x(Xchar*x,Xchar*y,Tchar*t) ...
            + v(Xchar*x,Xchar*y,Tchar*t).*u_y(Xchar*x,Xchar*y,Tchar*t))... 
     - mu*(u_xx(Xchar*x,Xchar*y,Tchar*t) + u_yy(Xchar*x,Xchar*y,Tchar*t)));


fy = @(x,y,t)  1/Fchar * (p_y(Xchar*x,Xchar*y,Tchar*t)...
               + rho*(v_t(Xchar*x,Xchar*y,Tchar*t)...
            + u(Xchar*x,Xchar*y,Tchar*t).*v_x(Xchar*x,Xchar*y,Tchar*t)...
            + v(Xchar*x,Xchar*y,Tchar*t).*v_y(Xchar*x,Xchar*y,Tchar*t))...
     - mu*(v_xx(Xchar*x,Xchar*y,Tchar*t) + v_yy(Xchar*x,Xchar*y,Tchar*t)));

f = {fx ; fy}; % Store in Cell Array for Passing to Solver

% Initialize Err vector (for convergence tests)
uErr = zeros(j+1,1);
vErr = zeros(j+1,1);
pErr = zeros(j+1,1);
S = cell(j+1,3);
grd = cell(j+1,6);
N = zeros(j+1,6);

% Get Approximation(s) ---------------------------------------------------
for i = 1:j+1
    D = [Nx(i),Ny(i),Wleft,Wmid,Wright,Hlow,Hmid,Hup]; % Domain Params.
    t = 0:dt(i):Tf;

    [err,s,g,n] = SolveH_NS(D,t,prms,IC,u_coeffs,v_coeffs,BCu,BCv,f,flags,Soln);
    uErr(i) = err(1); vErr(i) = err(2); pErr(i) = err(3);
    S(i,:) = s;
    grd(i,:) = g;
    N(i,:) = n;
end

save('data.mat')

% Create Convergence Table 
uROC_ex = EOC_ex(uErr');
uROC_ex = [NaN; uROC_ex'];
vROC_ex = EOC_ex(vErr');
vROC_ex = [NaN; vROC_ex'];
pROC_ex = EOC_ex(pErr');
pROC_ex = [NaN; pROC_ex'];

% uROC_aprx = EOC_aprx(S(:,1),grd(:,1),grd(:,2));
% uROC_aprx = [NaN; NaN; uROC_aprx];
% vROC_aprx = EOC_aprx(S(:,2),grd(:,3),grd(:,4));
% vROC_aprx = [NaN; NaN; vROC_aprx];
% pROC_aprx = EOC_aprx(S(:,3),grd(:,5),grd(:,6));
% pROC_aprx = [NaN; NaN; pROC_aprx];

% [uROC_aprx,vROC_aprx,pROC_aprx] = apprx_conv(S,grd,N);


Nx = Nx';
dt = dt';
table(Nx, dt, uErr, uROC_ex,vErr, vROC_ex, pErr, pROC_ex)

% table(Nx, dt, uErr, uROC_ex, uROC_aprx, vErr, vROC_ex,vROC_aprx, pErr, pROC_ex)



% Plot Errors
if plotFlag == 2 || 3
    figure(2)
    M = 10;
    if strcmpi(testFlag,'temp') == 1
        Nt = Tf./dt;
        line2 = dt.^(2)*exp(5);
        loglog(dt, uErr, 'bs', 'MarkerSize', M,'DisplayName', 'Error in U')
        hold on
        loglog(dt, vErr, 'r^', 'MarkerSize', M,'DisplayName', 'Error in V') 
        hold on
        loglog(dt, pErr, 'k*', 'MarkerSize', M,'DisplayName', 'Error in P') 
        hold on
        loglog(dt, line2, '-.','DisplayName', '\Delta t^2')
        xlabel('\Delta t','fontsize', 14)
        legend('fontsize', 16)
        ylabel('Error','fontsize', 14)
        axis('tight')
        title(['Error at $t = $', num2str(Tf),' with $\beta =$ ',...
            num2str(beta), ', $h = $ ',num2str(h),', $\omega = $', num2str(w)],...
            'interpreter', 'latex', 'fontsize', 18)
        set(gca,'LooseInset',get(gca,'TightInset'));
    else
        line = Nx.^(-2)*exp(2);
        loglog(Nx, uErr, 'bs', 'MarkerSize', M,'DisplayName', 'Error in U')
        hold on
        loglog(Nx, vErr, 'r^', 'MarkerSize', M,'DisplayName', 'Error in V')
        hold on
        loglog(Nx, pErr, 'k*', 'MarkerSize', M,'DisplayName', 'Error in P')
        hold on
        loglog(Nx, line, '-.','DisplayName', 'N^{-2}')
        xlabel('N','fontsize', 14)
        legend('fontsize', 16)
        ylabel('Error','fontsize', 14)
        axis('tight')
        title(['Error at $t = $', num2str(Tf),' with $\beta =$ ',...
            num2str(beta), ', $h = $ ',num2str(h),', $\omega = $', num2str(w)],...
            'interpreter', 'latex', 'fontsize', 18)
        set(gca,'LooseInset',get(gca,'TightInset'));
    end
end





