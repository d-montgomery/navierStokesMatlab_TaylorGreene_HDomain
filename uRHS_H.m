function[F] = uRHS_H(U0,NLU,NLU0,P0,Qu,Quo,N,dt,Re,XU,YU,XV,YV,ab,BC,t)
% INPUTS:
% U0 = Previous velocity U
% NLU = Nonlinear terms at t = t_n
% NLU0 = Nonlinear terms at t = t_{n-1}
% P0 = Pressure at t = t_n
% Qu = Body forcing terms at t = t_{n+1}
% Quo = Body forcing terms at t = t_n
% N = [NxL, NxM, Nx, NyL, NyM, Ny]
% dt = temporal step size
% Re = Reynolds number
% xu and yu is where u(x,y) is stored on grid
% xv and yv is where v(x,y) is stored on grid
% ab are the coefficients for the inlet/outlet BC's
% BC = cell array with vectors containing boundary conditions
% t = t^{n+1}

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% --- Matrix Sizes ---
% Au is NNXY x NNXY matrix, Fu is an NNXY x 1 vector
NNXY = (Nx+1)*(Ny+2) - (NxM - NxL -1) * ( Ny - NyM + NyL);
% Size of sml sub matrices = Nsub x Nsub
Nsml = (NxL + 1) + (Nx + 1 - NxM); 

F = zeros( NNXY,1 );

% --- Boundary Conditions ------------------------------------------------
% Unpack coeffs. for BCs: u_coeffs = [a5 a6; b5 b6];
a5 = ab(1,1); a6 = ab(1,2);
b5 = ab(2,1); b6 = ab(2,2);

% RHS for BC
g1 = BC{1,1}(XU(Ny+2,:),t);   g2 = BC{1,2}(XU(1,:),t); 
g3 = BC{1,3}(XU(Ny+2,:),t);   g4 = BC{1,4}(XU(1,:),t); 
g5 = BC{1,5}(YU(:,1),t);  g6 = BC{1,6}(YU(:,Nx+1),t);

w1 = BC{2,1}(YU(:,NxL+1),t);  w2 = BC{2,2}(YU(:,NxL+1),t); 
w3 = BC{2,3}(YU(:,NxM+1),t);  w4 = BC{2,4}(YU(:,NxM+1),t); 
w5 = BC{2,5}(XU(NyM+2,:),t);   w6 = BC{2,6}(XU(NyL+1,:),t);

% Top Left BC (\partial \Omega_1)
for j = 1:NxL+1
    F(NNXY-Nsml+j) = g1(j); 
end

% Bottom Left BC (\partial \Omega_2)
for j = 1:NxL+1
    F(j) = g2(j); 
end

% Top Right BC (\partial \Omega_3)
j_ctr = NxL+1;
for j = NxM+1 : Nx+1
    j_ctr = j_ctr + 1;
    F(NNXY-Nsml+j_ctr) = g3(j);
end

% Bottom Right BC (\partial \Omega_4)
j_ctr = NxL+1;
for j = NxM+1 : Nx+1
    j_ctr = j_ctr + 1;
    F(j_ctr) = g4(j); 
end

% Left and Right BC's (\partial \Omega_i, for i = 5, 6)
for i = 2 : Ny+1
    if i < NyL + 1 % Lower Legs
        ROW_L = (i-1)*Nsml + 1;
        ROW_R = (i)*Nsml;
    elseif i > NyL && i < NyM + 3 % Middle of H
        ROW_L = NyL*Nsml + (i - NyL - 1) * (Nx+1) + 1;
        ROW_R = NyL*Nsml + (i - NyL) * (Nx+1);
    else % Upper Legs
        ROW_L = NyL*Nsml + (NyM+2 - NyL)*(Nx+1) + (i-NyM - 3)*Nsml + 1;
        ROW_R = NyL*Nsml + (NyM+2 - NyL)*(Nx+1) + (i-NyM - 2)*Nsml;
    end
    F(ROW_L) = g5(i); % Left \partial \Omega_5
    F(ROW_R) = g6(i); % Right \partial \Omega_6
end

% Left and Right Walls for Upper Legs (\Gamma_i, i = 1,3)
for i = NyM + 2 : Ny + 1
    if i == NyM + 2
        ROW_L = NyL*Nsml + (NyM+1 - NyL)*(Nx+1) + NxL + 1;
        ROW_R = NyL*Nsml + (NyM+1 - NyL)*(Nx+1) + NxM + 1;
    else
        ROW_L = NyL*Nsml + (NyM+2 - NyL)*(Nx+1) ...
                + (i-NyM - 3)*Nsml + NxL + 1;
        ROW_R = NyL*Nsml + (NyM+2 - NyL)*(Nx+1) ...
                + (i-NyM - 3)*Nsml + NxL + 2;
    end

    F(ROW_L) = w1(i); % Left \Gamma_1
    F(ROW_R) = w3(i); % Left \Gamma_3
end

% Left and Right Walls for Lower Legs (\Gamma_i, i = 2,4)
for i = 2:NyL+1
    if i < NyL + 1
        ROW_L = (i-1)*Nsml + NxL + 1;
        ROW_R = (i-1)*Nsml + NxL + 2;
    else
        ROW_L = (i-1)*Nsml + NxL + 1;
        ROW_R = (i-1)*Nsml + NxM + 1;
    end

    F(ROW_L) = w2(i); % Left \Gamma_2
    F(ROW_R) = w4(i); % Left \Gamma_4
end

% Upper Mid-H Wall \Gamma_5
for j = NxL + 2 : NxM
    ROW = NyL*Nsml + (NyM + 1 - NyL) * (Nx+1) + j;
    F(ROW) = w5(j); % Left \Gamma_5
end

% Lower Mid-H Wall \Gamma_6
for j = NxL + 2 : NxM
    ROW = NyL*Nsml + j;
    F(ROW) = w6(j); % Left \Gamma_5
end
% --- End of Boundary Conditions ----------------------------------------


% ---- Fill Inner Cells in Lower Legs of H ------------------------------
for i = 2:NyL-1
    j_ctr = NxL+2;
    
for j = [2:NxL, NxM+2:Nx]
    
    if j < NxM+2
        ROW = Nsml*(i-1)+j;
    else 
        j_ctr = j_ctr + 1;
        ROW = Nsml*(i-1)+j_ctr;
    end
    
    % Cell height/width
    dx = XV(i,j+1) - XV(i,j);
    dy = YV(i,j) - YV(i-1,j);
    
    % Location of Nodes
    XE = XU(i,j+1);  XP = XU(i,j);  XW = XU(i,j-1);
    YN = YU(i+1,j);  YP = YU(i,j);  YS = YU(i-1,j);
    
    % Velocities           
    UP = U0(i,j);
    UE = U0(i,j+1); UN = U0(i+1,j); UW = U0(i,j-1); US = U0(i-1,j);


    % Diffusion Terms
    dw = - 1/Re * (UP - UW)*dy/(XP-XW)/2;
    de = 1/Re * (UE - UP)*dy/(XE-XP)/2;
    ds = -1/Re * (UP - US)*dx/(YP-YS)/2;
    dn = 1/Re * (UN - UP)*dx/(YN-YP)/2;
    D = de + dn + dw + ds;

    % Pressure Terms
    Pe = P0(i-1,j); Pw = P0(i-1,j-1);
    PR =  (Pe-Pw)* dy;
    F( ROW ) = U0(i,j)*dx*dy/dt - 1.5*NLU(i,j) +...
                              0.5*NLU0(i,j) + D - PR + ...
                             + (Qu(i,j) + Quo(i,j))*dx*dy/2;
  
end  
end

% ---- Fill Inner cells at Lower Interface yu_i, i = NyL, NyL+1 ---
for i = NyL : NyL+1
    j_ctr = NxL+2; % Starting matrix value right H when i == NyL only 
    
    for j = [2:NxL, NxM+2:Nx]
        if (i == NyL) && j > NxM+1
            j_ctr = j_ctr + 1;
            ROW = Nsml*(i-1)+j_ctr;
        else
            ROW = Nsml*(i-1)+j;
        end
        
        % Cell height/width
        dx = XV(i,j+1) - XV(i,j);
        dy = YV(i,j) - YV(i-1,j);

        % Location of Nodes
        XE = XU(i,j+1);  XP = XU(i,j);  XW = XU(i,j-1);
        YN = YU(i+1,j);  YP = YU(i,j);  YS = YU(i-1,j);

        % Velocities           
        UP = U0(i,j);
        UE = U0(i,j+1); UN = U0(i+1,j); UW = U0(i,j-1); US = U0(i-1,j);


        % Diffusion Terms
        dw = -1/Re * (UP - UW)*dy/(XP-XW)/2;
        de = 1/Re * (UE - UP)*dy/(XE-XP)/2;
        ds = -1/Re * (UP - US)*dx/(YP-YS)/2;
        dn = 1/Re * (UN - UP)*dx/(YN-YP)/2;
        D = de + dn + dw + ds;

        % Pressure Terms
        Pe = P0(i-1,j); Pw = P0(i-1,j-1);
        PR =  (Pe-Pw)* dy;

        F( ROW ) = U0(i,j)*dx*dy/dt - 1.5*NLU(i,j) +...
                                  0.5*NLU0(i,j) + D - PR + ...
                                 + (Qu(i,j) + Quo(i,j))*dx*dy/2;

    end 
end


% ---- Fill inner cells in middle of H for i = NyL+2, ..., NyM+1 ---
for i = NyL+2 : NyM+1
    for j = 2:Nx
        
        ROW = Nsml*NyL + (i-NyL-1)*(Nx+1) + j;
        
        % Cell height/width
        dx = XV(i,j+1) - XV(i,j);
        dy = YV(i,j) - YV(i-1,j);

        % Location of Nodes
        XE = XU(i,j+1);  XP = XU(i,j);  XW = XU(i,j-1);
        YN = YU(i+1,j);  YP = YU(i,j);  YS = YU(i-1,j);

        % Velocities           
        UP = U0(i,j);
        UE = U0(i,j+1); UN = U0(i+1,j); UW = U0(i,j-1); US = U0(i-1,j);


        % Diffusion Terms
        dw = -1/Re * (UP - UW)*dy/(XP-XW)/2;
        de = 1/Re * (UE - UP)*dy/(XE-XP)/2;
        ds = -1/Re * (UP - US)*dx/(YP-YS)/2;
        dn = 1/Re * (UN - UP)*dx/(YN-YP)/2;
        D = de + dn + dw + ds;

        % Pressure Terms
        Pe = P0(i-1,j); Pw = P0(i-1,j-1);
        PR =  (Pe-Pw)* dy;

        F( ROW ) = U0(i,j)*dx*dy/dt - 1.5*NLU(i,j) +...
                                  0.5*NLU0(i,j) + D - PR + ...
                                 + (Qu(i,j) + Quo(i,j))*dx*dy/2;

    end 
end

% ---- Fill Inner cells at Upper Interface yu_i, i = NyM+2, NyM+3 ----
for i =  NyM+2 : NyM+3
    j_ctr = NxL+2; % Starting matrix value right H when i == NyM+3 only 

    for j = [2:NxL, NxM+2:Nx]
        if (i == NyM+3) && (j > NxM+1)
            j_ctr = j_ctr + 1;
            ROW = Nsml*NyL + (NyM+2-NyL)*(Nx+1) + j_ctr;
        else
            ROW = Nsml*NyL + (i-NyL-1)*(Nx+1) + j;
        end
        
        % Cell height/width
        dx = XV(i,j+1) - XV(i,j);
        dy = YV(i,j) - YV(i-1,j);

        % Location of Nodes
        XE = XU(i,j+1);  XP = XU(i,j);  XW = XU(i,j-1);
        YN = YU(i+1,j);  YP = YU(i,j);  YS = YU(i-1,j);

        % Velocities           
        UP = U0(i,j);
        UE = U0(i,j+1); UN = U0(i+1,j); UW = U0(i,j-1); US = U0(i-1,j);


        % Diffusion Terms
        dw = -1/Re * (UP - UW)*dy/(XP-XW)/2;
        de = 1/Re * (UE - UP)*dy/(XE-XP)/2;
        ds = -1/Re * (UP - US)*dx/(YP-YS)/2;
        dn = 1/Re * (UN - UP)*dx/(YN-YP)/2;
        D = de + dn + dw + ds;

        % Pressure Terms
        Pe = P0(i-1,j); Pw = P0(i-1,j-1);
        PR =  (Pe-Pw)* dy;

        F( ROW ) = U0(i,j)*dx*dy/dt - 1.5*NLU(i,j) +...
                                  0.5*NLU0(i,j) + D - PR + ...
                                 + (Qu(i,j) + Quo(i,j))*dx*dy/2;
    end 
end

% ---- Fill Inner Cells in Upper Legs of H ------------------------------
for i = NyM+4 : Ny+1
    j_ctr = NxL+2; % Starting matrix value right H when i == NyM+3 only 

    for j = [2:NxL, NxM+2:Nx]
        if j < NxM+2
            ROW = Nsml*NyL + (NyM+2-NyL)*(Nx+1) + (i-NyM-3)*Nsml + j;
        else 
            j_ctr = j_ctr + 1;
            ROW = Nsml*NyL + (NyM+2-NyL)*(Nx+1) + (i-NyM-3)*Nsml + j_ctr;
        end
        
        % Cell height/width
        dx = XV(i,j+1) - XV(i,j);
        dy = YV(i,j) - YV(i-1,j);

        % Location of Nodes
        XE = XU(i,j+1);  XP = XU(i,j);  XW = XU(i,j-1);
        YN = YU(i+1,j);  YP = YU(i,j);  YS = YU(i-1,j);

        % Velocities           
        UP = U0(i,j);
        UE = U0(i,j+1); UN = U0(i+1,j); UW = U0(i,j-1); US = U0(i-1,j);


        % Diffusion Terms
        dw = -1/Re * (UP - UW)*dy/(XP-XW)/2;
        de = 1/Re * (UE - UP)*dy/(XE-XP)/2;
        ds = -1/Re * (UP - US)*dx/(YP-YS)/2;
        dn = 1/Re * (UN - UP)*dx/(YN-YP)/2;
        D = de + dn + dw + ds;

        % Pressure Terms
        Pe = P0(i-1,j); Pw = P0(i-1,j-1);
        PR =  (Pe-Pw)* dy;

        F( ROW ) = U0(i,j)*dx*dy/dt - 1.5*NLU(i,j) +...
                                  0.5*NLU0(i,j) + D - PR + ...
                                 + (Qu(i,j) + Quo(i,j))*dx*dy/2;
    end 
end