function[F] = vRHS_H(V0,NLV,NLV0,P0,Qv,Qvo,N,dt,Re,XU,YU,XV,YV,cd,BC,t)
% INPUTS:
% U0 = Previous velocity U
% NLU = Nonlinear terms at t = t_n
% NLU0 = Nonlinear terms at t = t_{n-1}
% P0 = Pressure at t = t_n
% Qv = Body forcing terms at t = t_{n+1}
% Qvo = Body forcing terms at t = t_n
% N = [NxL, NxM, Nx, NyL, NyM, Ny]
% dt = temporal step size
% Re = Reynolds number
% xu and yu is where u(x,y) is stored on grid
% xv and yv is where v(x,y) is stored on grid
% cd are the coefficients for the inlet/outlet BC's
% BC = cell array with vectors containing boundary conditions
% t = t^{n+1}


% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% --- Matrix Sizes ---
% Av is NNXY x NNXY matrix
NNXY = (Nx+2)*(Ny+1) - (NxM - NxL -2) * ( Ny - NyM + NyL);
% Size of sml sub matrices = Nsub x Nsub
Nsml = (NxL + 2) + (Nx + 2 - NxM); 

F = zeros( NNXY,1 );

% --- Boundary Conditions ------------------------------------------------
% Unpack coeffs. for BCs: v_coeffs = [c1 c2 c3 c4; d1 d2 d3 d4];
c1 = cd(1,1); c2 = cd(1,2); c3 = cd(1,3); c4 = cd(1,4);
d1 = cd(2,1); d2 = cd(2,2); d3 = cd(2,3); d4 = cd(2,4);

% RHS for BC
h1 = BC{1,1}(XV(Ny+1,:),t);   h2 = BC{1,2}(XV(1,:),t); 
h3 = BC{1,3}(XV(Ny+1,:),t);   h4 = BC{1,4}(XV(1,:),t); 
h5 = BC{1,5}(YV(:,1),t);  h6 = BC{1,6}(YV(:,Nx+2),t);

w1 = BC{2,1}(YV(:,NxL+2),t);  w2 = BC{2,2}(YV(:,NxL+2),t); 
w3 = BC{2,3}(YV(:,NxM+1),t);  w4 = BC{2,4}(YV(:,NxM+1),t); 
w5 = BC{2,5}(XV(NyM+1,:),t);   w6 = BC{2,6}(XV(NyL+1,:),t);

% Typical Velocities
VinL = min( h1(2:NxL+1) );
VinR = min( h3(NxM+2:Nx+1) );


% Top Left BC (\partial \Omega_1)
for j = 1:NxL+2
    F(NNXY-Nsml+j) = h1(j); 
end

% Bottom Left BC (\partial \Omega_2)
for j = 1:NxL+2
    F(j) = h2(j);
end

% Top Right BC (\partial \Omega_3)
j_ctr = NxL+2; % indexing for right bc starts at NxL + 4
for j = NxM+1 : Nx+2
    j_ctr = j_ctr + 1;
    F(NNXY-Nsml+j_ctr) = h3(j);
end

% Bottom Right BC (\partial \Omega_4)
j_ctr = NxL+2;
for j = NxM+1 : Nx+2
    j_ctr = j_ctr + 1;
    F(j_ctr) = h4(j);
end

% Left and Right BC's (\partial \Omega_i, for i = 5, 6)
for i = 2 : Ny
    if i < NyL + 1 % Lower Legs i = 1:NyL
        ROW_L = (i-1)*Nsml + 1;
        ROW_R = (i)*Nsml;
    elseif (i > NyL) && (i < NyM + 2) % Middle of H i = NyL+1:NyM+1
        ROW_L = NyL*Nsml + (i - NyL - 1) * (Nx+2) + 1;
        ROW_R = NyL*Nsml + (i - NyL) * (Nx+2);
    else % Upper Legs i = NyM+2:Ny+1
        ROW_L = NyL*Nsml + (NyM+1 - NyL)*(Nx+2) + (i-NyM - 2)*Nsml + 1;
        ROW_R = NyL*Nsml + (NyM+1 - NyL)*(Nx+2) + (i-NyM - 1)*Nsml;
    end
    F(ROW_L) = h5(i); % Left \partial \Omega_5
    F(ROW_R) = h6(i); % Right \partial \Omega_6
end

% Left and Right Walls for Upper Legs (\Gamma_i, i = 1,3)
for i = NyM + 2 : Ny 
    ROW_L = NyL*Nsml + (NyM+1 - NyL)*(Nx+2) ...
                + (i-NyM - 2)*Nsml + NxL + 2;
    ROW_R = ROW_L + 1;

    F(ROW_L) = w1(i); % Left \Gamma_1
    F(ROW_R) = w3(i); % Left \Gamma_3
end

% Left and Right Walls for Lower Legs (\Gamma_i, i = 2,4)
for i = 2:NyL
    ROW_L = (i-1)*Nsml + NxL + 2;
    ROW_R = ROW_L + 1;

    F(ROW_L) = w2(i); % Left \Gamma_2
    F(ROW_R) = w4(i); % Left \Gamma_4
end

% Upper Mid-H Wall \Gamma_5
for j = NxL + 2 : NxM + 1
    ROW = NyL*Nsml + (NyM - NyL) * (Nx+2) + j;
    F(ROW) = w5(j); % Left \Gamma_5
end

% Lower Mid-H Wall \Gamma_6
for j = NxL + 2 : NxM + 1
    ROW = NyL*Nsml + j;
    F(ROW) = w6(j); % Left \Gamma_5
end
% --- End of Boundary Conditions ----------------------------------------



% ---- Fill Inner Cells in Lower Legs of H ------------------------------
for i = 2:NyL-1
    j_ctr = NxL+3;
    
for j = [2 : NxL+1, NxM+2 : Nx+1]
    
    if j < NxM+2
        ROW = Nsml*(i-1)+j;
    else 
        j_ctr = j_ctr + 1;
        ROW = Nsml*(i-1)+j_ctr;
    end
    % Cell height/width
    dx = XU(i,j) - XU(i,j-1);
    dy = YU(i+1,j)-YU(i,j);
    
    % Location of Nodes
    XE = XV(i,j+1);  XP = XV(i,j);  XW = XV(i,j-1);
    YN = YV(i+1,j);  YP = YV(i,j);  YS = YV(i-1,j);
    
    % Velocities
    VP = V0(i,j);
    VE = V0(i,j+1); VN = V0(i+1,j); VW = V0(i,j-1); VS = V0(i-1,j);
   
    % Diffusion Terms
    dw = -1/Re * ( VP - VW )*dy/(XP-XW)/2;
    de = 1/Re * ( VE - VP )*dy/(XE-XP)/2;
    ds = -1/Re * ( VP - VS )*dx/(YP-YS)/2;
    dn = 1/Re * ( VN - VP )*dx/(YN-YP)/2;
    D = dw + de + ds + dn;

    % Pressure
    Pn = P0(i,j-1); Ps = P0(i-1,j-1);
    PR = (Pn-Ps) * dx;

    F( ROW ) = V0(i,j)*dx*dy/dt - 1.5*NLV(i,j) + ...
                              + 0.5*NLV0(i,j) + D - PR + ...
                              + 0.5*(Qv(i,j) + Qvo(i,j))*dy*dx;
  
end  
end


% ---- Fill Inner cells at Lower Interface yv_i, i = NyL, NyL+1 ---
for i = NyL : NyL+1
    j_ctr = NxL+3;
    
for j = [2 : NxL+1, NxM+2 : Nx+1]
    
    if (i == NyL) && (j > NxM+1)
        j_ctr = j_ctr + 1;
        ROW = Nsml*(i-1)+j_ctr;
    else
        ROW = Nsml*(i-1)+j;
    end
    
    % Cell height/width
    dx = XU(i,j) - XU(i,j-1);
    dy = YU(i+1,j)-YU(i,j);
    
    % Location of Nodes
    XE = XV(i,j+1);  XP = XV(i,j);  XW = XV(i,j-1);
    YN = YV(i+1,j);  YP = YV(i,j);  YS = YV(i-1,j);
    
    % Velocities
    VP = V0(i,j);
    VE = V0(i,j+1); VN = V0(i+1,j); VW = V0(i,j-1); VS = V0(i-1,j);
   
    % Diffusion Terms
    dw = -1/Re * ( VP - VW )*dy/(XP-XW)/2;
    de = 1/Re * ( VE - VP )*dy/(XE-XP)/2;
    ds = -1/Re * ( VP - VS )*dx/(YP-YS)/2;
    dn = 1/Re * ( VN - VP )*dx/(YN-YP)/2;
    D = dw + de + ds + dn;

    % Pressure
    Pn = P0(i,j-1); Ps = P0(i-1,j-1);
    PR = (Pn-Ps) * dx;

    F( ROW ) = V0(i,j)*dx*dy/dt - 1.5*NLV(i,j) + ...
                              + 0.5*NLV0(i,j) + D - PR + ...
                              + 0.5*(Qv(i,j) + Qvo(i,j))*dy*dx;
  
end  
end


% ---- Fill inner cells in middle of H for i = NyL +2, ..., NyM
for i = NyL+2 : NyM
for j = 2:Nx+1
    
    ROW = NyL*Nsml + (i - NyL-1)*(Nx+2) + j;
    
    % Cell height/width
    dx = XU(i,j) - XU(i,j-1);
    dy = YU(i+1,j)-YU(i,j);
    
    % Location of Nodes
    XE = XV(i,j+1);  XP = XV(i,j);  XW = XV(i,j-1);
    YN = YV(i+1,j);  YP = YV(i,j);  YS = YV(i-1,j);
    
    % Velocities
    VP = V0(i,j);
    VE = V0(i,j+1); VN = V0(i+1,j); VW = V0(i,j-1); VS = V0(i-1,j);
   
    % Diffusion Terms
    dw = -1/Re * ( VP - VW )*dy/(XP-XW)/2;
    de = 1/Re * ( VE - VP )*dy/(XE-XP)/2;
    ds = -1/Re * ( VP - VS )*dx/(YP-YS)/2;
    dn = 1/Re * ( VN - VP )*dx/(YN-YP)/2;
    D = dw + de + ds + dn;

    % Pressure
    Pn = P0(i,j-1); Ps = P0(i-1,j-1);
    PR = (Pn-Ps) * dx;

    F( ROW ) = V0(i,j)*dx*dy/dt - 1.5*NLV(i,j) + ...
                              + 0.5*NLV0(i,j) + D - PR + ...
                              + 0.5*(Qv(i,j) + Qvo(i,j))*dy*dx;
  
end  
end


% ---- Fill Inner cells at Upper Interface yv_i, i = NyM+1, NyM+2 ----
for i = NyM+1 : NyM+2
    j_ctr = NxL+3;
    
for j = [2 : NxL+1, NxM+2 : Nx+1]
    
    if (i == NyM+2) && (j > NxM+1)
        j_ctr = j_ctr + 1;
        ROW = Nsml*NyL + (NyM+1 - NyL) * (Nx+2) + j_ctr;
    else
        ROW = Nsml*NyL + (i-1 - NyL) * (Nx+2)+j;
    end
    
    % Cell height/width
    dx = XU(i,j) - XU(i,j-1);
    dy = YU(i+1,j)-YU(i,j);
    
    % Location of Nodes
    XE = XV(i,j+1);  XP = XV(i,j);  XW = XV(i,j-1);
    YN = YV(i+1,j);  YP = YV(i,j);  YS = YV(i-1,j);
    
    % Velocities
    VP = V0(i,j);
    VE = V0(i,j+1); VN = V0(i+1,j); VW = V0(i,j-1); VS = V0(i-1,j);
   
    % Diffusion Terms
    dw = -1/Re * ( VP - VW )*dy/(XP-XW)/2;
    de = 1/Re * ( VE - VP )*dy/(XE-XP)/2;
    ds = -1/Re * ( VP - VS )*dx/(YP-YS)/2;
    dn = 1/Re * ( VN - VP )*dx/(YN-YP)/2;
    D = dw + de + ds + dn;

    % Pressure
    Pn = P0(i,j-1); Ps = P0(i-1,j-1);
    PR = (Pn-Ps) * dx;

    F( ROW ) = V0(i,j)*dx*dy/dt - 1.5*NLV(i,j) + ...
                              + 0.5*NLV0(i,j) + D - PR + ...
                              + 0.5*(Qv(i,j) + Qvo(i,j))*dy*dx;
  
end  
end


% ---- Fill Inner Cells in Upper Legs of H ------------------------------
for i = NyM+3 : Ny
    j_ctr = NxL+3;
    
for j = [2 : NxL+1, NxM+2 : Nx+1]
    
    if j < NxM+2
        ROW = Nsml*NyL + (NyM+1 - NyL)*(Nx+2) + (i - NyM-2)*Nsml + j;
    else 
        j_ctr = j_ctr + 1;
        ROW = Nsml*NyL + (NyM+1 - NyL)*(Nx+2) + (i - NyM-2)*Nsml + j_ctr;
    end
    
    % Cell height/width
    dx = XU(i,j) - XU(i,j-1);
    dy = YU(i+1,j)-YU(i,j);
    
    % Location of Nodes
    XE = XV(i,j+1);  XP = XV(i,j);  XW = XV(i,j-1);
    YN = YV(i+1,j);  YP = YV(i,j);  YS = YV(i-1,j);
    
    % Velocities
    VP = V0(i,j);
    VE = V0(i,j+1); VN = V0(i+1,j); VW = V0(i,j-1); VS = V0(i-1,j);
   
    % Diffusion Terms
    dw = -1/Re * ( VP - VW )*dy/(XP-XW)/2;
    de = 1/Re * ( VE - VP )*dy/(XE-XP)/2;
    ds = -1/Re * ( VP - VS )*dx/(YP-YS)/2;
    dn = 1/Re * ( VN - VP )*dx/(YN-YP)/2;
    D = dw + de + ds + dn;

    % Pressure
    Pn = P0(i,j-1); Ps = P0(i-1,j-1);
    PR = (Pn-Ps) * dx;

    F( ROW ) = V0(i,j)*dx*dy/dt - 1.5*NLV(i,j) + ...
                              + 0.5*NLV0(i,j) + D - PR + ...
                              + 0.5*(Qv(i,j) + Qvo(i,j))*dy*dx;
  
end  
end
