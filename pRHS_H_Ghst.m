function[Fp] = pRHS_H_Ghst(N,U,V,dt,XU,YV)

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% --- Matrix Sizes ---
% Ap is NNXY x NNXY matrix without augmentation
NNXY = (Nx)*(Ny) - (NxM - NxL) * ( Ny - NyM + NyL);
% Size of sml sub matrices = Nsml x Nsml
Nsml = (NxL) + (Nx - NxM); 

Fp = zeros(NNXY,1);

% ---- Fill Cells in Lower Legs of H ------------------------------
for i = 1:NyL
    j_ctr = NxL;
for j = [1 : NxL, NxM+1 : Nx]
    if j < NxL+1
        ROW_ind = Nsml*(i-1)+j;
    else 
        j_ctr = j_ctr + 1;
        ROW_ind = Nsml*(i-1)+j_ctr;
    end
    
    dx = XU(i+1,j+1) - XU(i,j);  
    dy = YV(i+1,j+1) - YV(i,j+1);
   
    Fp( ROW_ind ) = 1/dt*( ( U(i+1,j+1) - U(i+1,j) )/dx ...
                          +  ( V(i+1,j+1) - V(i,j+1) )/dy );
    
end  
end

% ---- Fill Cells in Injury Channel of H ------------------------------
for i = NyL+1 : NyM   
for j = 1:Nx
    ROW_ind = Nsml*(NyL) + (i - NyL -1) * Nx + j;
    
    dx = XU(i+1,j+1) - XU(i,j);  
    dy = YV(i+1,j+1) - YV(i,j+1);
   
    Fp( ROW_ind ) = 1/dt*( ( U(i+1,j+1) - U(i+1,j) )/dx ...
                          +  ( V(i+1,j+1) - V(i,j+1) )/dy );
    
end  
end


% ---- Fill Cells in Upper Legs of H ------------------------------
for i = NyM+1:Ny
    j_ctr = NxL;
for j = [1 : NxL, NxM+1 : Nx]
    if j < NxL+1
        ROW_ind = Nsml*NyL + Nx*(NyM-NyL) + (i-NyM-1)*Nsml + j;
    else 
        j_ctr = j_ctr + 1;
        ROW_ind = Nsml*NyL + Nx*(NyM-NyL) + (i-NyM-1)*Nsml + j_ctr;
    end
    
    dx = XU(i+1,j+1) - XU(i,j);  
    dy = YV(i+1,j+1) - YV(i,j+1);
   
    Fp( ROW_ind ) = 1/dt*( ( U(i+1,j+1) - U(i+1,j) )/dx ...
                          +  ( V(i+1,j+1) - V(i,j+1) )/dy );
    
end  
end

% % Augmentation
% Fp(NNXY+1) = 0;





