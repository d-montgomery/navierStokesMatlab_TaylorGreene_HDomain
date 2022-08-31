function [PHI] = reshapePHI_Ghst(NxL,NxM,Nx,NyL,NyM,Ny,TMPP)

% --- Reshape PHI ---
PHI = zeros(Ny,Nx);

% Size of sml sub matrices = Nsub x Nsub
Nsml = (NxL) + (Nx - NxM); 

% Lower Legs
j = [1:NxL, NxM+1:Nx];
for i = 1:NyL
    j_ctr = (1 : Nsml) + (i-1)*Nsml;
    PHI(i,j) = TMPP( j_ctr ); 
end

% Middle of H
j = 1:Nx;
for i = NyL+1 : NyM
    j_ctr = NyL*Nsml + j + (i-NyL-1)*(Nx);
    PHI(i,j) = TMPP( j_ctr ); 
end

% Upper Legs
j = [1:NxL, NxM+1:Nx];
for i = NyM+1 : Ny
    j_ctr = NyL*Nsml + (NyM - NyL)*(Nx) + (1 : Nsml) + (i-NyM-1)*Nsml;
    PHI(i,j) = TMPP( j_ctr ); 
end
