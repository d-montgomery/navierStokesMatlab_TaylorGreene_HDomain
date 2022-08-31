function [PHI] = reshapePHI(NxL,NxM,Nx,NyL,NyM,Ny,TMPP)

% --- Reshape PHI ---
PHI = zeros(Ny+2,Nx+2);

% Size of sml sub matrices = Nsub x Nsub
Nsml = (NxL + 2) + (Nx + 2 - NxM); 

% Lower Legs
j = [1:NxL+2, NxM+1:Nx+2];
for i = 1:NyL
    j_ctr = (1 : Nsml) + (i-1)*Nsml;
    PHI(i,j) = TMPP( j_ctr ); 
end

% Middle of H
j = 1:Nx+2;
for i = NyL+1 : NyM+2
    j_ctr = NyL*Nsml + j + (i-NyL-1)*(Nx+2);
    PHI(i,j) = TMPP( j_ctr ); 
end

% Upper Legs
j = [1:NxL+2, NxM+1:Nx+2];
for i = NyM+3 : Ny+2
    j_ctr = NyL*Nsml + (NyM+2 - NyL)*(Nx+2) + (1 : Nsml) + (i-NyM-3)*Nsml;
    PHI(i,j) = TMPP( j_ctr ); 
end
