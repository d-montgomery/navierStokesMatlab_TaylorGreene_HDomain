function [U,V,P] = reshapeFinal(N,TMPU,TMPV,TMPP)

% Reshape final U, V, P for convergence analysis. Note that we remove the
% boundaries to avoid division by zero.

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% --- Reshape TMPU to U ---

% Lower Legs
j = [2:NxL, NxM+2:Nx];
Nsml = length(j);
for i = 3:NyL+2
    i_ctr = (1 : Nsml) + (i-3)*Nsml;
    U(i_ctr) = TMPU( i, j )'; 
end

% Middle of H
j = 2:Nx;
Nbig = length(j);
for i = NyL+3 : NyM
    i_ctr = (NyL-1)*Nsml + (1:Nbig) + (i-NyL-3)*Nbig;
   U(i_ctr) = TMPU(i, j)'; 
end

% Upper Legs
j = [2:NxL, NxM+2:Nx];
for i = NyM+1 : Ny+1
    i_ctr = (NyL-1)*Nsml + (NyM+1 - NyL-3)*Nbig + (1 : Nsml) + (i-NyM-1)*Nsml;
   U(i_ctr) = TMPU(i,j)'; 
end




% --- Reshape TMPV to V ---


% Lower Legs
j = [3:NxL, NxM+3:Nx];
Nsml = length(j);
for i = 2:NyL+1
    i_ctr = (1 : Nsml) + (i-2)*Nsml;
    V(i_ctr) = TMPV( i,j )'; 
end

% Middle of H
j = 3:Nx;
Nbig = length(j);
for i = NyL+2 : NyM
    i_ctr = NyL*Nsml + (1:Nbig) + (i-NyL-2)*Nbig;
    V(i_ctr) = TMPV( i,j )'; 
end

% Upper Legs
j = [3:NxL, NxM+3:Nx];
for i = NyM+1 : Ny
    i_ctr = NyL*Nsml + (NyM+1 - NyL-2)*Nbig + (1 : Nsml) + (i-NyM-1)*Nsml;
    V(i_ctr) = TMPV( i,j )'; 
end



% --- Reshape TMPP to P ---
% Lower Legs
j = [2:NxL-1, NxM+2:Nx-1];
Nsml = length(j); 
for i = 3:NyL+2
    i_ctr = (1 : Nsml) + (i-3)*Nsml;
    P(i_ctr) = TMPP( i,j )'; 
end

% Middle of H
j = 3:Nx;
Nbig = length(j);
for i = NyL+3 : NyM
    i_ctr = NyL*Nsml + (1:Nbig) + (i-NyL-3)*Nbig;
    P(i_ctr) = TMPP( i,j )'; 
end

% Upper Legs
j = [2:NxL-1, NxM+2:Nx-1];
Nsml = length(j); 
for i = NyM+1 : Ny
    i_ctr = NyL*Nsml + (NyM+1 - NyL-3)*Nbig + (1 : Nsml) + (i-NyM-1)*Nsml;
    P(i_ctr) = TMPP( i,j )'; 
end