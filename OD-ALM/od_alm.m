function [V,Sigma,angle,rechgvec,iter] = od_alm(A,U,memory)
% od_alm compute orthogonal tensor decomposition by the augmented Lagrangian method
%
% od_alm(A,U,memory) computes the orthogonal low rank approximation of A.
% The input A is a double tensor. The initialization U is a cell.
% The number of terms is determined by U automatically.
% The subproblems are solved by L-BFGS, where memory is the level of memory
% 
% [Sigma;V] is the output. Sigma is a vector, and V is a cell
% angle is the vector of the angle of each iteration
% rechgvec is the vector of relative change 
% iter is the vector of iteration numbers of inner loops

%% set parameters
Tol = 1e-4;    %% angle tolerance
maxiter = 24;  %% maximal iteration numbers
incCoe = 10;   %% increasing coefficient of penlty parameters

%%
dims = size(A);
N = ndims(A);
R = size(U{1},2);
tA = tensor(A);
Up = U;
x_pre = cell2vec(Up,N,R,dims);
Angle = zeros(R,R,N);
for n = 1:N
    Angle(:,:,n) = Up{n}'*Up{n};
end
Ipro = Angle;

angle = [];
rechgvec = [];
iter = [];
c = 1;             %% penlty parameter
C = c*ones(R);
lambda = 0;        %%  Lagrange multipliers
Lambda = lambda*ones(R);
kk = 1;
while (1)
    
    [Up,Lambda,Ipro,num] = alm_LBFGS(A,Up,C,Lambda,Ipro,memory);
    x = cell2vec(Up,N,R,dims);
    rechange = norm(x-x_pre)/norm(x_pre);
    x_pre = x;
    for n = 1:N
        normvec = diag(Ipro(:,:,n));
        normvec = sqrt(normvec);
        normvec = 1./normvec;
        normmatrix = normvec*normvec';
        Angle(:,:,n) = normmatrix.*Ipro(:,:,n);
        Angle(:,:,n) = Ipro(:,:,n);
    end
    Angle = abs(Angle);
    minmatrix = min(Angle,[],3);
    minmatrix(logical(eye(R)))= 0;
    maxinner = max(max(minmatrix));
    angle = [angle maxinner];
    rechgvec = [rechgvec rechange];
    iter = [iter num];
    if maxinner < Tol ||  kk>maxiter
        break;
    end
    kk = kk +1;
    C = incCoe*C;
end

V = Torth(Up);
U_mttkrp = mttkrp(tA,V,N);
Sigma = sum(V{N}.* U_mttkrp);
