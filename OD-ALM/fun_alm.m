function [f,g] = fun_alm(U,tA,N,R,normA,Lambda,C,dims)
% Compute the augmented Lagrangian function value and the gradient
% at point x
% tA is a tensor

%% compute function value f
ktnew = ktensor(U);
U_mttkrp = mttkrp(tA,U,N);
iprod = sum(sum(U{N}.* U_mttkrp));
Left = normA^2 + norm(ktnew)^2 - 2 * iprod;
Left = Left/2;

UtU = zeros(R,R,N);
for n = 1:N
    UtU(:,:,n) = U{n}'*U{n};
end
innerpro = prod(UtU,3);
Lag = Lambda.*innerpro;
Middle = sum(sum(Lag))/2;

innerpro_s = innerpro.^2;
weiinner = C.*innerpro_s;  
Right = sum(sum(weiinner))/4;
f = Left + Middle + Right;

%% compute gradient g
sumlength = sum(dims);
g = zeros(sumlength*R,1);
num = 0;
for n = 1:N 
    V = mttkrp(tA,U,n);
    Gram = prod(UtU(:,:,[1:n-1 n+1:N]),3);     
    VnTVn = UtU(:,:,n);
    M = Gram.^2;
    M = VnTVn.*M;
    M = C.*M;  
    M = M + Gram + Gram.*Lambda;   
    Gra = U{n}*M;
    Gra = Gra - V;   
 
    inc = dims(n)*R;
    num2 = num + inc;
    g(num+1:num2) = Gra(:);
    num = num2;
end

