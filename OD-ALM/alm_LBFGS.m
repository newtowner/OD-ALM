function [Up,Lambdap,Iprop,iter] = alm_LBFGS(A,U,C,Lambda,Ipro,m)

%m = 20;
%% basic value
dims = size(A);
N = ndims(A);
R = size(U{1},2);
Up = U;
Iprop = Ipro;
tA = tensor(A);
normA = norm(tA);

%% compute penlty matrix 
M_Ipro = prod(Ipro,3);
delta_s = diag(M_Ipro);        % norm squared of each rank-one component
indelta = 1./delta_s;
M_coe = indelta*indelta';
C = C.*M_coe;
C(logical(eye(R)))= 0;

%% preprocess input, such that each vector of the same term has the same norm
meandelta = delta_s.^(1/N);
for n = 1:N
    norm_mode = diag(Ipro(:,:,n));       % norm squared of each vector of mode-n
    ratio = sqrt(meandelta./norm_mode);
    Up{n} = bsxfun(@times, Up{n}, ratio'); % Matlab 2016
end


%% Set algorithm parameters 
maxiters = 500;
gtol = 1e-4;
rtol = 1e-8;
%% Initialize
[f,g_pre] = fun_alm(Up,tA,N,R,normA,Lambda,C,dims);
num_var = sum(dims)*R;
x_pre = cell2vec(Up,N,R,dims);
S = zeros(num_var,m);
Y = zeros(num_var,m);
rho = zeros(1,m);
alpha = rho;
%% main loop
for k = 1:maxiters
    if k == 1
        stp = 1.0;
        dir = -g_pre;
    else
        dir = -r;
    end
 
    [Up,fnew,g,stp,~,~] = cvsrch_alm(Up,tA,N,R,normA,Lambda,C,dims,f,g_pre,stp,dir);
    Gra = norm(g)/num_var;
    x = cell2vec(Up,N,R,dims);
    rechange = norm(x-x_pre)/norm(x_pre);
    if (Gra < gtol || rechange < rtol || k == maxiters) && k > 1
        iter = k;
        %Gra;
        for n = 1:N
            Iprop(:,:,n) = Up{n}'*Up{n};
        end
        innerU = prod(Iprop,3);
        Lambdap = innerU.*C;
        Lambdap = Lambda + Lambdap;
        break;
    end
    f = fnew;   
    s = x - x_pre;
    S = [s S(:,1:m-1)];
    x_pre = x;
    y = g - g_pre;
    Y = [y Y(:,1:m-1)];
    g_pre = g;    
    sy = y'*s;
    yy = y'*y;
    gamma = sy/yy;
    rho = [1/sy rho(:,1:m-1)];
    
    r = g;
    for n = 1:min(m,k)
        s = S(:,n);
        alpha(n) = rho(n)*(s'*r);
        r = r - alpha(n)*Y(:,n);
    end
    r = gamma * r;
    for n = min(m,k):-1:1
        y = Y(:,n);
        beta = rho(n)*(y'*r);
        r = r + (alpha(n)-beta)*S(:,n);
    end
 
end
    