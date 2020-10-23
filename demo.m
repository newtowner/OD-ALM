%% data initialization
N = 4;
dims = [20 16 10 32];

%%Hilbert tensor
dT = zeros(dims);
for i1 = 1:dims(1)
    for i2 = 1:dims(2)
        for i3 = 1:dims(3)
            for i4 = 1:dims(4)
                dT(i1,i2,i3,i4) = 1/(i1+i2+i3+i4-3);
            end
        end
    end
end
T = tensor(dT);
normT = norm(T);


%% initialization 
OR = 5;

multirank = OR*ones(N,1);
D = hosvd(T,1,'ranks',multirank);
U1 = D.U;


tic;
[D,~,out] = cp_als(T,OR,'init',U1,'printitn',0,'tol',1e-6);
toc;
Lambda = D.lambda;
U2 = D.U;
% absorb lambda into the last dimension
U2{end} = U2{end}*diag(Lambda);
rerr = 1-out.fit;


%%


tic;
[Up,lambda,angle,rechgvec,iter] = od_alm(dT,U2,20);
toc;
recover = ktensor(lambda',Up);
dr = double(recover);
diff = dr - dT;
err = norm(diff(:))/normT;



