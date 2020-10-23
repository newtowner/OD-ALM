function V = Torth(U)

N = length(U);
R = size(U{1},2);
V = U;

for n = 1:N
    vnorm = sqrt(sum(U{n}.*U{n}));
    inorm = 1./vnorm;
    V{n} = bsxfun(@times, U{n}, inorm); %2016
end

for r = 2:R
    P = zeros(N,r-1);
    for n = 1:N
        M = V{n}(:,1:r-1);
        v = V{n}(:,r);
        P(n,:) = abs(v'*M);
    end
    Ind = zeros(N,r-1);
    for s = 1:r-1
        [~,m] = min(P(:,s));
        Ind(m,s) = 1;
    end
    for n = 1:N
        row = Ind(n,:);
        num = sum(row);
        if num > 0
            B = V{n}(:,row>0);
            v = V{n}(:,r);
            A = B'*B;
            b = B'*v;
            x = A\b;
            v = v - B*x;
            eta = norm(v);
            V{n}(:,r) = v/eta;
        end
    end
end
    
    
    
    