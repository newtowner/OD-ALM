function x = cell2vec(U,N,R,dims)

Num = sum(dims);
x = zeros(Num*R,1);
num = 0;
for n = 1:N
    inc = dims(n)*R;
    num2 = num + inc;
    x(num+1:num2) = U{n}(:);
    num = num2;
end