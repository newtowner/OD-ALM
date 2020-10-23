function U = vec2cell(x,N,R,dims)

U = cell(N,1);
num = 0;
for n = 1:N
    inc = dims(n)*R;
    num2 = num + inc;
    tmp = x(num+1:num2);
    U{n} = reshape(tmp,[],R);
    num = num2;
end