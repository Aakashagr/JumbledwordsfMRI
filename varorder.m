% n = length of string
function neworder = varorder(n)
cnt = 0; neworder = [];
for i = 1:n
    for j = n:-1:i
        cnt = cnt + 1;
        if j == n; neworder(cnt) = i;
        else, neworder(cnt) = neworder(cnt-1)+j+1; end
    end
end

ordermax = max(neworder);

for i = 1:n-1
    for j = n:-1:i+1
        cnt = cnt + 1;        
        if j == n; neworder(cnt) = ordermax + i;
        else, neworder(cnt) = neworder(cnt-1)+j; end
    end
end


