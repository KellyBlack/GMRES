%\textasciitilde 

N = 20;
A = zeros(N);
for i=1:N,
 A(i,i) = sqrt((i+1)/i);
end

for i=2:N,
 A(i,i-1) = sqrt((i-1)/i);
end


N = 10;
m = 5;
B = zeros(N);
r = 2+m;
for i=1:N,
    B(i,i) = sqrt(r);
    r = (r*(2+m)-1)/r;
end

for i=2:N,
 B(i,i-1) = -1.0/B(i-1,i-1);
end

