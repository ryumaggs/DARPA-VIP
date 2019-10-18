function X = Generate_Meyer(n,w1,w2,T1,T2)

%% Meyer staircase
w     = @(t) w1+(w2-w1)*(t-1)/(n-1);
T     = @(t) T1+(T2-T1)*(t-1)/(n-2);

Supp  = [1 1+w1];

for k = 2 : n
    % start Supp(k-1,1)
    % translate T(k-1)
    % width w(k)
    Supp = [Supp ; Supp(k-1,1)+T(k-1) Supp(k-1,1)+T(k-1)+w(k)];
end

Supp = round(Supp);
D    = max(max(Supp));
X    = zeros(D,n);

for k = 1 : n
    X(Supp(k,1):Supp(k,2),k) = 1;
end

