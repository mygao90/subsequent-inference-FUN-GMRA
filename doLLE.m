function res = doLLE(X,k,d)
Y = lle(X,k,d);
N =size(Y,2);
res = norm(Y(:,N-1)-Y(:,N));