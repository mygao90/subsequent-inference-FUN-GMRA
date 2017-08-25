function [dsuper, dsub, dGMRA] = exactdistr(m,nrand) 
global Nmat theta0 n K;
dsuper = zeros(size(Nmat,1),1);
dsub = zeros(size(Nmat,1),1);
dGMRA =[];
pthetamat = auxiliary(m,nrand) ;
for i=1:size(Nmat,1)
   phat = Nmat(i,:)/sum(Nmat(i,:));
   dsuper(i) = infodistSimplex(phat,HWpvec(theta0));
   dsub(i) = infodistHW( mhde(phat) , theta0 );
   Matr = [pthetamat;phat;HWpvec(theta0)].';
   %dGMRA(i) = doLLE(Matr,K,1);
   dGMRA(i,:)= doGMRA(sqrt(Matr),1);
end
return 
