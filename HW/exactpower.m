function power = exactpower(tvec)
global Nmat n theta0 alpha thetaAseq
pvec0 = zeros(size(Nmat,1),1) ; 
for i= 1:size(Nmat,1) 
    pvec0(i) = mnpdf(Nmat(i,:),HWpvec(theta0));
end
[~,otv] = sort(tvec);
utv = unique(tvec(otv));
pmf = zeros(size(utv,1),1) ; 
cdf = zeros(size(utv,1),1) ; 
for i =1:size(utv,1)
    tv = utv(i);
    pmf(i) = sum(pvec0(tvec==tv));
    cdf(i) = sum(pmf(1:i));
end 
  cvindex = find(cdf>(1-alpha));
  cvi= cvindex(1);
  cv = utv(cvi);
  prej = sum(pvec0(tvec>cv));
 
  pcv = sum(pvec0(tvec==cv));
  rprob = alpha-prej;
  power = 1:size(thetaAseq,1);
for A= 1:size(thetaAseq,1) 
    pvecA= 1:size(Nmat,1) ; 
    for i = 1:size(Nmat,1) 
        pvecA(i) = mnpdf(Nmat(i,:),HWpvec(thetaAseq(A)));
    power(A) = sum(pvecA(tvec>cv)) + (rprob/pcv)*sum(pvecA(tvec==cv));
    end  
end
return
