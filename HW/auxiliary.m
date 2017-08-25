function pthetamat = auxiliary(m,nrand)
thetai = unifrnd(0,1,m,1);
pthetamat = zeros(m,3);
for i=1:m
    pthetamat(i,:) = HWpvec(thetai(i));
    pthetamat(i,:) = histc( randsample(3,nrand,true,pthetamat(i,:)), 1:1:3)/nrand;
end
return