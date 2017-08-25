 % enumerate all possible values for the observation Nvec ~ Multinomial(n,[3],p)
  % assumes global value for "n"
function Nmat = makeNmat(n)
nrow = (n+1)*(n+2)/2;
Nmat = zeros(nrow,3);
i=0;
for N1=0:n
    for N2=0:(n-N1)
      N3 = n-N1-N2 ;
      Nvec = [N1 N2 N3];
      i = i+1 ;
      Nmat(i,:) = Nvec;
    end
end

return
