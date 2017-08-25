function X= mhde(phat)

num = phat(1)^(1/2) - phat(3)^(1/2);
denom = 2 * ( 1 + phat(2) - 2*phat(1)^(1/2)*phat(3)^(1/2) )^(1/2);
if num==0
   X = 1/2;

elseif (1/2+num/denom) < 0
    X = 0;
  
elseif (1/2+num/denom) > 1 
    X = 1;

else
    X = 1/2+num/denom;
end
return 

