global n Nmat theta0 alpha thetaAseq d K;
n = 10;
theta0= 0.3;
alpha =0.05;
d=1;
thetaAseq = transpose(0:0.01:1);
m=25000;
K=4;
nrand =20; 
B=100;

Nmat = makeNmat(n);
[dsuper, dsub, dGMRA] = exactdistr(m,nrand);
powersuper = exactpower(dsuper); 
powersub = exactpower(dsub);
maxdif =[];
J = size(dGMRA,2);
powerGMRA = []; 


% for i=1:B 
%     [dsuper, dsub, dGMRA] = exactdistr(m,nrand);
%     for j =1:J
%         powerGMRAsum(:,j) = powerGMRAsum(:,j) + exactpower(dGMRA(:,j));
%         %maxdif(j)= max(abs(powersub-powerGMRA));
%     end
% end
for k =1:J
    powerGMRA(:,k) = exactpower(dGMRA(:,k));
   % maxdif(k)= max(abs(powersub-powerGMRA));
end


plot(thetaAseq,powerGMRA(:,5),'b--o', thetaAseq, powersub,'r')
print('twolines','-dpng')

plot(thetaAseq,powerGMRA(:,5).'-powersub)
print('diff','-dpng')


