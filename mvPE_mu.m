function [Out_MVPE,ppp]=mvPE_mu(x,m,t)

CH=size(x,1);
ly=size(x,2); % Length of each signal

for i=1:CH   
 y=x(i,:);  
permlist = perms(1:m);
c(1:length(permlist))=0;
    
 for j=1:ly-t*(m-1)
     [a,iv]=sort(y(j:t:j+t*(m-1)));
     for jj=1:length(permlist)
         if (abs(permlist(jj,:)-iv))==0
             c(jj) = c(jj) + var(y(j:t:j+t*(m-1)),1) ;
         end
     end
 end

hist = c;
 
p(i,:) = c/(sum(c)*CH);
end

pp=sum(p,1);

ppp=pp(pp~=0);  

Out_MVPE=-sum(ppp .* log(ppp));
Out_MVPE=Out_MVPE/log(factorial(m));

