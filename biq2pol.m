function  [p1,p2] = biq2pol(a,b)

%BIQ2POL Converte o conjunto de funções biquadráticas em um polinômio
%	  Sintaxe: [num,den] = biq2pol(numbiq,denbiq)
%         numbiq = numerador da função de transferência no formato biquadrática
%         denbiq = denominador da função de transferência no formato biquadrática
%	  num = numerador da função de transferência no formato polinomial
%	  den = denominador da função de transferência no formato polinomial
%
%	Copyright (c) 2015 by LaPeSC/Delmar B. Carvalho.


h=a(1,:);
k=b(1,:);
[l,c]=size(a);
if l>1
   for i=2:l
       h= conv(h,a(i,:)); 
       k= conv(k,b(i,:));
   end
end
[i,j]=find(k);  %encontra a posicao do primeiro elemento nao zero
p1=h(j(1):length(h));
p2=k(j(1):length(k));

return
