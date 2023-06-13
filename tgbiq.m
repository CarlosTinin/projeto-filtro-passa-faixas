function [p] = tgbiq(ni,di,w)

%TGBIQ	Calcula o atraso de grupo para secoes biquadraticas
%	Funcao utilizada por mag_fase
%	Sintaxe: [tg]=tgbiq_v2(num,den,w)
%	      num = numerador da funcao de transferencia
%       den = denominador da funcao de transferencia
%       w   = vetor de frequencias em rad/s
%	tg = - d0(w)/dw
%	Copyright (c) 2015 by LaPeSC/Delmar B. Carvalho.
% Versao: 2
% Data: 12/06/2022

[l,c]=size(ni);
if c > 3
 [n,d]=pol2biq(ni,di);
else
 n=ni;
 d=di;
end
[l,c]=size(n);
ns=l; %numero de secoes
if ns==1
   w2=w.^2;
   f1=n(1)*w2; %aw^2
   f2=n(3)-f1; %c-aw^2
   f3=n(2)^2.*w2; %b^2w^2
   fn=((n(2).*f2)+(2*n(2).*f1))./(f3+f2.^2);
   g1=d(1)*w2;
   g2=d(3)-g1;
   g3=d(2)^2.*w2;
   gd=((d(2).*g2)+(2*d(2).*g1))./(g3+g2.^2);
   p= -fn+gd;
else
   p=w.*0;
   for i=1:ns
       p=p+ tgbiq(n(i,:),d(i,:),w);
   end
end
return

