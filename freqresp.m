function r=freqresp(p1,p2,p3)

%--------------------------------------------
%FREQRESP	funcao que calcula a resposta em frequencia
%		Sintaxe: mag_complexa= freqresp(num,den,w)
%		Parametros:
%		p1 =  numerador 
%		p2 =  denominador 
%		p3 = vetor de frequencias em radianos
%
%
%	Copyright (c) 2015 by LaPeSC/Delmar B. Carvalho.



s=j.*p3;
num=polyval(p1,s);
den=polyval(p2,s);
r=num./den;

return
