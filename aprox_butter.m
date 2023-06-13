function [n,d]=aprox_butter(n,Amax,wp)

% Aproximacao de Butterworth
% Retorna o numerador e denominador da funcao T(s) de Butterworth
% Sintaxe: [n,d]=aprox_butter(ordem,Amax,wp)
% Amax = maxima atenuacao na banda de passagem
% wp = frequencia limite da banda de passagem

Hr(1)=1;
Hr(n+1)=1;

for k=2:n
	x=k-1;
	Hr(k)=(cos(((x-1)*pi)/(2*n))/sin((x*pi)/(2*n)))*Hr(x);
end
Hn=Hr(n+1:-1:1);		% funcao atenuacao normalizada
E=sqrt(10^(0.1*Amax)-1);	% maxima distorcao na banda de passagem
fdes=E^(1/n)/wp;		%fator de desnormalização
vpot=[n:-1:0];
Hd=zeros(size(Hn));
for x=1:length(vpot)
	Hd(x)=Hn(x)*(fdes^vpot(x));
end
H=Hd/Hd(1);
d=H;
nn=zeros(size(d));
nn(length(d))=d(length(d));
n=nn;

return















