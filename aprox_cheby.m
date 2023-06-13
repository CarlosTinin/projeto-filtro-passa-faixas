function [n,d,rh,rv]=aprox_cheby(ordem,Amax,wp)

% Aproximacao de Chebyshev
% Retorna o numerador e denominador da funcao T(s) de Chebyshev
% Sintaxe: [n,d]=aprox_cheby(ordem,Amax,wp)
% Amax = maxima atenuacao na banda de passagem
% wp = frequencia limite da banda de passagem


E=sqrt(10^(0.1*Amax)-1);	% maxima distorcao na banda de passagem
%iE=1/E;
%iordem=1/ordem;
argh=(1/ordem)*asinh(1/E);
rh=sinh(argh);
rv=cosh(argh);
v2=ordem/2;
iparte=floor(v2);
resto=v2-iparte;
if resto==0.0
	limite=ordem/2;
	ispar=1;
else
	limite=(ordem+1)/2;
	ispar=0;
end

for k=1:limite
	sigma(k)=sin((((2*k)-1)*pi)/(2*ordem))*(sinh(argh));
	omega(k)=cos((((2*k)-1)*pi)/(2*ordem))*(cosh(argh));
end
%ip=-sigma+j*omega;
for k=1:limite
  if floor(omega(k)/sigma(k))> 0.0
        in(k)=-sigma(k)-j*omega(k);
	ip(k)=-sigma(k)+j*omega(k);
  else
	ip(k)=-sigma(k);
  end
end

rz=[in ip]';
pc=real(poly(rz));
a0=pc(length(pc));
if ispar==1
	K=a0*(10^(-0.05*Amax));
else
	K=a0;
end
num=zeros(size(pc));
num(length(pc))=K;
den=pc;
% desnormalizacao para wp => s_=s/wp
%E=sqrt(10^(0.1*Amax)-1);	% maxima distorcao na banda de passagem
%fdes=E^(1/ordem)/wp;
fdes=1/wp;		%fator de desnormalização
vpot=[length(den)-1:-1:0];
Hd=zeros(size(den));
for x=1:length(vpot)
	Hd(x)=den(x)*(fdes^vpot(x));
end
n=num/Hd(1);
d=Hd/Hd(1);
return















