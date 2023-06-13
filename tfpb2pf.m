function [np,dp]=tfpb2pf_v2(ni,di,Amax,B,wo,op)

% Transformacao em frequencia PB -> PF (versao 25/05/2023)
% Requer o pacote signal e control
% Retorna o numerador e denominador da funcao T(s)
% Sintaxe: [n,d]=tfpb2rf(n,d,B,wo)
% n = coeficientes do polinomio do numerador
% d = coeficientes do polinomio do denominador
% B = largura de banda (w2-w1)		[ w3\  |w1  w2|  /w4]
% wo = frequencia central
% op = ´B' para Butterworth
%      'C' para Chebyshev
%


rd=roots(di);
ordem=length(rd);
if ((op == 'C') && (rem(ordem,2)==1))
  rz=roots([1 B -wo^2]);
  w1=rz(2);
  w2=wo^2/w1;
  [n,d]=cheby1(ordem,Amax,[w1 w2],'s');
  np=[zeros(length(d)-length(n),1)' n];
  dp=d;
else
  [nbq,dbq]=pol2biq(ni,di); %converte polinomio em biquadratica
  [Nbiq,c]=size(nbq);
  n=zeros(Nbiq,c+2);
  d=zeros(Nbiq,c+2);
  x=1/B;
  y=x*wo^2;
  for ii=1:Nbiq
    if dbq(ii,1) == 0 % nao apresenta fator quadratico = 1a ordem
      %n(ii,1) = 0
      n(ii,4) = nbq(ii,3)/x;
      %n(ii,3) = 0;
      d(ii,3) = 1;
      d(ii,4) = dbq(ii,2)/x;
      d(ii,5) = y/x;
    else
      n(ii,3) = nbq(ii,3)/x^2;
      d(ii,1) = 1;
      d(ii,2) = dbq(ii,2)/x;
      d(ii,3) = (dbq(ii,3)/x^2)+(2*y/x);
      d(ii,4) = (dbq(ii,2)*y)/x^2;
      d(ii,5) = y^2/x^2;
    endif
  endfor
  [nx,dx]=biq2pol(n,d);
  mo=mag_fase(ni,di,0); %  avalia a magnitude na origem
  mt=mag_fase(nx,dx,wo); % avalia a magnitude em w0

  if mt/mo == 1
    np=nx;dp=dx;
  else
    np=nx.*(mo/mt);dp=dx;
  endif
endif

return












