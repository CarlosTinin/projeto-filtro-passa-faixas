function [p1,p2,p3,p4] = mag_fase(a,b,c,d)

%MAG_FASE Calcula os parametros referente a uma funcao de transf
%	  cujos coeficientes estao nos vetores(matrizes) num e den
%         em ordem decrescente de potencia.
%         Sintaxe: [mag,fase,tg,tp] = mag_fase_v2(num,den,w,op)
%         num = numerador da funcao de transferencia
%         den = denominador da funcao de transferencia
%         w   = vetor de frequencias em rad/s
%
%         op='m' retorna apenas o vetor magnitude
%	        op='f' retorna apenas o vetor fase
%         op='g' retorna apenas o vetor atraso de grupo
%         op='p' retorna apenas o vetor atraso de fase
%
%	  Se acrescentado o sub-indice 'i' a qq opcao retorna valores
%	  para cada secao. Ex: magnitude individual = 'mi'.
%
%	  Se omitido o parametro op retorna todos os parametros/secao:
%   [magnitude, fase, atraso de grupo, atraso de fase]=mag_fase_v2(n,d,w);
%
%	Copyright (c) 2015 by LaPeSC/Delmar B. Carvalho.
% Versao: 2.0
% Data: 12/06/2022

if nargin < 3
	disp(' ')
        disp('Erro...parametros insuficientes')
        disp(' ')
 %       break;
        return;
end
if nargin==3
   d='t';
end
[ln,cn]=size(a);
[ld,cd]=size(b);
if ( cn ~= cd)
	disp(' ')
	disp('Erro: os vetores num,den devem ter mesma dimensao')
	disp(' ')
%   break;
	return;
end
if d=='gi'
  for i=1:ln
      tg(:,i)   = tgbiq(a(i,:),b(i,:),c);
  end
  p1=tg;
elseif d=='tg'| d=='g'
  tg   = tgbiq(a,b,c);
  p1=tg;
elseif d=='pi'
  for i=1:ln
      h(:,i)=freqresp(a(i,:),b(i,:),c);
      fase(:,i) = unwrap(angle(h(:,i)));
      tp(:,i)   = - fase(:,i)./c;
  end
  p1=abs(tp);
elseif d=='tp'| d=='p'
   fase=mag_fase(a,b,c,'ft');
   %if sum((fase>0))
   %   fase=fase-(ceil(abs(min(fase))));
   %end
   tp=-fase./c;
   %tp=c.*0;
   %for i=1:ln
   %    h(:,i)=freqs(a(i,:),b(i,:),c);
   %   fase(:,i) = unwrap(angle(h(:,i)));
   %   tp   =tp+ (- fase(:,i)./c);
   %end

  p1=tp;
elseif d=='mi'
  for i=1:ln
      h(:,i)=freqresp(a(i,:),b(i,:),c);
      m(:,i)=abs(h(:,i));
  end
  p1=m;
elseif d=='mt'| d=='m'
  if ln > 1
      [npol,dpol]=biq2pol(a,b);
  else
        npol=a;dpol=b;
  end
  h=freqresp(npol,dpol,c);
  p1=abs(h);
elseif d=='fi'
  for i=1:ln
      h(:,i)=freqresp(a(i,:),b(i,:),c);
      fase(:,i) = unwrap(angle(h(:,i)));
  end
  p1=fase;
elseif d=='ft'| d=='f'
  if ln > 1
        [npol,dpol]=biq2pol(a,b);
  else
        npol=a;dpol=b;
  end
  h=freqresp(npol,dpol,c);  %obtem a resp. em freq.
  p1=unwrap(angle(h));      %obtem a fase
elseif d=='mf'             %calcula mag e fase total da biq.
  if ln > 1
        [npol,dpol]=biq2pol(a,b);
  else
        npol=a;dpol=b;
  end
  h=freqresp(npol,dpol,c);  %obtem a resp. em freq.
  p1=abs(h);
  p2=unwrap(angle(h));   %obtem a fase

elseif d=='ef'
   [mag,fase]=mag_fase(a,b,c,'mf');
   dw=(c(length(c))-c(1))/(length(c)-1);
   E=(mag.^2).*(dw/pi);
   Et=sum(E);
   tg=mag_fase(a,b,c,'tg');
   tgp=(E'/Et)*tg;
   tp=mag_fase(a,b,c,'tp');
   tpp=(E'/Et)*tp;
   tx=erimp(a,b,'X');          %obtem o instante de maxima amplitude
   p1=[fase+(mean(tg)*c) fase+(mean(tp)*c) fase+(tx*c) fase+(tgp*c) fase+(tpp*c)];
elseif d=='t'
  p1=mag_fase(a,b,c,'m');
  p2=mag_fase(a,b,c,'f');
  p3=mag_fase(a,b,c,'g');
  p4=mag_fase(a,b,c,'p');
end

return
