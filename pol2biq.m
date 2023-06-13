function  [n,d] = pol2biq(a,b)

%POL2BIQ Converte o conjunto de funções polinômiais em biquadráticas 
%	  Sintaxe: [numbiq,denbiq] = pol2biq(numpol,denpol)
%         numbiq = numerador da função de transferência no formato biquadrática
%         denbiq = denominador da função de transferência no formato biquadrática
%	  numpol = numerador da função de transferência no formato polinomial
%	  denpol = denominador da função de transferência no formato polinomial
%		versão: 2
%	Copyright (c) 2015 by LaPeSC/Delmar B. Carvalho.

rza=roots(a);
rzb=roots(b);
[la,c]=size(rza);
[lb,c]=size(rzb);
ttrz=(lb/2) - fix(lb/2);
	if (ttrz > 0.0)
		nrzimpar=1;
		nrzpar=0;
	else
		nrzimpar=0;
		nrzpar=1;
	endif
	if (nrzimpar==1)
		n(1,:)=[0 0 abs(rzb(1))];
		d(1,:)=[0 1 abs(rzb(1))];
		k=2;
		offsetn=2;
		offsetd=2;
	else
		k=1;
		offsetn=1;
		offsetd=1;
		
	endif
	if (la==0)
	  ln=length(b);
	  cte=nthroot(b(ln),ln);
	  for i=k:lb-1
	    n(i,:)=[0 0 cte];
	    offsetn=offsetn+2;
	    if (offsetn>lb-1)
	      break;
	    endif  
	  endfor
	elseif (la == lb)
	    for i=k:la-1
		n(i,:)=real(poly(rza(offsetn:offsetn+1)));
		offsetn=offsetn+2;
		if (offsetn>la-1)
		   break;
		endif
	    endfor	
	elseif (la < lb) 
	  ganhok=max(a);
	  cte=nthroot(ganhok,la);
	  for i=k:la
	    n(i,:)=[0 cte 0];
	  endfor
	endif  
%monta biquadraticas denominador		
	for i=k:lb-1
	  d(i,:)=real(poly(rzb(offsetd:offsetd+1)));
	  offsetd=offsetd+2;
	  if (offsetd>lb-1)
	    break;
	  endif
	endfor	
	
return
