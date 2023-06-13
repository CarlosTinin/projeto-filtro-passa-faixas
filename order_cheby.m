function order=order_cheby(Amax,Amin,wp,ws)

% Aproximacao de Chebyshev
% Retorna o ordem de um filtro LP Chebyshev
% Sintaxe: order=order_cheby(Amax,Amin,wp,ws)

w=ws/wp;
t1=10^(0.1*Amin)-1;
t2=10^(0.1*Amax)-1;
order=ceil(acosh(sqrt(t1/t2))/acosh(w));

return















