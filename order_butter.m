function order=order_butter(Amax,Amin,wp,ws)

% Aproximacao de Butterworth
% Retorna o ordem de um filtro LP Butterworth
% Sintaxe: order=order_butter(Amax,Amin,wp,ws)

w=ws/wp;
t1=10^(0.1*Amin)-1;
t2=10^(0.1*Amax)-1;
order=ceil(log10(t1/t2)/(2*log10(w)));

return















