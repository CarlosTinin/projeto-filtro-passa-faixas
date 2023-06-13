function hd=desenha_gabpf(Amax,Amin,w1,w2,w3,w4,w)

% Desenha gabarito PB
% Use hold on antes de chamar a funcao
% Sintaxe: hd=desenha_gabpb(Amax,Amin,wp,ws,w)


h1=line([w(1) w3],[-Amin -Amin]);
set(h1,'color',[1 0 0]);
h2=line([w3 w3],[-Amin 0]);
set(h2,'color',[1 0 0]);
h3=line([w1 w2],[-Amax -Amax]);
set(h3,'color',[1 0 0]);
h4=line([w1 w1],[-Amax -Amin]);
set(h4,'color',[1 0 0]);
h5=line([w2 w2],[-Amax -Amin]);
set(h5,'color',[1 0 0]);
h6=line([w4 w4],[0 -Amin]);
set(h6,'color',[1 0 0]);
h7=line([w4 w(length(w))],[-Amin -Amin]);
set(h7,'color',[1 0 0]);

return















