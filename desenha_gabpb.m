function hd=desenha_gabpb(Amax,Amin,wp,ws,w)

% Desenha gabarito PB
% Use hold on antes de chamar a funcao
% Sintaxe: hd=desenha_gabpb(Amax,Amin,wp,ws,w)

g=find(w<wp);
vamx=ones(length(g))*(-Amax);
h1=line([w(1) wp],[-Amax -Amax]);
set(h1,'color',[0 0 0]);
h2=line([wp wp],[-Amax -Amin]);
set(h2,'color',[0 0 0]);
h3=line([ws ws],[0 -Amin]);
set(h3,'color',[0 0 0]);
h4=line([ws w(length(w))],[-Amin -Amin]);
set(h4,'color',[0 0 0]);

return















