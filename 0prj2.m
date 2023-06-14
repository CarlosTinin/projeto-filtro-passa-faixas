clear all;
close all;

Amax = 3.01;
Amin = 35; % dB

%  begin - LETRA A
w0 = 2 * pi * 225; % Frequência central em rad

B = 2 * pi * 270; % Largura de banda em rad

w1_eq = [1 B -power(w0, 2)]; % Polinômio que define os valores para w1, já em rad;
w1_roots = roots(w1_eq); % Calcular as raízes do polinômio
w1 = w1_roots(2); % Selecionar raiz positiva

w2 = B + w1;

w3 = 2 * pi * 9;

w4 = power(w0, 2) / w3;

wp = 1; % rad
ws = (w4 - w3) / (w2 - w1); % rad

%  end - LETRA A

%  begin - LETRA B e C

% Determinacao da ordem Butterworth
Nb = order_butter(Amax, Amin, wp, ws);

% Determinacao da ordem Chebyshev
Nc = order_cheby(Amax, Amin, wp, ws);

% Determinacao de T(s) normalizado Butterworth
[nb, db] = aprox_butter(Nb, Amax, wp);

% Determinacao de T(s) normalizado Chebyshev
[nc, dc] = aprox_cheby(Nc, Amax, wp);

% Tranforma PB - > PF
[nbpf, dbpf] = tfpb2pf(nb, db, Amax, B, w0, 'B');
[ncpf, dcpf] = tfpb2pf(nc, dc, Amax, B, w0, 'C');

% Vetor de frequencias para avaliar a magnitude e fase dos filtro normalizados
w = linspace(0.1, 10, 1000); % uma decada acima de wp

% Obtencao da magnitude e fase dos filtros normalizados
[mb, fb] = mag_fase(nb, db, w);
[mc, fc] = mag_fase(nc, dc, w);

% Obtencao da magnitude e fase dos filtros desnormalizados
wd = linspace(w0 / 10, w0 * 10, 1000); %vetor de freq desnormalizados [Hz]
[mbpf, fbpf] = mag_fase(nbpf, dbpf, wd);
[mcpf, fcpf] = mag_fase(ncpf, dcpf, wd);

% Graficos de magnitude comparativo entre filtros normalizados
% Grafico de Butterworth
subplot(2, 1, 1)
semilogx(w, 20 * log10(mb)), grid
xlabel("\\omega [rad/s]")
ylabel("|T(s)|_{dB}");
title(["Magnitude do filtro LP normalizado aproximação de Butterworth (N=" num2str(Nb) ")"])
hold on
desenha_gabpb(Amax, Amin, wp, ws, w);

% Grafico de Chebyshev
subplot(2, 1, 2)
semilogx(w, 20 * log10(mc)), grid
xlabel("\\omega [rad/s]")
ylabel("|T(s)|_{dB}");
title(["Magnitude do filtro LP normalizado aproximação de Chebyshev (N=" num2str(Nc) ")"])
hold on
desenha_gabpb(Amax, Amin, wp, ws, w);

%  end - LETRA B e C

%  begin - LETRA D

tfb = tf(nb, db);
tfc = tf(nc, dc);

% end - LETRA D

% begin - LETRA F

% Grafico da magnitude dos dois filtros com o gabarito normalizado
figure;

subplot(211);
semilogx(w, 20 * log10(mb), w, 20 * log10(mc))
grid on
xlabel("\\omega [rad/s]");
ylabel("|T(\\omega)|_{dB}");
title(["Magnitude de T(s) para filtros normalizados: Nb=" num2str(Nb) ", Nc= " num2str(Nc)])
hold
desenha_gabpb(Amax, Amin, wp, ws, w);
legend("Butterworth", "Chebyshev")

subplot(212);
semilogx(w, fb, w, fc);
grid;
xlabel("\\omega [rad/s]");
ylabel("Fase (\\deg)");
title(["Fase de T(s) para filtros normalizados: Nb=" num2str(Nb) ", Nc= " num2str(Nc)])
legend('Butterworth', 'Chebyshev');

%  end - LETRA F

%  begin - LETRA E

figure;
grid on;
subplot(211);

pzmap(tfb);
axis([-1 0]);
title('Polos e Zeros - Filtro Butterworth normalizado');
xlabel("\\sigma - Eixo Real");
ylabel("j\\omega - Eixo Imaginário");
legend("Pólos");

subplot(212);
pzmap(tfc);
axis([-1 0]);
title('Polos e Zeros - Filtro Chebyshev normalizado');
xlabel("\\sigma - Eixo Real");
ylabel("j\\omega - Eixo Imaginário");
legend("Pólos");

%  end - LETRA E

% begin - LETRA G

tpfb = tf(nbpf, dbpf);
tpfc = tf(ncpf, dcpf);

% função transferência do PF desnormalizado

% end - LETRA G

%  begin - LETRA H

% Graficos de magnitude comparativo entre filtros desnormalizados
% Grafico de Butterworth
figure
subplot(2, 1, 1)
semilogx(wd / (2 * pi), 20 * log10(mbpf)), grid
xlabel("Frequência [Hz]")
ylabel("|T(s)|_{dB}");
title(["Magnitude do filtro PF aproximação de Butterworth (N=" num2str(Nb) ")"])
hold on
desenha_gabpf(Amax, Amin, w1 / (2 * pi), w2 / (2 * pi), w3 / (2 * pi), w4 / (2 * pi), wd ./ (2 * pi));

% Grafico de Chebyshev
subplot(2, 1, 2)
semilogx(wd / (2 * pi), 20 * log10(mcpf)), grid
xlabel("Frequência [Hz]")
ylabel("|T(s)|_{dB}");
title(["Magnitude do filtro PF aproximação de Chebyshev (N=" num2str(Nc) ")"])
hold on
desenha_gabpf(Amax, Amin, w1 / (2 * pi), w2 / (2 * pi), w3 / (2 * pi), w4 / (2 * pi), wd ./ (2 * pi));

%  end - LETRA H

%  begin - LETRA I

figure;
hold on;

subplot(211);
semilogx(wd / (2 * pi), fbpf);
axis([w1 w2]);
grid on;
xlabel("\\omega [rad/s]");
ylabel("Fase (\\deg)");
title(["Fase de T(s) Butterworth desnormalizado"])

subplot(212);
semilogx(wd / (2 * pi), fcpf);
axis([w1 w2]);
grid on;
xlabel("\\omega [rad/s]");
ylabel("Fase (\\deg)");
title(["Fase de T(s) Chebyshev desnormalizado"])

% ! somente banda passante

%  end - LETRA I

%  begin - LETRA J

figure;
[y, t] = impulse(tpfb);

subplot(211);
plot(t, y);
grid;
title('Resposta ao Impulso do Filtro Butterworth desnormalizado');
xlabel('Time [s]');
ylabel('Amplitude');

[y, t] = impulse(tpfc);

subplot(212);
plot(t, y);
grid;

title('Resposta ao Impulso do Filtro Chebyshev desnormalizado');
xlabel('Time [s]');
ylabel('Amplitude');

%  end - LETRA J

% Criacao do sinal a ser usado no lsim
% Um sinal senoidal de 5 rad/s + 1 rad/s + 10 rad/s

% Frequencia de amostragem para definir o vetor temporal
% No minimo 10 x a maior frequencia do sinal contido na banda
Ws = 10 * w2; % freq de amostragem no minimo 10x a maxima freq na banda
Ts = (2 * pi) / Ws; % periodo de amostragem
L = 2048; % numero de amostras
t = (0:L - 1) * Ts; % vetor de tempo

sin1 = sin(2 * pi * 9 * t);
sin2 = sin(w0 * t);
sin3 = sin(2 * pi * 1000 * t);
x = sin1 + sin2 + sin3; % sinal composto de 3 senoides

% Obtem a transformada de Fourier do sinal
X = fft(x);
mX = abs(X / L); %magnitude da transformada
mdX = mX(1:L / 2 +1); %magnitude deslocada
mdX(2:length(mdX) - 1) = 2 * mdX(2:length(mdX) - 1);

% vetor de frequencias
w = (Ws / L) * (0:(L / 2));

% Funcao de tranferencia dos filtros PF
Sb = tf(nbpf, dbpf); %filtro PF Butterworth
Sc = tf(ncpf, dcpf); %filtro PF Chebyshev

Sb = tpfb;
Sb = tpfc;

% Simulacao da TF com o sinal x(t) para ambas aproximacoes
yb = lsim(Sb, x, t);
yc = lsim(Sc, x, t);

%Obtem a tranformada do sinal apos a simulacao

Yb = fft(yb);
mYb = abs(Yb / L); %magnitude da transformada
mdYb = mYb(1:L / 2 +1); %magnitude deslocada
mdYb(2:length(mdYb) - 1) = 2 * mdYb(2:length(mdYb) - 1);

Yc = fft(yc);
mYc = abs(Yc / L); %magnitude da transformada
mdYc = mYc(1:L / 2 +1); %magnitude deslocada
mdYc(2:length(mdYc) - 1) = 2 * mdYc(2:length(mdYc) - 1);

% Graficos da FFT dos sinais

mdX(mdX <= 0.05) = NaN;
mdYb(mdYb <= 0.05) = NaN;
mdYc(mdYc <= 0.05) = NaN;

figure
subplot(3, 1, 1)
stem(w, mdX)
title("Amplitude espectral do sinal x(t)");
xlabel("\\omega [rad/s]");
ylabel("|X(\\omega)|");

subplot(3, 1, 2)
stem(w, mdYb)
title("Amplitude espectral do sinal y(t) apos o filtro PF Butterworth");
xlabel("\\omega [rad/s]");
ylabel("|Y(\\omega)|");

subplot(3, 1, 3)
stem(w, mdYc)
title("Amplitude espectral do sinal y(t) apos o filtro PF Chebyshev");
xlabel("\\omega [rad/s]");
ylabel("|Y(\\omega)|");

% Sinais de entrada (x(t)) e saida (y(t)) apos o processo de filtragem
% Aproximacao de Butterworth
figure
subplot(4, 1, 1)
plot(Ws * t(100:300), sin1(100:300))
xlabel("tempo [s]");
ylabel("x_1 (t)");
title("Sinal x_1(t)=sen(0.5\\omega_0t)");
subplot(4, 1, 2)
plot(Ws * t(100:300), sin2(100:300))
xlabel("tempo [s]");
ylabel("x_2 (t)");
title("Sinal x_2(t)=sen(\\omega_0t)");
subplot(4, 1, 3)
plot(Ws * t(100:300), sin3(100:300))
xlabel("tempo [s]");
ylabel("x_3 (t)");
title("Sinal x_3(t)=sen(3\\omega_0t)");
subplot(4, 1, 4)
plot(Ws * t(100:300), x(100:300))
xlabel("tempo [s]");
ylabel("x (t)");
title("Sinal x(t)=sen(0.5\\omega_0t)+sen(\\omega_0t)+sen(3\\omega_0t)");

% Sinais de entrada (x(t)) e saida (y(t)) apos o processo de filtragem
% Aproximacoes de Butterworth e Chebyshev
figure
subplot(2, 1, 1)
plot(Ws * t(100:300), sin2(100:300), Ws * t(100:300), yb(100:300))
xlabel("tempo [s]");
ylabel("x_2(t) e y(t)");
title("Sinais x_2(t) e y(t) após o processo de filtragem usando a aproximação de Butterworth");

subplot(2, 1, 2)
plot(Ws * t(100:300), sin2(100:300), Ws * t(100:300), yc(100:300))
xlabel("tempo [s]");
ylabel("x_2(t) e y(t)");
title("Sinais x_2(t) e y(t) após o processo de filtragem usando a aproximação de Chebyshev");
