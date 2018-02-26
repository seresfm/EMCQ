
%Ejecucion del algortimo para EMCA paper
%by Sergio Gonzalez

clearvars
clc
seed = [479 122 489 169 363 44 265 37 67 124 34 112 390 46 119 69 99 450 77 191];
N = 100;
NoMBS = 2;
K = 500;

p = 0.1           
r1 = 0.5;
r2 = 0.2;
tau_type = 9
 


disp ( ' seed   K(long. max. buffer)  Pa1_h   Pa_plus_h  Pa1_l  Pa_plus_l   L_T  D_T  Throughput Psuccess Pfree  Pcollision   Prxbussy Pp');

for i = 1: size(seed,2)

[ Pa1_h, Pa_plus_h, Pa1_l, Pa_plus_l, L_T, D_T, Throughput, Psuccess, Pfree, Pcollision, Prxbussy, Pp] = EMCA2(seed(i), p, N, NoMBS, r1, r2, tau_type, K);

format short g
Results = [seed(i) K Pa1_h Pa_plus_h Pa1_l Pa_plus_l  L_T D_T  Throughput Psuccess Pfree Pcollision Prxbussy Pp];
disp (Results)

end

p = 0.15
disp ( ' seed   K(long. max. buffer)  Pa1_h   Pa_plus_h  Pa1_l  Pa_plus_l   L_T  D_T  Throughput Psuccess Pfree  Pcollision   Prxbussy Pp');

for i = 1: size(seed,2)


[ Pa1_h, Pa_plus_h, Pa1_l, Pa_plus_l, L_T, D_T, Throughput, Psuccess, Pfree, Pcollision, Prxbussy, Pp] = EMCA2(seed(i), p, N, NoMBS, r1, r2, tau_type, K);

format short g
Results = [seed(i) K Pa1_h Pa_plus_h Pa1_l Pa_plus_l  L_T D_T  Throughput Psuccess Pfree Pcollision Prxbussy Pp];
disp (Results)

end
p = 0.2
disp ( ' seed   K(long. max. buffer)  Pa1_h   Pa_plus_h  Pa1_l  Pa_plus_l   L_T  D_T  Throughput Psuccess Pfree  Pcollision   Prxbussy Pp');

for i = 1: size(seed,2)


[ Pa1_h, Pa_plus_h, Pa1_l, Pa_plus_l, L_T, D_T, Throughput, Psuccess, Pfree, Pcollision, Prxbussy, Pp] = EMCA2(seed(i), p, N, NoMBS, r1, r2, tau_type, K);

format short g
Results = [seed(i) K Pa1_h Pa_plus_h Pa1_l Pa_plus_l  L_T D_T  Throughput Psuccess Pfree Pcollision Prxbussy Pp];
disp (Results)

end
p = 0.25
disp ( ' seed   K(long. max. buffer)  Pa1_h   Pa_plus_h  Pa1_l  Pa_plus_l   L_T  D_T  Throughput Psuccess Pfree  Pcollision   Prxbussy Pp');

for i = 1: size(seed,2)


[ Pa1_h, Pa_plus_h, Pa1_l, Pa_plus_l, L_T, D_T, Throughput, Psuccess, Pfree, Pcollision, Prxbussy, Pp] = EMCA2(seed(i), p, N, NoMBS, r1, r2, tau_type, K);

format short g
Results = [seed(i) K Pa1_h Pa_plus_h Pa1_l Pa_plus_l  L_T D_T  Throughput Psuccess Pfree Pcollision Prxbussy Pp];
disp (Results)

end
p = 0.3
disp ( ' seed   K(long. max. buffer)  Pa1_h   Pa_plus_h  Pa1_l  Pa_plus_l   L_T  D_T  Throughput Psuccess Pfree  Pcollision   Prxbussy Pp');

for i = 1: size(seed,2)


[ Pa1_h, Pa_plus_h, Pa1_l, Pa_plus_l, L_T, D_T, Throughput, Psuccess, Pfree, Pcollision, Prxbussy, Pp] = EMCA2(seed(i), p, N, NoMBS, r1, r2, tau_type, K);

format short g
Results = [seed(i) K Pa1_h Pa_plus_h Pa1_l Pa_plus_l  L_T D_T  Throughput Psuccess Pfree Pcollision Prxbussy Pp];
disp (Results)

end
p = 0.35
disp ( ' seed   K(long. max. buffer)  Pa1_h   Pa_plus_h  Pa1_l  Pa_plus_l   L_T  D_T  Throughput Psuccess Pfree  Pcollision   Prxbussy Pp');

for i = 1: size(seed,2)


[ Pa1_h, Pa_plus_h, Pa1_l, Pa_plus_l, L_T, D_T, Throughput, Psuccess, Pfree, Pcollision, Prxbussy, Pp] = EMCA2(seed(i), p, N, NoMBS, r1, r2, tau_type, K);

format short g
Results = [seed(i) K Pa1_h Pa_plus_h Pa1_l Pa_plus_l  L_T D_T  Throughput Psuccess Pfree Pcollision Prxbussy Pp];
disp (Results)

end
p = 0.4
disp ( ' seed   K(long. max. buffer)  Pa1_h   Pa_plus_h  Pa1_l  Pa_plus_l   L_T  D_T  Throughput Psuccess Pfree  Pcollision   Prxbussy Pp');

for i = 1: size(seed,2)


[ Pa1_h, Pa_plus_h, Pa1_l, Pa_plus_l, L_T, D_T, Throughput, Psuccess, Pfree, Pcollision, Prxbussy, Pp] = EMCA2(seed(i), p, N, NoMBS, r1, r2, tau_type, K);

format short g
Results = [seed(i) K Pa1_h Pa_plus_h Pa1_l Pa_plus_l  L_T D_T  Throughput Psuccess Pfree Pcollision Prxbussy Pp];
disp (Results)

end
p = 0.45
disp ( ' seed   K(long. max. buffer)  Pa1_h   Pa_plus_h  Pa1_l  Pa_plus_l   L_T  D_T  Throughput Psuccess Pfree  Pcollision   Prxbussy Pp');

for i = 1: size(seed,2)


[ Pa1_h, Pa_plus_h, Pa1_l, Pa_plus_l, L_T, D_T, Throughput, Psuccess, Pfree, Pcollision, Prxbussy, Pp] = EMCA2(seed(i), p, N, NoMBS, r1, r2, tau_type, K);

format short g
Results = [seed(i) K Pa1_h Pa_plus_h Pa1_l Pa_plus_l  L_T D_T  Throughput Psuccess Pfree Pcollision Prxbussy Pp];
disp (Results)

end
p = 0.5
disp ( ' seed   K(long. max. buffer)  Pa1_h   Pa_plus_h  Pa1_l  Pa_plus_l   L_T  D_T  Throughput Psuccess Pfree  Pcollision   Prxbussy Pp');

for i = 1: size(seed,2)


[ Pa1_h, Pa_plus_h, Pa1_l, Pa_plus_l, L_T, D_T, Throughput, Psuccess, Pfree, Pcollision, Prxbussy, Pp] = EMCA2(seed(i), p, N, NoMBS, r1, r2, tau_type, K);

format short g
Results = [seed(i) K Pa1_h Pa_plus_h Pa1_l Pa_plus_l  L_T D_T  Throughput Psuccess Pfree Pcollision Prxbussy Pp];
disp (Results)

end
p = 0.55
disp ( ' seed   K(long. max. buffer)  Pa1_h   Pa_plus_h  Pa1_l  Pa_plus_l   L_T  D_T  Throughput Psuccess Pfree  Pcollision   Prxbussy Pp');

for i = 1: size(seed,2)


[ Pa1_h, Pa_plus_h, Pa1_l, Pa_plus_l, L_T, D_T, Throughput, Psuccess, Pfree, Pcollision, Prxbussy, Pp] = EMCA2(seed(i), p, N, NoMBS, r1, r2, tau_type, K);

format short g
Results = [seed(i) K Pa1_h Pa_plus_h Pa1_l Pa_plus_l  L_T D_T  Throughput Psuccess Pfree Pcollision Prxbussy Pp];
disp (Results)

end
p = 0.6
disp ( ' seed   K(long. max. buffer)  Pa1_h   Pa_plus_h  Pa1_l  Pa_plus_l   L_T  D_T  Throughput Psuccess Pfree  Pcollision   Prxbussy Pp');

for i = 1: size(seed,2)


[ Pa1_h, Pa_plus_h, Pa1_l, Pa_plus_l, L_T, D_T, Throughput, Psuccess, Pfree, Pcollision, Prxbussy, Pp] = EMCA2(seed(i), p, N, NoMBS, r1, r2, tau_type, K);

format short g
Results = [seed(i) K Pa1_h Pa_plus_h Pa1_l Pa_plus_l  L_T D_T  Throughput Psuccess Pfree Pcollision Prxbussy Pp];
disp (Results)

end
p = 0.65
disp ( ' seed   K(long. max. buffer)  Pa1_h   Pa_plus_h  Pa1_l  Pa_plus_l   L_T  D_T  Throughput Psuccess Pfree  Pcollision   Prxbussy Pp');

for i = 1: size(seed,2)


[ Pa1_h, Pa_plus_h, Pa1_l, Pa_plus_l, L_T, D_T, Throughput, Psuccess, Pfree, Pcollision, Prxbussy, Pp] = EMCA2(seed(i), p, N, NoMBS, r1, r2, tau_type, K);

format short g
Results = [seed(i) K Pa1_h Pa_plus_h Pa1_l Pa_plus_l  L_T D_T  Throughput Psuccess Pfree Pcollision Prxbussy Pp];
disp (Results)

end
p = 0.7
disp ( ' seed   K(long. max. buffer)  Pa1_h   Pa_plus_h  Pa1_l  Pa_plus_l   L_T  D_T  Throughput Psuccess Pfree  Pcollision   Prxbussy Pp');

for i = 1: size(seed,2)


[ Pa1_h, Pa_plus_h, Pa1_l, Pa_plus_l, L_T, D_T, Throughput, Psuccess, Pfree, Pcollision, Prxbussy, Pp] = EMCA2(seed(i), p, N, NoMBS, r1, r2, tau_type, K);

format short g
Results = [seed(i) K Pa1_h Pa_plus_h Pa1_l Pa_plus_l  L_T D_T  Throughput Psuccess Pfree Pcollision Prxbussy Pp];
disp (Results)

end
p = 0.75
disp ( ' seed   K(long. max. buffer)  Pa1_h   Pa_plus_h  Pa1_l  Pa_plus_l   L_T  D_T  Throughput Psuccess Pfree  Pcollision   Prxbussy Pp');

for i = 1: size(seed,2)


[ Pa1_h, Pa_plus_h, Pa1_l, Pa_plus_l, L_T, D_T, Throughput, Psuccess, Pfree, Pcollision, Prxbussy, Pp] = EMCA2(seed(i), p, N, NoMBS, r1, r2, tau_type, K);

format short g
Results = [seed(i) K Pa1_h Pa_plus_h Pa1_l Pa_plus_l  L_T D_T  Throughput Psuccess Pfree Pcollision Prxbussy Pp];
disp (Results)

end
p = 0.8
disp ( ' seed   K(long. max. buffer)  Pa1_h   Pa_plus_h  Pa1_l  Pa_plus_l   L_T  D_T  Throughput Psuccess Pfree  Pcollision   Prxbussy Pp');

for i = 1: size(seed,2)


[ Pa1_h, Pa_plus_h, Pa1_l, Pa_plus_l, L_T, D_T, Throughput, Psuccess, Pfree, Pcollision, Prxbussy, Pp] = EMCA2(seed(i), p, N, NoMBS, r1, r2, tau_type, K);

format short g
Results = [seed(i) K Pa1_h Pa_plus_h Pa1_l Pa_plus_l  L_T D_T  Throughput Psuccess Pfree Pcollision Prxbussy Pp];
disp (Results)

end
p = 0.85
disp ( ' seed   K(long. max. buffer)  Pa1_h   Pa_plus_h  Pa1_l  Pa_plus_l   L_T  D_T  Throughput Psuccess Pfree  Pcollision   Prxbussy Pp');

for i = 1: size(seed,2)


[ Pa1_h, Pa_plus_h, Pa1_l, Pa_plus_l, L_T, D_T, Throughput, Psuccess, Pfree, Pcollision, Prxbussy, Pp] = EMCA2(seed(i), p, N, NoMBS, r1, r2, tau_type, K);

format short g
Results = [seed(i) K Pa1_h Pa_plus_h Pa1_l Pa_plus_l  L_T D_T  Throughput Psuccess Pfree Pcollision Prxbussy Pp];
disp (Results)

end
p = 0.9
disp ( ' seed   K(long. max. buffer)  Pa1_h   Pa_plus_h  Pa1_l  Pa_plus_l   L_T  D_T  Throughput Psuccess Pfree  Pcollision   Prxbussy Pp');

for i = 1: size(seed,2)


[ Pa1_h, Pa_plus_h, Pa1_l, Pa_plus_l, L_T, D_T, Throughput, Psuccess, Pfree, Pcollision, Prxbussy, Pp] = EMCA2(seed(i), p, N, NoMBS, r1, r2, tau_type, K);

format short g
Results = [seed(i) K Pa1_h Pa_plus_h Pa1_l Pa_plus_l  L_T D_T  Throughput Psuccess Pfree Pcollision Prxbussy Pp];
disp (Results)

end
p = 0.95
disp ( ' seed   K(long. max. buffer)  Pa1_h   Pa_plus_h  Pa1_l  Pa_plus_l   L_T  D_T  Throughput Psuccess Pfree  Pcollision   Prxbussy Pp');

for i = 1: size(seed,2)


[ Pa1_h, Pa_plus_h, Pa1_l, Pa_plus_l, L_T, D_T, Throughput, Psuccess, Pfree, Pcollision, Prxbussy, Pp] = EMCA2(seed(i), p, N, NoMBS, r1, r2, tau_type, K);

format short g
Results = [seed(i) K Pa1_h Pa_plus_h Pa1_l Pa_plus_l  L_T D_T  Throughput Psuccess Pfree Pcollision Prxbussy Pp];
disp (Results)

end
p = 1
disp ( ' seed   K(long. max. buffer)  Pa1_h   Pa_plus_h  Pa1_l  Pa_plus_l   L_T  D_T  Throughput Psuccess Pfree  Pcollision   Prxbussy Pp');

for i = 1: size(seed,2)


[ Pa1_h, Pa_plus_h, Pa1_l, Pa_plus_l, L_T, D_T, Throughput, Psuccess, Pfree, Pcollision, Prxbussy, Pp] = EMCA2(seed(i), p, N, NoMBS, r1, r2, tau_type, K);

format short g
Results = [seed(i) K Pa1_h Pa_plus_h Pa1_l Pa_plus_l  L_T D_T  Throughput Psuccess Pfree Pcollision Prxbussy Pp];
disp (Results)

end