clf
hold on
%%%%%%%%%%%%%%%
%subproblem parallelism
% n = 128
ser = 0.000124;


%for frames=1
t1 = [0.0001    0.0002    0.0002    0.0003    0.0027    0.0117];
p1 = [1         2         4         8         16        32    ];
sp1 = ser./t1;

t2 = [0.0001349	0.0002	0.0004089	0.01406	0.02433];
p2 = [2         4         8         16        32    ];
sp2 = 2*ser./t2;

t4 = [0.0002389	0.0005	0.002793	0.01615];
p4 = [4         8         16        32    ];
sp4 = 4*ser./t4;

t8 = [0.0002732	0.01049	0.01138	0.03823];
p8 = [8         16      32      64      ];
sp8 = 8*ser./t8;


t16 = [0.000526	0.01642	0.02522];
p16 = [16       32      64      ];
sp16 = 16*ser./t16;

%plot(p1,t1, 'r',p2, sp2, 'g', p4, sp4, 'b', p8, sp8, 'r-.', p16, sp16, 'b--' )
semilogy(p1,t1, 'r',p2, sp2, 'g', p4, sp4, 'b', p8, sp8, 'r-.', p16, sp16, 'b--' )
legend('1 Frame', '2 Frames', '4 Frames', '8 Frames', '16 Frames')

title('Negative Speedup from Parallelizing Subproblems, n = 128','FontWeight', 'bold')
xlabel('Number of processes')
ylabel('Speedup')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%time frame parallelism
clf
hold on


p = [1         2         4         8         16];
t64=[7.80E-005 8.01E-005 9.58E-005 0.0001571 0.001164];
ser = 0.000124;
t128 = [0.000124 0.0001349 0.0002389 0.0002732 0.000526];
f = [1 2 4 8 16];
sp128 = ser.*f./t128;
ser = 7.80E-005;
sp64 = ser.*f./t64;

semilogy(p, sp64, 'r', p, sp128, 'b--')

title('Speedup from Parallelizing by Time Frame','FontWeight', 'bold')
xlabel('Number of processes')
ylabel('Speedup')
legend('n=64', 'n=128')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%process/frame
%n=128
clf
hold on
ser = 0.000124;
pperf = [1         2         4         8         16];
f1 =    [0.000124 0.0001531 0.000164 0.0002611 0.002719];
f2 =    [0.0001349 0.0002	0.0004089	0.01406	0.02433];
f4 =    [0.0002389	0.0005	0.002793	0.01615];
f8 =    [0.0002732	0.01049	0.01138	0.03823];
f16 =   [0.000526	0.01642	0.02522];

spf1 = ser./f1;
spf2 = 2*ser./f2;
spf4 = 4*ser./f4;
spf8 = 8*ser./f8;
spf16 = 16*ser./f16;


semilogy(pperf,spf1, 'r',pperf, spf2, 'g', pperf(1:4), spf4, 'b', pperf(1:4), spf8, 'r-.', pperf(1:3), spf16, 'b--' )

legend('1 Frame', '2 Frames', '4 Frames', '8 Frames', '16 Frames')

title('Negative Speedup by Processes per Frame, n = 128','FontWeight', 'bold')
xlabel('Number of processes per Frame')
ylabel('Speedup')





