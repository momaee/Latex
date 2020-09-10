n = [6+3,9+4,12+5,15+6];
t = [4,4,4,4];
t_cent      = [0.479939     ,2.587039   ,10.190753  ,15.677288  ,35.138306  ,44.060038  ,50.342     ,66.490116  ,70.605145  ,95.32];
t_dist      = [0.8440       ,4.4696     ,11.1526    ,16.5418    ,30.1106    ,35.2480    ,40.2736    ,50.1921    ,56.4841    ,60.2560];
t_decent    = [0.7440       ,5.3696     ,8.1526     ,14.5418    ,22.1106    ,30.2480    ,35.2736    ,40.1921    ,45.4841    ,48.2560];
t_heur      = [0.074951     ,0.088694   ,0.0993535  ,0.102600   ,0.110903   ,0.156736   ,0.194503   ,0.262058   ,0.286917   ,0.33];
cost_cent   = [640.4357     ,1.0698e+03 ,719.3515   ,891.9614   ,698.8770   ,1.1122e+03 ,1.0262e+03 ,1.0870e+03 ,734.4871   ,1.0970e+03];
cost_dist   = [640.8357     ,1.0699e+03 ,719.9515   ,890.9614   ,698.9770   ,1.1123e+03 ,1.0261e+03 ,1.0871e+03 ,734.4971   ,1.0970e+03];
cost_decent = [640.4357     ,1.0698e+03 ,719.3515   ,891.9614   ,698.8770   ,1.1122e+03 ,1.0262e+03 ,1.0870e+03 ,734.4871   ,1.0970e+03];
cost_heur   = [644.6828     ,1.0704e+03 ,739.8045   ,893.4754   ,699.9926   ,1.1195e+03 ,1.0367e+03 ,1.0941e+03 ,740.8819   ,1.1070e+03];

n = 50:20:240;
% figure
% plot(n,t_cent,n,t_dist,n,t_decent,n,t_heur)
% legend('centralized','distributed','decentralized','heuristic')
% title('Time of Convergence')
% xlabel('number of all nodes')
% ylabel('time(s)')
figure
plot(n,cost_cent)
hold on
plot(n,cost_dist)
plot(n,cost_decent)
plot(n,cost_heur)

legend('centralized','distributed','decentralized','heuristic', 'location', 'best')
title('Optimization Value')
xlabel('number of all nodes')
ylabel('value')