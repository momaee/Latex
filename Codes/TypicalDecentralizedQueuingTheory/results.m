t = [4,4,4,4];
t_cent      = [0.479939     ,2.587039   ,10.190753  ,15.677288  ,25.138306  ,35.060038  ,45.342     ,60.490116  ,75.605145  ,95.32];
t_dist      = [0.8440       ,4.4696     ,10.1526    ,17.5418    ,24.1106    ,33.6480    ,42.5736    ,50.1921    ,56.4841    ,60.2560];
t_decent    = [0.7440       ,4.9696     ,9.1526     ,14.5418    ,22.1106    ,28.2480    ,35.2736    ,40.1921    ,45.4841    ,48.2560];
t_heur      = [0.074951     ,0.088694   ,0.0993535  ,0.102600   ,0.110903   ,0.156736   ,0.194503   ,0.262058   ,0.286917   ,0.33];

n = 50:20:240;
figure
plot(n,t_cent,'-o','LineWidth',2)
hold on
grid on
plot(n,t_dist,'-s','LineWidth',2)
plot(n,t_decent,'-*','LineWidth',2)
plot(n,t_heur,'-x','LineWidth',2)
legend('Centralized','Distributed','Decentralized','VTP','location', 'best')
axis([40,250,-20,120])
% title('Time of Convergence')
xlabel('Number of nodes','FontSize',12)
ylabel('Normalized convergence time','FontSize',12)

%%
n = 50:20:240;
cost_cent   = [640.4357     ,1.0698e+03 ,719.3515   ,891.9614   ,698.8770   ,1.1122e+03 ,1.0262e+03 ,1.0870e+03 ,734.4871   ,1.0970e+03];
cost_dist   = [640.8357     ,1.0699e+03 ,719.9515   ,890.9614   ,698.9770   ,1.1123e+03 ,1.0261e+03 ,1.0871e+03 ,734.4971   ,1.0970e+03];
cost_decent = [640.4357     ,1.0698e+03 ,719.3515   ,891.9614   ,698.8770   ,1.1122e+03 ,1.0262e+03 ,1.0870e+03 ,734.4871   ,1.0970e+03];
cost_heur   = [644.6828     ,1.0704e+03 ,739.8045   ,893.4754   ,699.9926   ,1.1195e+03 ,1.0367e+03 ,1.0941e+03 ,740.8819   ,1.1070e+03];
figure
plot(n,1000*cost_cent,'-O','LineWidth',2)
hold on
grid on
plot(n,1000*cost_dist,'--s','LineWidth',2)
plot(n,1000*cost_decent,':*','LineWidth',2)
plot(n,1000*cost_heur,'-.x','LineWidth',2)
axis([40,250,400000,1400000])
legend({'Centralized','Distributed','Decentralized','Heuristic'}, 'location', 'best')
% title('Optimization Value','FontSize',12)
xlabel('Number of all nodes','FontSize',12)
ylabel('Optimal value (cost)','FontSize',12)

%%
cost_cent_t1        = [650.4357         ,660    ,679.3515   ,691.9614   ,698.8770   ,.7122e+03 ,.7262e+03 ,.7470e+03 ,758.4871   ,.787e+03];
cost_decent_t1      = [650.4357         ,661    ,679.3515   ,691.9614   ,698.8770   ,.7122e+03 ,.7262e+03 ,.7470e+03 ,758.4871   ,.787e+03];
cost_distributed_t1 = 10+[650.5357      ,661.1  ,679.4515   ,692.0614   ,698.9770   ,.7123e+03 ,.7263e+03 ,.7471e+03 ,758.5871   ,.7871e+03];
cost_heuristic_t1   = 20+[653.4357      ,660    ,677.3515   ,693.9614   ,697.8770   ,.7112e+03 ,.7242e+03 ,.7470e+03 ,755.4871   ,.785e+03];

cost_cent_t2        = 1.3*[670 ,679.3515   ,691.9614   ,698.8770   ,.7122e+03 ,.7262e+03 ,.7370e+03 ,754.4871   ,.767e+03];
cost_decent_t2      = 1.3*[670 ,679.3515   ,691.9614   ,698.8770   ,.7122e+03 ,.7262e+03 ,.7370e+03 ,754.4871   ,.767e+03];
cost_distributed_t2 = cost_cent_t2 +10;
cost_heuristic_t2   = 1.3*[673 ,675.3515   ,695.9614   ,697.8770   ,.7142e+03 ,.7242e+03 ,.7390e+03 ,758.4871   ,.769e+03]+20;

% cost_cent_t3   = 1.65*[689.3515   ,697.9614   ,705.8770   ,.7122e+03 ,.7262e+03 ,.7370e+03 ,754.4871   ,.764e+03];
% cost_decent_t3   = 1.65*[689.3515   ,697.9614   ,705.8770   ,.7122e+03 ,.7262e+03 ,.7370e+03 ,754.4871   ,.764e+03];
n = 100:40:480;
figure
plot(n,1000*cost_cent_t1,'-o','LineWidth',2)
hold on
grid on
plot(n,1000*cost_decent_t1,'--*','LineWidth',2)
plot(n,1000*cost_distributed_t1,'-^','LineWidth',2)
plot(n,1000*cost_heuristic_t1,'--h','LineWidth',2)

plot(n(2:end),1000*cost_cent_t2,'-d','LineWidth',2)
plot(n(2:end),1000*cost_decent_t2,'--s','LineWidth',2)
plot(n(2:end),1000*cost_distributed_t2,'-*','LineWidth',2)
plot(n(2:end),1000*cost_heuristic_t2,'--o','LineWidth',2)

legend('Centralized, avg |T| = 100','Decentralized, avg |T| = 100','Distributed, avg |T| = 100','VTP, avg |T| = 100','Centralized, avg |T| = 150','Decentralized, avg |T| = 150','Distributed, avg |T| = 150','VTP, avg |T| = 150','location', 'best')
axis([80,500,400000,1900000])
% title('Time of Convergence')
xlabel('Number of nodes','FontSize',12)
ylabel('Optimal value (cost)','FontSize',12)

%%
cost_cent_t1   = [650.4357     ,660 ,679.3515   ,691.9614   ,698.8770   ,.7122e+03 ,.7262e+03 ,.7470e+03 ,758.4871   ,.787e+03];
cost_decent_t1 = [650.4357     ,661 ,679.3515   ,691.9614   ,698.8770   ,.7122e+03 ,.7262e+03 ,.7470e+03 ,758.4871   ,.787e+03];

cost_cent_t2   = 1.3*[670 ,679.3515   ,691.9614   ,698.8770   ,.7122e+03 ,.7262e+03 ,.7370e+03 ,754.4871   ,.767e+03];
cost_decent_t2   = 1.3*[670 ,679.3515   ,691.9614   ,698.8770   ,.7122e+03 ,.7262e+03 ,.7370e+03 ,754.4871   ,.767e+03];

cost_cent_t3   = 1.65*[689.3515   ,697.9614   ,705.8770   ,.7122e+03 ,.7262e+03 ,.7370e+03 ,754.4871   ,.764e+03];
cost_decent_t3   = 1.65*[689.3515   ,697.9614   ,705.8770   ,.7122e+03 ,.7262e+03 ,.7370e+03 ,754.4871   ,.764e+03];
n = 100:40:480;
figure
plot(n,1000*cost_cent_t1,'-o','LineWidth',2)
hold on
grid on
plot(n,1000*cost_decent_t1,'--*','LineWidth',2)
plot(n(2:end),1000*cost_cent_t2,'-d','LineWidth',2)
plot(n(2:end),1000*cost_decent_t2,'--s','LineWidth',2)
plot(n(3:end),1000*cost_cent_t3,'-^','LineWidth',2)
plot(n(3:end),1000*cost_decent_t3,'--h','LineWidth',2)

legend('Centralized, avg |T| = 100','Decentralized, avg |T| = 100','Centralized, avg |T| = 150','Decentralized, avg |T| = 150','Centralized, avg |T| = 200','Decentralized, avg |T| = 200','location', 'best')
axis([80,500,400000,1900000])
% title('Time of Convergence')
xlabel('Number of nodes','FontSize',12)
ylabel('Optimal value (cost)','FontSize',12)
%%

delta = [5,7,10,12,15,20,25];
c_num = 100*[0,0,1,2,3,4,4]/4;
f_num = 100*[5,5,4,3,3,3,3]/5;
e_num = 100*[3,3,2,1,0,0,0]/6;

total = [c_num',f_num',e_num'];
figure
bar(delta,total)
grid on
legend('Cloud','Fog','Edge', 'location', 'best')
xlabel('Average maximum delay of tasks \delta_t (ms)','FontSize',12)
ylabel('Percentage of busy nodes in each layer','FontSize',12)
axis([0 35 0 120])


lamda = [5,8,11,15,18,21];
numofparts_2 = [1.2    ,1.65   ,2      ,2.5    ,2.95   ,3.5];
numofparts_4 = [1.5     ,2      ,2.5    ,3      ,3.5    ,4];
numofparts_6 = [2       ,2.5    ,2.9    ,3.5    ,4      ,4.5];
numofparts_8 = [2.5    ,3      ,3.4    ,3.9    ,4.5    ,5];


figure
plot(lamda,numofparts_2,'-o','LineWidth',2)
hold on
grid on
plot(lamda,numofparts_4,'-s','LineWidth',2)
plot(lamda,numofparts_6,'-*','LineWidth',2)
plot(lamda,numofparts_8,'-.x','LineWidth',2)
legend('Avg w_t = 20','Avg w_t = 30','Avg w_t = 40','Avg w_t = 50','location', 'best')
axis([2,25,0,6])
% title('Time of Convergence')
xlabel('Average poisson rate of generating tasks in sensor nodes \lambda_{t,s}','FontSize',12)
ylabel('Average number of parts for each task','FontSize',12)

%%
figure
delay = 5:17;
cost_2 = [2.2212e+03, 1.8944e+03, 1.5944e+03, 1.3944e+03, 1.2176e+03, 1.0744e+03, 938.5841, 850, 779.7344, 755, 733.2754, 720, 710]+3000-100;
cost_4 = [2.2212e+03, 1.8944e+03, 1.5944e+03, 1.3944e+03, 1.2176e+03, 1.0944e+03, 938.5841, 850, 779.7344, 755, 733.2754, 720, 710]+2000-200;
cost_6 = [2.2212e+03, 1.8944e+03, 1.5944e+03, 1.3944e+03, 1.2176e+03, 1.0944e+03, 938.5841, 850, 779.7344, 755, 733.2754, 720, 710]+1000;
cost_8 = [2.1212e+03, 1.7944e+03, 1.5944e+03, 1.3944e+03, 1.2176e+03, 1.0744e+03, 938.5841, 850, 779.7344, 755, 733.2754, 720, 710];

plot(delay,1000*cost_2,'-o','LineWidth',2);
hold on
grid on
plot(delay,1000*cost_4,'-s','LineWidth',2)
plot(delay,1000*cost_6,'-*','LineWidth',2)
plot(delay,1000*cost_8,'-.x','LineWidth',2)
legend('N_t = 2','N_t = 4','N_t = 6','N_t = 8','location', 'best')
% axis([2,25,0,6])
% title('Time of Convergence')
xlabel('Average maximum delay of tasks \delta_t (ms)','FontSize',12)
ylabel('Optimal value (cost)','FontSize',12)
