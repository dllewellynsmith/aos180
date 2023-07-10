%%
plot(x,uLF5(1,:),'LineWidth',2)
legend('Initial Condition')
title('Experiment 4 Initial Condition')
xlabel('x(m)')
ylabel('u')
%%
plot(x,uexactLF5(201,:),'LineWidth',2)
hold on
plot(x,uLF5(201,:),'MarkerSize',12)
plot(x,uRG5(201,:),'LineWidth',2)
legend('Exact Solution','Leapfrog','Runge-Kutta')
title('Experiment 4 Final State')
xlabel('x(m)')
ylabel('u')
%%
semilogy(ampLF5)
hold on
semilogy(ampRG5)
%plot(exp(-4/900*(1:1001)))
legend('Leapfrog','RG3','Fitted Exponential Function')
title('Experiment 4 Amplitude')
xlabel('t (s)')
ylabel('Amplitude')
%%
loglog(tRG5,errorRG5)
hold on
loglog(tLF5,errorLF5)
legend('RG3','Leapfrog')
title('Experiment 4 Error')
xlabel('t (s)')
ylabel('\epsilon')