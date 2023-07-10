% [t1,error1,amp1,u1,uexact1] = FTBS(2,50);
[t2,error2,amp2,u2,uexact2] = FTBS(1,50);
% [t3,error3,amp3,u3,uexact3] = FTBS(0.5,50);
% [t4,error4,amp4,u4,uexact4] = FTBS(4,50);
% [t5,error5,amp5,u5,uexact5] = FTBS(8,50);
% [t6,error6,amp6,u6,uexact6] = FTBS(2,10);
% [t7,error7,amp7,u7,uexact7] = FTBS(1,10);
% [t8,error8,amp8,u8,uexact8] = FTBS(2,50,1,1);
% x = linspace(0,50,51);
waveplot(u2,uexact2)