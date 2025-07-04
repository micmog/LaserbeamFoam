deltaH = rhoM*(cpL - cpS)*T + rhoM*Hf;

//Implicit Term coefficient

SI =  - deltaH*DFbyDT(T, Tl, Ts)*alphaM / runTime.deltaT();

//Explicit Term - Divergence part included in TEqn

SE = alphaM * deltaH*(gT.oldTime() - gT)/ runTime.deltaT() - SI*invLiquidFraction(gT, Tl, Ts);

