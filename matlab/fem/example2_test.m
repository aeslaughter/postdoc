% Timestep 1
Khat = [173.208 , 41.5   ,  85.3333, 42.2917;
         41.5   , 87.875 ,  42.2917,  0.0   ;
         85.333 , 42.2917, 6.66668e7 , 3.33334e7;
         42.2917,  0.0   ,  3.33334e7, 6.66668e7];
     
Fhat = [17116.7; 8583.33; 3.0e10; 3.0e10];
Khat\Fhat


Khat = [173.208 , 41.5   ,  0, 0;
         41.5   , 87.875 ,  0,  0.0   ;
        0 , 0, 3.33335e7 ,  1.66667e7;
         0,  0.0   ,  1.66667e7, 3.33334e7];
     
Fhat = [17116.7; 8583.33; 2.75e10; 2.75e10];

Khat\Fhat


Khat = [173.208 , 41.5   ,  0, 0;
         41.5   , 87.875 ,  0,  0.0   ;
         0 , 0, 1 , 0;
         0,  0.0   ,0, 1];
     
Fhat = [49210.6;193438; 300; 300];

Khat\Fhat