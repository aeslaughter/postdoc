%% Solving dT/dh

syms h x f cf cs Te hf hl
EQ1 = f*((cf - cs)  *(x - Te) + hf) + cs*x - h;
EQ2 = solve(EQ1,x);
EQ3 = diff(EQ2,sym('h'));
pretty(EQ3)


EQ1 = (1-f) * ((cf - cs) * (x - Te) + hf) - hl - h;
EQ2 = solve(EQ1,x);
EQ3 = diff(EQ2,sym('h'));
pretty(EQ3)