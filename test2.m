x = linspace(0,5,100);
y = erfc(x);
plot(x,y);
xlabel('x');
ylabel('erfc(x)');
title('Complementary Error Function');
