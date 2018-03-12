syms x;
y0 = x^2 * exp(-3 * x^2) + (x / 40)^2;
y1 = diff(y0);
y2 = diff(y1);
y0 = inline(y0);
y1 = inline(y1);
y2 = inline(y2);
