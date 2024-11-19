fx = @(V) V;        %derivada em ordem ao tempo de x
fv = @(X) -k*X/m;

for
    r1x = fx(v(i))
    riv = fv(x(i))

    r2x = fx(v(i)+h/2*r1x)
    r2v = fv(x(i)+h/2*riv)

    r3x = fx(v(i)+0*h*r1x+h/2*r2x)
    r3v = fv(x(i)+0*h*r1v+h/2*r2v)

    r4x = fx(v(i)+0*h*r1x+0*h*r2x+1*h*r3x)
    r4v = fv(x(i)+0*h*r1v+0*h*r2v+1*h*r3v)

    x(i+1) = x(i)+h*(1/6*r1x+1/3*r2x+1/3*r3x+1/6*r4x)
    v(i+1) = v(i)+h*(1/6*r1v+1/3*r2v+1/3*r3v+1/6*r4v)
end