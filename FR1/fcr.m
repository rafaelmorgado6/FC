function F = fcr(xv,xold,vold,const)
% const(1), const(2) e const(3) estão definidas no programa principal.
% xold é x(k) e vold é vx(k).
% xv(1) é x(k+1) e xv(2) é vx(k+1).
F(1)= xv(1)-xold-const(1)*(xv(2)+vold);
F(2)= xv(2)-vold+const(2)*(xold+xv(1)+const(3)*(xold^3+xv(1)^3));
end