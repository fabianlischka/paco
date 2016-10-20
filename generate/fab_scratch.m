% to generate the velocity field graph as in Felix, Fig. 1
x1 = 0:(1.0/16):1
x2 = x1
[X1,X2] = meshgrid(x1,x2)
b1 = sin(pi*X1).*sin(pi*X2).*(X2-0.5)X`
b2 = sin(pi*X1).*sin(pi*X2).*(0.5-X1)
quiver(X1,X2,b1,b2)

delb1delx1 = pi*cos(pi*X1).*sin(pi*X2).*(X2-0.5)
delb2delx2 = pi*sin(pi*X1).*cos(pi*X2).*(0.5-X1)
divb = delb1delx1 + delb2delx2
scatter(X1(:), X2(:), 50*(divb(:)+0.6))

divb2 = divergence(X1,X2,b1,b2)
scatter(X1(:), X2(:), 50*(divb(:)+0.6))

divm = pi*(1/2-X1).* sin(pi*X1).*cos(pi*X2) + pi*(X2-1/2).*cos(pi*X1).*sin(pi*X2)
scatter(X1(:), X2(:), 50*(abs(divb(:))+0.01), sign(divb(:)))
