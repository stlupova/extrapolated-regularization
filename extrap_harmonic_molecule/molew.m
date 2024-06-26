% this program reproduces results in Fig. 9, left
%   in``Extrapolated regularization...''
%   for the molecular surface, 5th order reg'n, delta=rho*h
%   values are in Table 10 of original version on arxiv

tic;
tstart = tic;

N = input('enter 1/h     ');  % assume N even  !!
h = 1/N;
w = 1 + 10*h;  % computational box is square, -w to w
npt = 2*N + 19;  % npt = number of interior points
theta = 70*pi/180;  % a parameter for quadrature
below = .9*cos(theta); % to reject pts not needed for quad
fourpi = 4*pi; rtpi = sqrt(pi);
rt2 = sqrt(2); rt3 = sqrt(3); rt6 = sqrt(6);
floor = .45; % parameter for quadrature

% set up grid, decide inside or outside
X = zeros(npt,npt,npt); Y = X; Z = X;
for i = 1:npt
for j = 1:npt
for k = 1:npt
X(i,j,k) = i*h-w; Y(i,j,k) = j*h-w; Z(i,j,k) = k*h-w;
end 
end
end

% functions defining molecular surface
%  the centers are xc1...xc4
rc = .5; rcsq = rc^2;
xc1 = [rt3/3; 0; -rt6/12];
xc2 = [-rt3/6; .5; -rt6/12];
xc3 = [-rt3/6; -.5; -rt6/12];
xc4 = [0; 0; rt6/4];
cc1 = [xc1(1) xc2(1) xc3(1) xc4(1)];
cc2 = [xc1(2) xc2(2) xc3(2) xc4(2)];
cc3 = [xc1(3) xc2(3) xc3(3) xc4(3)];
e1 = @(xr,yr,zr) exp(-((xr-xc1(1)).^2+(yr-xc1(2)).^2+(zr-xc1(3)).^2)/rcsq);
e2 = @(xr,yr,zr) exp(-((xr-xc2(1)).^2+(yr-xc2(2)).^2+(zr-xc2(3)).^2)/rcsq);
e3 = @(xr,yr,zr) exp(-((xr-xc3(1)).^2+(yr-xc3(2)).^2+(zr-xc3(3)).^2)/rcsq);
e4 = @(xr,yr,zr) exp(-((xr-xc4(1)).^2+(yr-xc4(2)).^2+(zr-xc4(3)).^2)/rcsq);
molefcn = @(xr,yr,zr) e1(xr,yr,zr)+e2(xr,yr,zr)+e3(xr,yr,zr)+e4(xr,yr,zr) - .6;
lev = molefcn(X,Y,Z);  % the level set function
chi = zeros(npt,npt,npt);   % chi = 1 inside, else = 0
chi(find(lev > 0)) = 1; 
clear lev;

% find irregular intervals, intersecting surface
cross1 = zeros(npt,npt,npt); cross2 = cross1; cross3 = cross1;
for i = 1:npt-1
cross1(i,:,:) = chi(i+1,:,:) - chi(i,:,:);
end
for j = 1:npt-1
cross2(:,j,:) = chi(:,j+1,:) - chi(:,j,:);
end
for k = 1:npt-1
cross3(:,:,k) = chi(:,:,k+1) - chi(:,:,k);
end
% cross1 = -1  if (i,j,k) inside and (i+1,j,k) outside
% cross1 =  1 if (i,j,k) outside and (i+1,j,k) inside
clear chi

list1 = find(cross1); list2 = find(cross2); list3 = find(cross3);
ints1 = length(list1); ints2 = length(list2); ints3 = length(list3);
ints = ints1 + ints2 + ints3; ints12 = ints1 + ints2;

% list endpoints of irregular intervals, eg (ind1,jnd1,knd1)
ind1 = zeros(ints1,1); jnd1 = ind1; knd1 = ind1;
ind2 = zeros(ints2,1); jnd2 = ind2; knd2 = ind2;
ind3 = zeros(ints3,1); jnd3 = ind3; knd3 = ind3;
ii = 0; jj = 0; kk = 0;
for i = 1:npt-1
    for j = 1:npt-1
        for k = 1:npt-1
            if cross1(i,j,k)~=0
            ii = ii + 1;
            ind1(ii) = i; jnd1(ii) = j; knd1(ii) = k;        
            end % if 
            if cross2(i,j,k)~=0
            jj = jj + 1;
            ind2(jj) = i; jnd2(jj) = j; knd2(jj) = k;
            end % if
            if cross3(i,j,k)~=0
            kk = kk + 1;
            ind3(kk) = i; jnd3(kk) = j; knd3(kk) = k;
            end % if
        end % k
    end % j
end % i

% find intersection points (xq,yq,zq) to be quadrature pts
%    and tag grid points that might be within h of surface

itag = zeros(npt,npt,npt);
iq = 0; xq = zeros(ints1,1); yq = zeros(ints2,1); zq = zeros(ints3,1);
eps = 1.e-12;

for k = 1:ints1
ii = ind1(k); jj = jnd1(k); kk = knd1(k);
xold = ii*h - w;  yin = jj*h - w; zin = kk*h - w;
d1 = molefcn(xold+h,yin,zin) - molefcn(xold-h,yin,zin);
d2 = molefcn(xold,yin+h,zin) - molefcn(xold,yin-h,zin);
d3 = molefcn(xold,yin,zin+h) - molefcn(xold,yin,zin-h);
prenl = abs(d1/sqrt(d1^2+d2^2+d3^2)); % approx comp't 1 of normal
itag(ii,jj,kk) = 1; itag(ii+1,jj,kk) = 1;
if prenl < below
continue; end
fcn = @(xin) molefcn(xin,yin,zin);
xcut = fzero(fcn,[xold-eps,xold+h+eps]);
iq = iq + 1;
xq(iq) = xcut; yq(iq) = yin; zq(iq) = zin;
if prenl < floor, continue; end
itag(ii-1,jj,kk) = 1; itag(ii+2,jj,kk) = 1;
end %k
quads1 = iq;

for k = 1:ints2
ii = ind2(k); jj = jnd2(k); kk = knd2(k);
yold = jj*h - w; xin = ii*h - w; zin = kk*h - w;
d1 = molefcn(xin+h,yold,zin) - molefcn(xin-h,yold,zin);
d2 = molefcn(xin,yold+h,zin) - molefcn(xin,yold-h,zin);
d3 = molefcn(xin,yold,zin+h) - molefcn(xin,yold,zin-h);
prenl = abs(d2/sqrt(d1^2+d2^2+d3^2));
itag(ii,jj,kk) = 1; itag(ii,jj+1,kk) = 1;
if prenl < below
continue; end
fcn = @(yin) molefcn(xin,yin,zin);
ycut = fzero(fcn,[yold-eps,yold+h+eps]);
iq = iq + 1;
xq(iq) = xin; yq(iq) = ycut; zq(iq) = zin;
if prenl < floor, continue; end
itag(ii,jj-1,kk) = 1; itag(ii,jj+2,kk) = 1;
end %k
quads12 = iq;

for k = 1:ints3
ii = ind3(k); jj = jnd3(k); kk = knd3(k);
zold = kk*h - w; xin = ii*h - w; yin = jj*h - w;
d1 = molefcn(xin+h,yin,zold) - molefcn(xin-h,yin,zold);
d2 = molefcn(xin,yin+h,zold) - molefcn(xin,yin-h,zold);
d3 = molefcn(xin,yin,zold+h) - molefcn(xin,yin,zold-h);
prenl = abs(d3/sqrt(d1^2+d2^2+d3^2));
itag(ii,jj,kk) = 1; itag(ii,jj,kk+1) = 1;
if prenl < below
continue; end
fcn = @(zin) molefcn(xin,yin,zin);
zcut = fzero(fcn,[zold-eps,zold+h+eps]);
iq = iq + 1;
xq(iq) = xin; yq(iq) = yin; zq(iq) = zcut;
if prenl < floor, continue; end
itag(ii,jj,kk-1) = 1; itag(ii,jj,kk+2) = 1;
end %k
quads = iq;

clear cross1 cross2 cross3

% find normal vector at quadrature pts
nx = zeros(quads,1); ny = nx; nz = nx;
for k = 1:quads
xr = xq(k); yr = yq(k); zr = zq(k);
ee = [e1(xr,yr,zr) e2(xr,yr,zr) e3(xr,yr,zr) e4(xr,yr,zr)];
d1 = xr - cc1; d2 = yr - cc2; d3 = zr - cc3;
pre1 = d1*ee'; pre2 = d2*ee'; pre3 = d3*ee';  % 2 minus signs
mag = sqrt(pre1^2 + pre2^2 + pre3^2);
nx(k) = pre1/mag; ny(k) = pre2/mag; nz(k) = pre3/mag;
end %k

% find target points yy, closest pts zz, signed distance bee
preouts = length(find(itag));
yy1 = zeros(preouts,1); yy2 = yy1; yy3 = yy1;
zz1 = yy1; zz2 = yy1; zz3 = yy1; bee = yy1;
eps = 1.e-12;
iout = 0;
for i = 1:npt-1
for j = 1:npt-1
for k = 1:npt-1
if itag(i,j,k)==0 continue; end %
yout1 = X(i,j,k); yout2 = Y(i,j,k); yout3 = Z(i,j,k);
yr1 = yout1; yr2 = yout2; yr3 = yout3;
% restrict to the first octant
if yr3 < 0, continue; end
if yr2 < 0, continue; end
if yr1 < 0, continue; end
% find dist to surface and closest pt using
%   Lagrange multipliers and Newton's method
%   phi = level set fcn, lam = Lagrange mult'r
%   bigf = (phi,grad(dist^2/2 + lam*grad phi)
%   find dbigf, solve bigf = 0 using Newton
xon1 = yr1; xon2 = yr2; xon3 = yr3; lam = 0;
newterr = 1; it = 0;
while newterr > eps
    it = it + 1;
    if it>20
        disp (' stuck')
        break; end
phi = molefcn(xon1,xon2,xon3);
ee = [e1(xon1,xon2,xon3) e2(xon1,xon2,xon3)...
      e3(xon1,xon2,xon3) e4(xon1,xon2,xon3) ] ;
d1 = xon1 - cc1; d2 = xon2 - cc2; d3 = xon3 - cc3;
d1 = -2*d1/rcsq; d2 = -2*d2/rcsq; d3 = -2*d3/rcsq;
phi1 = d1*ee'; phi2 = d2*ee'; phi3 = d3*ee';
bigf = [phi xon1-yr1+lam*phi1 ...
   xon2-yr2+lam*phi2 xon3-yr3+lam*phi3]';
phi12 = lam*(d1.*d2)*ee'; %really phi12 times lam !!!!
phi23 = lam*(d2.*d3)*ee';
phi13 = lam*(d1.*d3)*ee';
phi11 = lam*(d1.*d1)*ee';
phi22 = lam*(d2.*d2)*ee';
phi33 = lam*(d3.*d3)*ee';
eep = lam*[2 2 2 2]*ee';
phi11 = phi11 - eep; phi22 = phi22 - eep; 
phi33 = phi33 - eep;
dbigf = [0 phi1 phi2 phi3; ...
  phi1 1+phi11 phi12 phi13; ...
  phi2 phi12 1+phi22 phi23; ...
  phi3 phi13 phi23 1+phi33 ] ;
change = - dbigf\bigf;  %Newton!
lam = lam + change(1);
xon1 = xon1 + change(2);
xon2 = xon2 + change(3);
xon3 = xon3 + change(4);
newterr = norm(change);
end %while
gradmag = sqrt(phi1^2 + phi2^2 + phi3^2);
height = lam*gradmag;
if abs(height) > h, continue; end
iout = iout + 1;
yy1(iout) = yout1; yy2(iout) = yout2; yy3(iout) = yout3;
zz1(iout) = xon1; zz2(iout) = xon2; zz3(iout) = xon3;
bee(iout) = - height; %minus!!
end %k
end %j
end %i
outs = iout;

clear X Y Z

% weights for quadrature
apu = 2;
id = zeros(quads,1);
id(1:quads1) = ones;
id(quads1+1:quads12) = 2*ones;
id(quads12+1:quads) = 3*ones;
qw = zeros(quads,1);
for k = 1:quads
idir = id(k);
nn(1) = abs(nx(k)); nn(2) = abs(ny(k)); nn(3) = abs(nz(k));
if nn(idir) < below, continue; end
ww = 0; bb = zeros(3,1);
for i = 1:3
ww = (acos(nn(i))/theta)^2;
if ww < 1, bb(i) = exp(apu*ww/(ww-1)); end
end % i
bbtop = bb(idir);
area = 1/nn(idir);
psi = bbtop/(bb(1) + bb(2) + bb(3));
qw(k) = psi*area*h^2;
end %k


ee = exp(zq);
ff = (sin(xq) + sin(yq)).*ee;
gg = cos(xq).*ee.*nx + cos(yq).*ee.*ny + ff.*nz;

rho = [2 3 4]; del = rho*h;
% h45 = h^.8; del = rho*h45;

sgrand = zeros(quads,3);  % single layer integrand
sh = zeros(quads,3); dsh = sh;
sglint = zeros(1,3); dblint = sglint;
g0 = -(1/2)*pi^(-3/2)*del.^(-1); 
g0 = g0*(-fourpi);  % reg'zed G at 0

uexact = (sin(yy1) + sin(yy2)).*exp(yy3);
chi = ones(outs,1);  %  reusing name chi!!
ip = find(bee>0);
uexact(ip) = 0;
chi(ip) = 0;
ffz = (sin(zz1) + sin(zz2)).*exp(zz3);

finerr = zeros(outs,1); unew = finerr; vnew = finerr;

% compute the integrals with 3 choices of delta
for k = 1:outs
diff1 = xq - yy1(k);
diff2 = yq - yy2(k);
diff3 = zq - yy3(k);
rr = sqrt(diff1.^2 + diff2.^2 + diff3.^2);
temp = diff1.*nx + diff2.*ny + diff3.*nz;
rel = rr./del;
sh = erf(rel);
dsh = sh - (2/rtpi)*rel.*exp(-rel.^2);
sgrand = g0.*ones(quads,1);
non = rr > 1.e-14;
sgrand(non,:) = sh(non,:)./rr(non);
sgrand = sgrand.*gg;
temp = (ff - ffz(k)).*temp;
dgrand = zeros(quads,3);
dgrand(non,:) = temp(non).*dsh(non,:)./rr(non).^3;
for ido = 1:3 
sglint(ido) = - sum(sgrand(:,ido).*qw)/fourpi;
dblint(ido) = sum(dgrand(:,ido).*qw)/fourpi;
end %ido
ucomp = sglint;
vcomp = dblint + chi(k)*ffz(k);

% the extrapolation
lam = bee(k)./del; alam = abs(lam);
eye0 = - alam.*erfc(alam) + exp(-lam.^2)/rtpi;
eye2 = - (lam.^2-1/2).*exp(-lam.^2)/rtpi;
eye2 = (2/3)*(eye2 + alam.^3.*erfc(alam) );
co1 = eye0.*rho;
co2 = eye2.*rho.^3;

themat = [1 co1(1) co2(1); 1 co1(2) co2(2); 1 co1(3) co2(3)];
soln = themat\ucomp';
unew(k) = soln(1);
soln = themat\vcomp';
vnew(k) = soln(1);
finerr(k) = - unew(k) + vnew(k) - uexact(k);

end %k

l2err = norm(finerr)/sqrt(outs)
maxerr = max(abs(finerr))

uex2 = norm(uexact)/sqrt(outs);
uexmax = max(abs(uexact));

% outt = [err2; errmax; uex2; uexmax];


telapsed = toc(tstart)

