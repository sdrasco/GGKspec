function h=SeriesCartesianMetric(orbit,octupole_flag)
%
% h=SeriesCartesianMetric(orbit)
%
% Series version of CartesianMetric

% See also CARTESIANMETRIC KLUDGEDINSPIRAL KLUDGEDXOFL HPHX
%
% Steve Drasco
%

% start clock to keep track of computational cost
InitialCPUTime = cputime;

% series for x y and z
cosT = SeriesRe(orbit.eiT);
z = SeriesMultiply(orbit.r,cosT);
sinT = SeriesIm(orbit.eiT);
rsinT = SeriesMultiply(orbit.r,sinT);
cosphi = SeriesRe(orbit.eiphi);
sinphi = SeriesIm(orbit.eiphi);
x = SeriesMultiply(rsinT,cosphi);
y = SeriesMultiply(rsinT,sinphi);

% quadrupole moment I^{jk} 
Ixx = SeriesMultiply(x,x);
Iyy = SeriesMultiply(y,y);
Izz = SeriesMultiply(z,z);
Ixy = SeriesMultiply(x,y);
Ixz = SeriesMultiply(x,z);
Iyz = SeriesMultiply(y,z);

% first t-derivatives of I^{jk} 
dIxx = SeriesDifferentiate(Ixx);
dIyy = SeriesDifferentiate(Iyy);
dIzz = SeriesDifferentiate(Izz);
dIxy = SeriesDifferentiate(Ixy);
dIxz = SeriesDifferentiate(Ixz);
dIyz = SeriesDifferentiate(Iyz);

% second t-derivatives of I^{jk} 
ddIxx = SeriesDifferentiate(dIxx);
ddIyy = SeriesDifferentiate(dIyy);
ddIzz = SeriesDifferentiate(dIzz);
ddIxy = SeriesDifferentiate(dIxy);
ddIxz = SeriesDifferentiate(dIxz);
ddIyz = SeriesDifferentiate(dIyz);

% quadrupole term
h.quad.xx = ddIxx;
h.quad.xx.coefs = 2 * h.quad.xx.coefs;
h.quad.yy = ddIyy;
h.quad.yy.coefs = 2 * h.quad.yy.coefs;
h.quad.zz = ddIzz;
h.quad.zz.coefs = 2 * h.quad.zz.coefs;
h.quad.xy = ddIxy;
h.quad.xy.coefs = 2 * h.quad.xy.coefs;
h.quad.xz = ddIxz;
h.quad.xz.coefs = 2 * h.quad.xz.coefs;
h.quad.yz = ddIyz;
h.quad.yz.coefs = 2 * h.quad.yz.coefs;

if nargin > 1 % WARNING: doesn't pass tests
    
    % coordinate derivatives
    dx = SeriesDifferentiate(x);
    dy = SeriesDifferentiate(y);
    dz = SeriesDifferentiate(z);
    
    % second t-derivatives of S^{ijk}
    ddSxxx = SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(dx,Ixx)));
    ddSxyy = SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(dx,Iyy)));
    ddSxzz = SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(dx,Izz)));
    ddSxxy = SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(dx,Ixy)));
    ddSxxz = SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(dx,Ixz)));
    ddSxyz = SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(dx,Iyz)));
    ddSyxx = SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(dy,Ixx)));
    ddSyyy = SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(dy,Iyy)));
    ddSyzz = SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(dy,Izz)));
    ddSyxy = SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(dy,Ixy)));
    ddSyxz = SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(dy,Ixz)));
    ddSyyz = SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(dy,Iyz)));
    ddSzxx = SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(dz,Ixx)));
    ddSzyy = SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(dz,Iyy)));
    ddSzzz = SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(dz,Izz)));
    ddSzxy = SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(dz,Ixy)));
    ddSzxz = SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(dz,Ixz)));
    ddSzyz = SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(dz,Iyz)));
    
    % scale ddSijk by -2
    ddSxxx.coefs = -2 * ddSxxx.coefs;
    ddSxyy.coefs = -2 * ddSxyy.coefs;
    ddSxzz.coefs = -2 * ddSxzz.coefs;
    ddSxxy.coefs = -2 * ddSxxy.coefs;
    ddSxxz.coefs = -2 * ddSxxz.coefs;
    ddSxyz.coefs = -2 * ddSxyz.coefs;
    ddSyxx.coefs = -2 * ddSyxx.coefs;
    ddSyyy.coefs = -2 * ddSyyy.coefs;
    ddSyzz.coefs = -2 * ddSyzz.coefs;
    ddSyxy.coefs = -2 * ddSyxy.coefs;
    ddSyxz.coefs = -2 * ddSyxz.coefs;
    ddSyyz.coefs = -2 * ddSyyz.coefs;
    ddSzxx.coefs = -2 * ddSzxx.coefs;
    ddSzyy.coefs = -2 * ddSzyy.coefs;
    ddSzzz.coefs = -2 * ddSzzz.coefs;
    ddSzxy.coefs = -2 * ddSzxy.coefs;
    ddSzxz.coefs = -2 * ddSzxz.coefs;
    ddSzyz.coefs = -2 * ddSzyz.coefs;
    
    % third t-derivatives of M^{ijk}
    dddMxxx = SeriesDifferentiate(SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(x,Ixx))));
    dddMxyy = SeriesDifferentiate(SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(x,Iyy))));
    dddMxzz = SeriesDifferentiate(SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(x,Izz))));
    dddMxxy = SeriesDifferentiate(SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(x,Ixy))));
    dddMxxz = SeriesDifferentiate(SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(x,Ixz))));
    dddMxyz = SeriesDifferentiate(SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(x,Iyz))));
    dddMyxx = SeriesDifferentiate(SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(y,Ixx))));
    dddMyyy = SeriesDifferentiate(SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(y,Iyy))));
    dddMyzz = SeriesDifferentiate(SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(y,Izz))));
    dddMyxy = SeriesDifferentiate(SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(y,Ixy))));
    dddMyxz = SeriesDifferentiate(SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(y,Ixz))));
    dddMyyz = SeriesDifferentiate(SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(y,Iyz))));
    dddMzxx = SeriesDifferentiate(SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(z,Ixx))));
    dddMzyy = SeriesDifferentiate(SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(z,Iyy))));
    dddMzzz = SeriesDifferentiate(SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(z,Izz))));
    dddMzxy = SeriesDifferentiate(SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(z,Ixy))));
    dddMzxz = SeriesDifferentiate(SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(z,Ixz))));
    dddMzyz = SeriesDifferentiate(SeriesDifferentiate(SeriesDifferentiate(SeriesMultiply(z,Iyz))));

    % octupole term
    h.oct.xxx = SeriesAdd(ddSxxx,dddMxxx);
    h.oct.xyy = SeriesAdd(ddSxyy,dddMxyy);
    h.oct.xzz = SeriesAdd(ddSxzz,dddMxzz);
    h.oct.xxy = SeriesAdd(ddSxxy,dddMxxy);
    h.oct.xxz = SeriesAdd(ddSxxz,dddMxxz);
    h.oct.xyz = SeriesAdd(ddSxyz,dddMxyz);
    h.oct.yxx = SeriesAdd(ddSyxx,dddMyxx);
    h.oct.yyy = SeriesAdd(ddSyyy,dddMyyy);
    h.oct.yzz = SeriesAdd(ddSyzz,dddMyzz);
    h.oct.yxy = SeriesAdd(ddSyxy,dddMyxy);
    h.oct.yxz = SeriesAdd(ddSyxz,dddMyxz);
    h.oct.yyz = SeriesAdd(ddSyyz,dddMyyz);
    h.oct.zxx = SeriesAdd(ddSzxx,dddMzxx);
    h.oct.zyy = SeriesAdd(ddSzyy,dddMzyy);
    h.oct.zzz = SeriesAdd(ddSzzz,dddMzzz);
    h.oct.zxy = SeriesAdd(ddSzxy,dddMzxy);
    h.oct.zxz = SeriesAdd(ddSzxz,dddMzxz);
    h.oct.zyz = SeriesAdd(ddSzyz,dddMzyz);
    
end


% log the computational cost of this job
h.CPUsec = cputime - InitialCPUTime;
