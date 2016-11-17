function [ii,jj,ss] = fvm(vv,q,u,where)

ii = zeros(12,1);
jj = zeros(12,1);
ss = zeros(12,1);
count = 0;
for i=1:3,
    ip = mod(i,3) + 1; ipp = mod(ip, 3) + 1;
    unL = -((q(ip,2) + q(i,2) - 2*q(ipp,2))*u(1)...
        -(q(ip,1) + q(i,1) - 2*q(ipp,1))*u(2))/6;
    ii(count+1) = vv(i);  ss(count+1) = unL;
    ii(count+2) = vv(ip); ss(count+2) = -unL;
    if (unL > 0),
        jj(count+1) = vv(i);
        jj(count+2) = vv(i);
    else
        jj(count+1) = vv(ip);
        jj(count+2) = vv(ip);
    end;
    count = count+2;
    if (where(i)~=0 & where(ip)~=0),
        unL = ((q(ip,2) - q(i,2))*u(1) - (q(ip,1)-q(i,1))*u(2))/2;
        if (unL > 0),
            ii(count+1) = vv(i); jj(count+1) = vv(i); ss(count+1) = unL;
            ii(count+2) = vv(i); jj(count+2) = vv(i); ss(count+2) = unL;
            count = count+2;
        end;
    end;
end;
ii = ii(1:count); jj = jj(1:count); ss = ss(1:count);