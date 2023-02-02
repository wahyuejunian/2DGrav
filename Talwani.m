% Talwani polygon
function g=Talwani(x0,z0,x,z,n_sudut)

temp = 0;
for i = 1:n_sudut
    if i == n_sudut
        n2 = 1;
    else
        n2 = i+1;
    end
    % susut benda 4 sisi
    x1 = x(i)-x0;
    z1 = z(i)-z0;
    x2 = x(n2)-x0;
    z2 = z(n2)-z0;
	
    % r1 and r2
    r1 = x1*x1 + z1*z1;
    r2 = x2*x2 + z2*z2;
    
    dz = z2 - z1;
    if dz == 0
        dz = 1e-6;
    end
   
    an = (x2-x1)/dz;
    bn = (x1*z2-x2*z1)/dz;
	
   
    temp = temp + bn/(1+an^2)*((0.5*(log(r2)-log(r1)))-an*(atan2(z2,x2)-atan2(z1,x1)));
end

g = temp*1e+5;

end
