function dxdt = two_bd_cartesian_model(~,x,xi_sim)
r=x(1:3);
rmag=norm(r);
v=x(4:6);
dxdt=[v(1)*xi_sim(1);v(2)*xi_sim(2);v(3)*xi_sim(3);
      (r(1)*rmag^-3)*xi_sim(4);(r(2)*rmag^-3)*xi_sim(5);...
      (r(3)*rmag^-3)*xi_sim(6)];
end