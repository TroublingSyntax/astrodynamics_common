function Theta = compTheta(x,mu)
Theta(:,1)=x(:,4);
Theta(:,2)=x(:,5);
Theta(:,3)=x(:,6);
Theta(:,4)=mu./(norm(x(:,1:3))^3);
Theta(:,5)=mu./(norm(x(:,1:3))^3);
Theta(:,6)=mu./(norm(x(:,1:3))^3);
end