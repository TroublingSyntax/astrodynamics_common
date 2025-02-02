function v_s = EquitoSyn(v_r,n,v_t)
%RtoP This function transforms a vector from the
%     equitorial frame (x,y,z) to the synodic frame for the cr3bp.

v_s = zeros(3,length(v_r));

for i=1:length(v_r)
    F = [cos(n*v_t(i)) sin(n*v_t(i)) 0;
        -sin(n*v_t(i)) cos(n*v_t(i)) 0;
        0              0             1];
    v_s(:,i) = F * v_r(i,:)';
end

v_s = v_s';
end