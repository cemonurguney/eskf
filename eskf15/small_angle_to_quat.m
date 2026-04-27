function dq = small_angle_to_quat(dtheta)
%SMALL_ANGLE_TO_QUAT Küçük açı vektöründen küçük düzeltme quaternion'u üretir.
%
% dtheta: 3x1 küçük açı vektörü [rad]
% dq ~= [1; 0.5*dtheta]

    if ~isvector(dtheta) || numel(dtheta) ~= 3
        error('small_angle_to_quat girdisi 3 elemanlı vektör olmalıdır.');
    end

    dtheta = dtheta(:);

    dq = [1;
          0.5 * dtheta];

    dq = quat_normalize(dq);
end