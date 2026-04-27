function S = skew3(a)
%SKEW3 3x1 vektörden skew-symmetric matris üretir.

    if ~isvector(a) || numel(a) ~= 3
        error('skew3 girdisi 3 elemanlı bir vektör olmalıdır.');
    end

    a = a(:);

    ax = a(1);
    ay = a(2);
    az = a(3);

    S = [  0   -az   ay;
          az     0  -ax;
         -ay    ax    0 ];
end