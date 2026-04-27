function qn = quat_normalize(q)
%QUAT_NORMALIZE Quaternion'u normalize eder.

    if ~isvector(q) || numel(q) ~= 4
        error('quat_normalize girdisi 4 elemanlı quaternion olmalıdır.');
    end

    q = q(:);
    nq = norm(q);

    if nq < 1e-12
        error('Quaternion normu çok küçük, normalize edilemez.');
    end

    qn = q / nq;
end