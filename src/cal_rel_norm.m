function norm = cal_rel_norm(u, v, w, h, N)
% calculate the N norm of u
% para: u,  val
%       v,  exact
%       w,  weights
%       h,  length of each intervals
%       N,  L^N norm

norm = 0;
if(N == 2)
  norm1 = sum((u-v).^2 .* w * h/2);
  norm1 = sqrt(norm1);
  norm2 = sum(v.^2 .* w * h/2);
  norm2 = sqrt(norm2);
  norm = norm1/norm2;
elseif(N == 1)
  norm1 = sum(abs(u-v) .* w * h/2);
  norm2 = sum(abs(v) .* w * h/2);
  norm = norm1/norm2;
elseif(N == 'inf')
  for p = 1:length(u)
    if(abs(v(p)) > 1e-14)
      norm1 = abs((u(p)-v(p))/v(p));
    else
      norm1 = 0;
    end
    if(norm1 > norm) norm = norm1;
    end
  end
end

