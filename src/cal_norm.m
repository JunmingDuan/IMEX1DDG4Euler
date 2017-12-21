function norm = cal_norm(u, w, h, N)
% calculate the N norm of u
% para: u,  val
%       w,  weights
%       h,  length of each intervals
%       N,  L^N norm

norm = 0;
if(N == 2)
  norm = sum(u.^2 .* w * h/2);
  norm = sqrt(norm);
elseif(N == 1)
  norm = sum(abs(u) .* w * h/2);
elseif(N == 'inf')
  for p = 1:length(u)
    if(abs(u(p)) > norm) norm = abs(u(p));
    end
  end
end

