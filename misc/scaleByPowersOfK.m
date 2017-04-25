function [xout] = scaleByPowersOfK(xin,K)

% if size(xin,1) > 1
%     disp('need 1xN array')
%     return
% end
% 
[T,M] = size(xin);

Kmat = repmat(K.^([1:M]),T,1);

xout = xin./Kmat;


end