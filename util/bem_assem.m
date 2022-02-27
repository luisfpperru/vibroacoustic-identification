function P=bem_assem(edof,P,Pe,n,el)
% P=bem_assem(edof,P,Pe,n,el)
%----------------------------------------------------------------------
% PURPOSE
% Assemble element influence matrix Pe for acoustic problems into
% the global influence matrix P according to the topology matrix
% edof, the influenced node n, and the influencing element el.
%
% INPUT: edof: dof topology matrix
% P: global influence matrix
% Pe: element influence matrix
% n: influenced node
% el: influencing element
%
% OUTPUT: P: New global influence matrix
%----------------------------------------------------------------------
N=size(edof); 
if N(1,2)==2
t=abs(edof(:,1)-el);
[val,p]=min(t);
P(n,edof(p,2))=P(n,edof(p,2))+Pe;
elseif N(1,2)==5
t=abs(edof(:,1)-el);
[val,p]=min(t);
P(n,edof(p,2:5))=P(n,edof(p,2:5))+Pe;
end
%-------------------------------end------------------------------------

% function P=bem_assem(edof,P,Pe,n,el)
% % P=bem_assem(edof,P,Pe,n,el)
% %----------------------------------------------------------------------
% % PURPOSE
% % Assemble element influence matrix Pe for acoustic problems into
% % the global influence matrix P according to the topology matrix
% % edof, the influenced node n, and the influencing element el.
% %
% % INPUT: edof: dof topology matrix
% % P: global influence matrix
% % Pe: element influence matrix
% % n: influenced node
% % el: influencing element
% %
% % OUTPUT: P: New global influence matrix
% %----------------------------------------------------------------------
% N=size(edof); 
% if N(1,2)==2
% t=abs(edof(:,1)-el);
% [val,p]=min(t);
% P(n,edof(p,2))=P(n,edof(p,2))+Pe;
% elseif N(1,2)==5
% t=abs(edof(:,1)-el);
% [val,p]=min(t);
% P(n,edof(p,2:5))=P(n,edof(p,2:5))+Pe;
% end
% %-------------------------------end------------------------------------