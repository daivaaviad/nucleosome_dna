function [U, wopt] = seq2nucleosome_grad_descent(S)

% 
% This function finds an optimal DNA configuration on a nucleosome
% for a given sequence S, by minimising the cgDNA+ model energy
% subject to constraints on the phosphate positions.
% 
% Input:
% 
% S    - any chosen DNA sequence of length 147.
%        For example: S = randseq(147);
% 
% Output:
% 
% wopt - optimal configuration on a nucleosome for the sequence S;
% U    - potential (free) energy, needed to deform DNA sequence S
%        into the configuration wopt. 
% 
% If you find this code useful, please cite:
% 
% Rasa Giniunaite and Daiva Petkeviciute-Gerlach.
% Predicting the sequence-dependent configuration and
% energy of DNA in a nucleosome by coarse-grain modelling.
% Submitted (2022).
% 

addpath('cgDNA+')

% constraints for phosphate positions:
load('refintervalsMinMax', 'ref'); 

% initial configuration in cgDNA+ coordinates x, 
% central axis of the nucleosome aaxis, cpoint
% and the diad point coordinates r0
% plus orientation of the central base pair R0:
load('initial_point', 'aaxis', 'cpoint', 'r0', 'R0', 'x'); 

% cgDNA+ parameter set:
load('cgDNA+/cgDNAplus_paramset_RSnoMH', 'params');


[w, K] = constructSeqParms(S, params);
grad = K*(w-x);

ref1.indc = [];
ref1.indw = [];

ix = nz_entries(x, ref);
ix0 = ones(size(x));
beta = 0.0001;  % step size for the gradient (all coordinates) 
gamma = 0.005;  % step size for the gradient (only inter base pair coordinates) 
alpha = 0.8;    % multiplier to reduce the step size on inter base pair coordinates

U = 0.5*(x-w)'*K*(x-w);

for ii = 1:5000  % maximum number of iterations
  
  xt = x + beta * grad.*ix0 + gamma * grad.*ix;
  
  d = int(xt, aaxis, cpoint, ref, r0, R0);
  
  dallc = d.rc + d.hc + d.ac;
  dallw = d.rw + d.hw + d.aw;
  
  if sum(dallc + dallw)
     if sum(dallc)  
       iviolated = find(dallc);
       for j = 1:length(iviolated) 

         if (isempty(ref1.indc))||(sum(ref1.indc == ref.indc(iviolated(j)))==0)
            ref1.indc = [ref1.indc ref.indc(iviolated(j))];
         end    
         
       end
     end  
     if sum(dallw)
       iviolated = find(dallw);
       for j = 1:length(iviolated)

         if (isempty(ref1.indw))||(sum(ref1.indw == ref.indw(iviolated(j)))==0)
            ref1.indw = [ref1.indw ref.indw(iviolated(j))];
         end  
         
       end  
     end
     ix = nz_entries(x,ref1); 
     
     refall = [ref1.indc ref1.indw];
     ic = 74*24;
     i1 = min(refall(refall<74))-1;
     ix0(i1*24-17:ic) = ix0(i1*24-17:ic)*alpha;
     
     ic = 73*24;
     i2 = max(refall(refall>74))+1;
     ix0(ic+1:i2*24-18) = ix0(ic+1:i2*24-18)*alpha;
     
   else
     if norm(grad - K*(w-xt)) < 0.001 % stop if the algorithm converged
         U = 0.5*(x-w)'*K*(x-w)/147; 
         wopt = x;
         break
     else % iterate further   
       x = xt; 
       grad = K*(w-x);
     end  
   end    

end

end % function

%%%%%%%%%%%%%%%%%%

function ix = nz_entries(v,ref)

    %Selecting entries to be used for gradient descent
  
    [eta,w,etapW,wpW,u,v,etapC,wpC] = vector2shapes(v);

    u = zeros(size(u));
    v = zeros(size(v));

    etapW(ref.indw-1,:) = zeros(length(ref.indw),3);
    wpW(ref.indw-1,:) = zeros(length(ref.indw),3); 
    etapC(ref.indc,:) = zeros(length(ref.indc),3);
    wpC(ref.indc,:) = zeros(length(ref.indc),3);
    
    eta(ref.indw,:) = zeros(length(ref.indw),3);
    w(ref.indw,:) = zeros(length(ref.indw),3); 
    eta(ref.indc,:) = zeros(length(ref.indc),3);
    w(ref.indc,:) = zeros(length(ref.indc),3); 
  
    g = shapes2vector(eta,w,etapW,wpW,u,v,etapC,wpC);
    ix = abs(g)>0;
   
end

function d = int(x, aaxis, cpoint, ref, r0, R0)
  
  % Checking if the constrains are satisfied
  
  % ----- Absolute coordinates
  
  bp_level = frames_ref(x, ref, r0, R0);
  
  rpoint = r0; 
  rpaxis = cpoint + dot(rpoint-cpoint, aaxis) * aaxis;
  rvector = rpoint - rpaxis;
  eps = 0.01;
  
  % ----- Cylindrical coordinates of phosphates 

  % Crick strand
  for i = 1:length(ref.indc) %only touching phosphates

    p = bp_level(ref.indc(i)).rpc;
    rpc = norm(cross(cpoint - p, aaxis));  % d coordinate 
    hpc = abs(dot(aaxis,rpoint - p) + 40); % h coordinate   
    
    pvector = p - cpoint - dot(p - cpoint, aaxis) * aaxis;
    
    v12 = cross(rvector,pvector);
    c = sign(dot(v12, aaxis)) * norm(v12);
    a2pi = atan2(c, dot(rvector,pvector));
    if (ref.indc(i) < 50)&&(a2pi > 0) 
       a2pi = a2pi - 2*pi;
    elseif (ref.indc(i) > 100)&&(a2pi < 0)
       a2pi = a2pi + 2*pi;  
    end  
    apc = a2pi; % phi coordinate
    
    d.rc(i) = ~((rpc - ref.rc(1,i) > -eps)&&(ref.rc(2,i) - rpc > -eps));
    d.hc(i) = ~((hpc - ref.hc(1,i) > -eps)&&(ref.hc(2,i) - hpc > -eps));
    d.ac(i) = ~((apc - ref.ac(1,i) > -eps)&&(ref.ac(2,i) - apc > -eps));
  end    
   
  % Watson strand
  for i = 1:length(ref.indw) %only touching phosphates
     
    p = bp_level(ref.indw(i)).rpw;
    rpw = norm(cross(cpoint - p, aaxis));   % d coordinate 
    hpw = abs(dot(aaxis, rpoint - p) + 40); % h coordinate 
    
    pvector = p - cpoint - dot(p - cpoint, aaxis) * aaxis;
    
    v12 = cross(rvector,pvector);
    c = sign(dot(v12, aaxis)) * norm(v12);
    a2pi = atan2(c, dot(rvector,pvector));
    if (ref.indw(i) < 50)&&(a2pi > 0) 
       a2pi = a2pi - 2*pi;
    elseif (ref.indw(i) > 100)&&(a2pi < 0)
       a2pi = a2pi + 2*pi;  
    end  
    apw = a2pi; % phi coordinate
  
    d.rw(i) = ~((rpw - ref.rw(1,i) > -eps)&&(ref.rw(2,i) - rpw > -eps));
    d.hw(i) = ~((hpw - ref.hw(1,i) > -eps)&&(ref.hw(2,i) - hpw > -eps));
    d.aw(i) = ~((apw - ref.aw(1,i) > -eps)&&(ref.aw(2,i) - apw > -eps));
  end   

end

function bp_level = frames_ref(shapes, ref, r0, R0)  

    %reconstruction starting from the middle base-pair (the dyad)
    
    nbp = 147;
    
    % absolute coordinates of the middle base-pair 
    G = R0;
    q = r0;
   
    % relative coordinates of the oligomer
    [eta, w, ~, wpW, u, v, ~, wpC] = vector2shapes(shapes);
    
    bp_level = InitializeStruct(nbp) ;
    
    for i = 74:nbp %second half (forward reconstruction)
        
        %bp_level(i).r = q ;
        
        if sum(i == [ref.indc ref.indw]) %only touching phosphates
            
          % base pair:
          r = cay(eta(i,:));
          Gw = G * w(i,:)';  
 
          % complementary strand 
          bp_level(i).Rc = G * (sqrtm(r))'; 
          bp_level(i).rc = q - 0.5 * Gw;
        
          % main strand
          bp_level(i).Rw = bp_level(i).Rc * r;
          bp_level(i).rw = bp_level(i).rc + Gw;
          
        end 
          
        if i < nbp
            ru = cay(u(i,:));
            sqrtru = sqrtm(ru);
            H = G * sqrtru;
            % next base pair:
            G = G * ru;
            q = q + H * v(i,:)';   
        end  
    end
    
    G = R0;
    q = r0;
       
    for i = 74:-1:1 %first half (backward reconstruction)
        
        %bp_level(i).r = q ;
        
        if sum(i == [ref.indc ref.indw]) %only touching phosphates
            
          % base pair:
          r = cay(eta(i,:));
          Gw = G * w(i,:)';  
 
          % complementary strand 
          bp_level(i).Rc = G * (sqrtm(r))'; 
          bp_level(i).rc = q - 0.5 * Gw;
        
          % main strand
          bp_level(i).Rw = bp_level(i).Rc * r;
          bp_level(i).rw = bp_level(i).rc + Gw;
          
        end 
          
        if i > 1
            ru = cay(-u(i-1,:));
            sqrtru = sqrtm(ru);
            H = G * sqrtru;
            % next base pair:
            G = G * ru;
            q = q - H * v(i-1,:)';   
        end  
    end
    
    for i = 1:length(ref.indw)
      ii = ref.indw(i);  
      bp_level(ii).rpw  = bp_level(ii).rw + bp_level(ii).Rw*wpW(ii-1,:)'; 
    end  
      
    for i = 1:length(ref.indc)
      ii = ref.indc(i);  
      Rc = bp_level(ii).Rc*diag([1,-1,-1]) ;
      bp_level(ii).rpc  = bp_level(ii).rc + Rc*wpC(ii,:)';
    end  
    
end

function bp_level = InitializeStruct(nbp)

  bp_level  = struct('R', [],'r',[], ...
             'Rw', [], 'rw', [], ...
             'Rc', [], 'rc', [], ...
             'Rpw', [], 'rpw', [], ...
             'Rpc', [], 'rpc', cell(nbp,1)) ;
end

function Q = cay(k)

    I = eye(3) ;
    alpha = 1/10 ;
    k = alpha*k ;
    X = [   0   -k(3)  k(2) ;
           k(3)   0   -k(1) ;
          -k(2)  k(1)   0 ] ;
    Q = (I+X)/(I-X) ;

end


%-------------------------------------------------------------
