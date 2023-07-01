function roe = MMA(roe,volumeCon,varargin) % function handles passed in, input roe should be col vec

max_k = 5;

% constants
C0 = 0.5; % paper said 0.5 is reasonable value
C1 = 0.7; % paper said 0.7 is reasonable, 0.5 could be more stable
C2 = 1.3; % reasonable

% fmincon vars
lb = zeros(size(roe))+1E-3; % lb and ub will probably be ~0 and 1 because those are max and min densities
ub = ones(size(roe));
A = (100/volumeCon)*ones(size(roe')); % volumeCon is max allowable percentage of total domain vol used
b = sum(ones(size(roe)));
Aeq = [];
beq = [];

nargine = length(varargin);
threshold = 1E-3;

% vars per iteration
L = zeros(size(roe));
U = zeros(size(roe));
past_L = zeros(2,length(roe));
past_U = zeros(2,length(roe));
past_roe = zeros(2,length(roe));

% vars per objs/constrs
r = zeros(nargine,1); % residuals
p = zeros(length(roe),nargine); % each col is vert vector of p's for a given obj/constraint
q = zeros(length(roe),nargine); % each col is vert vector of q's for a given obj/constraint
approx_objs = cell(1,nargine); % this will store the sub-problem obj/constr functions in each col
hessians = cell(1,nargine); % this will store the hessians per obj/constr

for k = 1:max_k % idk what k should end at
    if k < 3 % first 2 iterations
        L = roe-C0*(ub-lb);
        U = roe+C0*(ub-lb);
    else
        idx = logical(max((roe-past_roe(1,:)').*(past_roe(1,:)'-past_roe(2,:)'),0)); % false if entry is neg, true if not
        diff1 = past_roe(1,:) - past_L(1,:); % can't index temp arrays so need to do this
        diff2 = past_U(1,:) - past_roe(1,:);
        
        % if (roe-past_roe(0,0))*(past_roe(0,0)-past_roe(0,1)) > 0
        L(idx) = roe(idx) - C2*diff1(idx)';
        U(idx) = roe(idx) + C2*diff2(idx)';
        
        % if (roe-past_roe(0,0))*(past_roe(0,0)-past_roe(0,1)) <= 0
        L(~idx) = roe(~idx) - C1*diff1(~idx)';
        U(~idx) = roe(~idx) + C1*diff2(~idx)';  
    end
    L(L<lb) = lb(L<lb);
    U(U>ub) = ub(U>ub);
    
    % calc sub-function for each obj and constr and fill r
    for j = 1:nargine
        f = varargin{j}; % whatever function you are on
        [objVal,objDir] = f(roe); % get obj and sensitivity
        idx = logical(objDir > 0); % find which elements have pos sensitivities
        p(idx,j) = (U(idx)-roe(idx)).^2.*objDir(idx); % pos sens --> p is set
        q(~idx,j) = -(roe(~idx)-L(~idx)).^2.*objDir(~idx); % neg sens --> q is set
        r(j) = objVal-sum((p(:,j)./(U-roe)+q(:,j)./(roe-L))); % must pass current roe to f and then subtract function summation from paper
        approx_objs{j} = @(x) MMA_f(x,L,U,p,q,r,j);
    end 
    
    % testing compliance
    approx_objs{1} = @(x) MMA_f(x,L,U,p,q,r,1);
    approx_objs{2} = @(x) MMA_f(x,L,U,p,q,r,2);
    
    
    % store 2 previous iterations of L and U and roe
    past_L(2,:) = past_L(1,:);
    past_L(1,:) = L';
    past_U(2,:) = past_U(1,:);
    past_U(1,:) = U';
    past_roe(2,:) = past_roe(1,:);
    past_roe(1,:) = roe';
    
    % run fmincon on subproblem w/ or w/out nonlinear constraints
    if nargine > 1
        % need to add hessian if there are nonlin constrts
        options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'HessianFcn',@(x,lambda) MMA_hess(x,lambda,L,U,p,q)); % ,'OutputFcn',@outFcn); % for graphics
        roe = fmincon(approx_objs{1},roe,A,b,Aeq,beq,L,U,approx_objs{2:end},options); % with arbitrary num of nlcons
    else 
        options = optimoptions(@fmincon,'SpecifyObjectiveGradient',true,'HessianFcn',@(x,lambda) MMA_hess(x,lambda,L,U,p,q)); % ,'OutputFcn',@outFcn);
        roe = fmincon(approx_objs{1},roe,A,b,Aeq,beq,L,U,[],options); % with no nlcons
    end
    
    % threshold to stop optimization
    if abs(max(roe'-past_roe(1,:))) < threshold
        display(strcat('Final k: ',num2str(k)))
        break;
    end
end
end