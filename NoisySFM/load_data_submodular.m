function [F,param_F] = load_data_submodular(dataset,seed)

% % load data for minimizing submodular function
% rand('state',seed);
% randn('state',seed);

rng(seed)

switch dataset
    
    case 1 % GENRMF LONG (seed in 1-10)
        
        W = read_genrmf(sprintf('graph_long_%d.txt',seed));
        p = size(W,1)-2;
        
        % randomize order
        rp = randperm(p);
        W = W([ 1 rp+1 p+2] ,[ 1 rp+1 p+2]);
        
        F = @submodular_fct_quad_cut;
        param_F.p = p;
        param_F.W = W;
        
        
    case 2 % GENRMF WIDE (seed in 1-10)
        
        
        W = read_genrmf(sprintf('graph_wide_%d.txt',seed));
        p = size(W,1)-2;
        
        % randomize order
        rp = randperm(p);
        W = W([ 1 rp+1 p+2] ,[ 1 rp+1 p+2]);
        
        F = @submodular_fct_quad_cut;
        param_F.p = p;
        param_F.W = W;
        
    case 3 % TWO MOONS
        
        % data parameters
        n1 = 200;
        n2 = 200;
        noise = .4;
        
        
        % build two moons
        r1 = 2 + noise * randn(n1,1);
        theta1 = -pi/2 + pi*rand(n1,1);
        X1 = [  -.5+r1.* cos(theta1), 1+ r1 .* sin(theta1)];
        
        r2 = 2 + noise * randn(n2,1);
        theta2 = pi/2 + pi*rand(n1,1);
        X2 = [  .5 + r2.* cos(theta2), -1+ r2 .* sin(theta2)];
        X = [ X1;X2];
        p = n1+n2;
        
        % randomize order
        rp = randperm(p);
        X = X(rp,:);
        
        
        % compute kernel matrix
        D = sqdist(X',X');
        K = exp(- 8 * D / median(D(:))) + eye(p) * 1e-4 * p ;
        
        
        % define_submodular_function
        F = @submodular_fct_logdetsym;
        param_F.p = p;
        param_F.K = K;
        param_F.FV = sum(log(diag(chol(K))))*2;
        
        % compute minimal value to add to block certain nodes to be in or out
s = zeros(p,1);
param_F.s = s;
for i=1:p, Fsingletons(i) = F(i,param_F); end
val = 2 * max(Fsingletons);

% define labelled points
rp =1:p;
nlabeled = 8;
valp = rp(1:nlabeled);
valn = rp(end-nlabeled+1:end);
s = zeros(p,1);
s(valp)=val;
s(valn)=-val;
param_F.s = s;


% randomize order
rp = randperm(p);
param_F.s = param_F.s(rp);
param_F.K = param_F.K(rp,rp);


        
    case 4 % SET COVER FROM BILMES' GROUP - single one
        
        [W] = read_bilmes('uniform-sent.800');
        
        % randomize order
        rp = randperm(size(W,1));
        W = W(rp,:);
        
        
        F = @submodular_fct_bilmes;
        param_F.p = size(W,1);
        param_F.lambda = 22;
        param_F.alpha = 0.5;
        param_F.W = W;
        
        
    case 5 % IWATA'S TEST FUNCTION
        
        p = 100;
        param_F.p = p;
        F = @submodular_fct_iwata;
        
        
        
    case 6 % SIMPLE IMAGES
        
        
        im = imread('gray3_labels.png');
        im = imresize(im,.5);
        valp = find( im(:,:,1)<100);
        valn = find( im(:,:,3)<200);
        
        im = imread('gray3.png');
        im = imresize(im,.5);
        im = sum(double(im)/255/3,3);
        im = im + randn(size(im))/20;
        
        
        [nx ny] = size(im);
        p = nx*ny;
        [XX,YY] = meshgrid(1:ny,1:nx);
        
        
        % compute kernel matrix
        D1 = sqrt(sqdist([XX(:) YY(:) ]', [XX(:) YY(:) ]' ));
        D2 = sqrt(sqdist(im(:)', im(:)' ));
        K = exp(- 32 * D1 / median(D1(:)) - 4 * D2 / median(D2(:))) + eye(p) * 1e-4 * p ;
        
        
        K = K - diag(diag(K));
        d = sum(K,2);
        F = @submodular_fct_cut;
        
        % randomize order
        rp = randperm(p);
        K = K(rp,rp);
        d = sum(K,2);
        
        
        param_F.W = K;
        param_F.d = d;
        param_F.p = p;
        
        % compute minimal value to add to block certain nodes
s = zeros(p,1);
param_F.s = s;
for i=1:p, Fsingletons(i) = F(i,param_F); end
val = 2 * max(Fsingletons);

% labelled
s = zeros(p,1);
s(valp)=val;
s(valn)=-val;
param_F.s = s;


    case 7
        
            p = 64;
    n = round(log2(p));
    % long
    a = round(2^(n/4));
    b = round(2^(n/2));
    % wide
    % a = round(2^(2*n/5));
    % b = round(2^(n/5));
    c1 = 0;
    c2 = 4;

    unix(sprintf('../genrmf.src/genrmf -a %d -b %d -c1 %d  -c2 %d -seed %d -out ''graph.txt''',a,b,c1,c2,seed));
    [W,tab ] = read_genrmf('graph.txt');
                p = size(W,1)-2;
        
        % randomize order
        rp = randperm(p);
        W = W([ 1 rp+1 p+2] ,[ 1 rp+1 p+2]);
        
        F = @submodular_fct_quad_cut;
        param_F.p = p;
        param_F.W = W;

end