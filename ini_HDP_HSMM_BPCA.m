function [struct,diminfo] = ini_HDP_HSMM_BPCA(data,L,seed,varargin)
% data:T\times q
    vars=["alpha_Gamma_a","alpha_Gamma_b","gamma_Gamma_a","gamma_Gamma_b","theta_alpha_threshold","du0","D"];
    values=[10,0.01,10,0.01,0.01,10,15];
    p=inputParser;
    for i=1:length(vars)
        addParameter(p,vars(i),values(i));
    end
    parse(p,varargin{:});
    hyperparameters.alpha_Gamma_a=p.Results.alpha_Gamma_a;
    hyperparameters.alpha_Gamma_b=p.Results.alpha_Gamma_b;
    hyperparameters.gamma_Gamma_a=p.Results.gamma_Gamma_a;
    hyperparameters.gamma_Gamma_b=p.Results.gamma_Gamma_b;
    hyperparameters.du0=p.Results.du0;
    hyperparameters.D=p.Results.D;
    [T,q]=size(data);
    v=q-1;
    diminfo.T=T;
    diminfo.q=q;
    diminfo.v=v;
    diminfo.data=data;
    diminfo.inihyperparameters=hyperparameters;
    diminfo.L=L;
    diminfo.seed=seed;
    diminfo.D=p.Results.D;
    diminfo.du0=p.Results.du0;
    rng(seed)
    datailen=floor(T/L);
%     meandata=mean(data);
%     covdata=cov(data);
    diminfo.inithetaparameters.beta=0.001;
    diminfo.inithetaparameters.alpha_a=0.00001;
    diminfo.inithetaparameters.alpha_b=0.00001*ones(v,1);
    diminfo.inithetaparameters.tau_c=0.001;
    diminfo.inithetaparameters.tau_d=0.001;
    diminfo.inithetaparameters.theta_alpha_threshold=p.Results.theta_alpha_threshold;

%     [~,GammaGMM,~] = main_EM_GMM(data,L,100,1e-6);
%     [~,GMMstatesequence]=max(GammaGMM,[],2);
%     [kmeansidx,~]=kmeans(data,L);
%     kmeansidx=truestates;
    for i=1:L
        datai=data(1+(i-1)*datailen:i*datailen,:);
%         meani=mean(datai);
%         covi=cov(datai);
%         struct.theta{i,1}.mu=meani;
%         datai=data(GMMstatesequence==i,:);
%         datai=data(kmeansidx==i,:);
        [leni,~]=size(datai);
        if leni >1
            struct.theta{i,1}.alpha_a=diminfo.inithetaparameters.alpha_a;
            struct.theta{i,1}.alpha_b=diminfo.inithetaparameters.alpha_b;
            for j=1:v
%                 struct.theta{i,1}.alpha(j,1)=gamrnd(struct.theta{i,1}.alpha_a,1/struct.theta{i,1}.alpha_b(j,1));
                struct.theta{i,1}.alpha(j,1)=struct.theta{i,1}.alpha_a/struct.theta{i,1}.alpha_b(j,1);
            end
%             mean_alpha=struct.theta{i,1}.alpha_a./struct.theta{i,1}.alpha_b;
            [structdatai,~] = ppca_ML(datai,v);
            struct.theta{i,1}.W_mean=structdatai.W;
            struct.theta{i,1}.W_sigma=diag(1./struct.theta{i,1}.alpha);
            struct.theta{i,1}.W=struct.theta{i,1}.W_mean;
            struct.theta{i,1}.tau_c=diminfo.inithetaparameters.tau_c;
            struct.theta{i,1}.tau_d=diminfo.inithetaparameters.tau_d;
%             struct.theta{i,1}.tau=gamrnd(struct.theta{i,1}.tau_c,1/struct.theta{i,1}.tau_d);
            struct.theta{i,1}.tau=1;
            struct.theta{i,1}.x_mean=mean(structdatai.x.mean);
            struct.theta{i,1}.x_sigma=structdatai.x.Sigma;
            struct.theta{i,1}.mu_mean=mean(datai)';
            struct.theta{i,1}.mu_sigma=diminfo.inithetaparameters.beta*diag(ones(q,1));
            struct.theta{i,1}.mu=mvnrnd(struct.theta{i,1}.mu_mean,struct.theta{i,1}.mu_sigma)';
            struct.theta{i,1}.real_v=v;
        else
            struct.theta{i,1}.alpha_a=diminfo.inithetaparameters.alpha_a;
            struct.theta{i,1}.alpha_b=diminfo.inithetaparameters.alpha_b;
            for j=1:v
%                 struct.theta{i,1}.alpha(j,1)=gamrnd(struct.theta{i,1}.alpha_a,1/struct.theta{i,1}.alpha_b(j,1));
                struct.theta{i,1}.alpha(j,1)=struct.theta{i,1}.alpha_a/struct.theta{i,1}.alpha_b(j,1);

            end
            [structdata,~] = ppca_ML(data,v);
            struct.theta{i,1}.W_mean=structdata.W;
            struct.theta{i,1}.W_sigma=diag(1./struct.theta{i,1}.alpha);
            for j=1:v
                struct.theta{i,1}.W(:,j)=mvnrnd(struct.theta{i,1}.W_mean(:,j),diag(1/struct.theta{i,1}.alpha(j,1)*ones(q,1)));
            end
            struct.theta{i,1}.tau_c=diminfo.inithetaparameters.tau_c;
            struct.theta{i,1}.tau_d=diminfo.inithetaparameters.tau_d;
%             struct.theta{i,1}.tau=gamrnd(struct.theta{i,1}.tau_c,1/struct.theta{i,1}.tau_d);
            struct.theta{i,1}.tau=1;
            struct.theta{i,1}.x_mean=mean(structdata.x.mean);
            struct.theta{i,1}.x_sigma=structdata.x.Sigma;
            struct.theta{i,1}.mu_mean=mean(data)';
            struct.theta{i,1}.mu_sigma=diminfo.inithetaparameters.beta*diag(ones(q,1));
            struct.theta{i,1}.mu=mvnrnd(struct.theta{i,1}.mu_mean,struct.theta{i,1}.mu_sigma)';
            struct.theta{i,1}.real_v=v;
        end

    end
    struct.hyperparameters=hyperparameters;
    alpha=gamrnd(hyperparameters.alpha_Gamma_a,1/(hyperparameters.alpha_Gamma_b));
    gamma=gamrnd(hyperparameters.gamma_Gamma_a,1/(hyperparameters.gamma_Gamma_b));
    beta=dirichlet(repmat(gamma/L,1,L),1);
    TPM=zeros(L,L);
    TPMraw=zeros(L,L);
    for i=1:L
        c=alpha*beta;
        TPMraw(i,:)=dirichlet(c,1);
        TPM(i,:)=TPMraw(i,:)/(1-TPMraw(i,i));
        TPM(i,i)=0;
    end
    struct.instantiation.gamma=gamma;
    struct.instantiation.alpha=alpha;
    struct.instantiation.beta=beta;
    struct.instantiation.TPM=TPM;
    struct.instantiation.TPMraw=TPMraw;
    struct.Dirichlet.beta=repmat(gamma/L,1,L);
    for i=1:L
        c=alpha*beta;
        struct.Dirichlet.TPMraw(i,:)=c;
    end
    struct.Dirichlet.DUM=repmat(diminfo.du0,L,diminfo.D);
    DUM=zeros(L,diminfo.D);
    for i=1:L
        c=struct.Dirichlet.DUM(i,:);
        DUM(i,:)=dirichlet(c,1);
    end
    struct.instantiation.DUM=DUM;
end

