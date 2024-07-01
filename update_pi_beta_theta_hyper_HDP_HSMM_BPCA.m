function [struct_update,diminfo] = update_pi_beta_theta_hyper_HDP_HSMM_BPCA(struct,diminfo,state_sequence,auxiliary_variable,theta_times,hyperparameters_times)
    struct_update=struct;
    du0=diminfo.inihyperparameters.du0;
    D=diminfo.inihyperparameters.D;
    L=diminfo.L;
    q=diminfo.q;
%     v=diminfo.v;
    inihyperarameters=diminfo.inihyperparameters;
    inithetaparameters=diminfo.inithetaparameters;
    theta_alpha_threshold=inithetaparameters.theta_alpha_threshold;
    data=diminfo.data;
    seed=diminfo.seed;
    rng(seed)
    % update beta
    mdot=sum(auxiliary_variable.m,1);
    struct_update.Dirichlet.beta=struct.instantiation.gamma/L+mdot;
    struct_update.instantiation.beta=dirichlet(struct_update.Dirichlet.beta,1);

    % update TPM
    TFM=auxiliary_variable.TFM;
    rho=auxiliary_variable.rho;
    for k=1:L
        struct_update.Dirichlet.TPMraw(k,:)=struct.instantiation.alpha*struct_update.instantiation.beta+TFM(k,:);
        struct_update.Dirichlet.TPMraw(k,k)=struct.instantiation.alpha*struct_update.instantiation.beta(k)+rho(k,1);
        struct_update.instantiation.TPMraw(k,:)=dirichlet(struct_update.Dirichlet.TPMraw(k,:),1);
        struct_update.instantiation.TPM(k,:)=struct_update.instantiation.TPMraw(k,:)/(1-struct_update.instantiation.TPMraw(k,k));
        struct_update.instantiation.TPM(k,k)=0;
    end

    % update state duration distributions
    for j=1:L
        dutable = duration_calculate(state_sequence,j,D);
        struct_update.Dirichlet.DUM(j,:)=repmat(du0,1,D)+dutable(:,2)';
        struct_update.instantiation.DUM(j,:)=dirichlet(struct_update.Dirichlet.DUM(j,:),1);
    end
 

    % update theta
    state_frequency=zeros(L,1);
    for k=1:L
        datak=data(state_sequence==k,:);
        [lenk,~]=size(datak);
        state_frequency(k,1)=lenk;
        if lenk>=1
            real_v=struct_update.theta{k,1}.real_v;
            for i=1:theta_times
                struct_update.theta{k,1}.x_sigma=inv(eye(real_v)+struct_update.theta{k,1}.tau*struct_update.theta{k,1}.W'*struct_update.theta{k,1}.W);
                struct_update.theta{k,1}.x_sigma=(struct_update.theta{k,1}.x_sigma+struct_update.theta{k,1}.x_sigma')/2;
%                 struct_update.theta{k,1}.x_sigma=diag(diag(struct_update.theta{k,1}.x_sigma));
                struct_update.theta{k,1}.x_mean=zeros(lenk,real_v);
                struct_update.theta{k,1}.x=zeros(lenk,real_v);
                for j=1:lenk
                    struct_update.theta{k,1}.x_mean(j,:)=(struct_update.theta{k,1}.tau*struct_update.theta{k,1}.x_sigma*struct_update.theta{k,1}.W'*(datak(j,:)-struct_update.theta{k,1}.mu')')';
                    struct_update.theta{k,1}.x(j,:)=mvnrnd(struct_update.theta{k,1}.x_mean(j,:),struct_update.theta{k,1}.x_sigma);
%                     struct_update.theta{k,1}.x(j,:)=struct_update.theta{k,1}.x_mean(j,:);
                end
                struct_update.theta{k,1}.mu_sigma=1/(inithetaparameters.beta+lenk*struct_update.theta{k,1}.tau)*eye(q);
                sum1=0;
                for j=1:lenk
                    sum1=sum1+datak(j,:)'-struct_update.theta{k,1}.W*struct_update.theta{k,1}.x(j,:)';
                end
                struct_update.theta{k,1}.mu_mean=struct_update.theta{k,1}.tau*struct_update.theta{k,1}.mu_sigma*sum1;
                struct_update.theta{k,1}.mu=mvnrnd(struct_update.theta{k,1}.mu_mean',struct_update.theta{k,1}.mu_sigma)';
%                 sum2=0;
%                 for j=1:lenk
%                     sum2=sum2+struct_update.theta{k,1}.x(j,:)'*struct_update.theta{k,1}.x(j,:);
%                 end
%                 struct_update.theta{k,1}.W_sigma=inv(diag(struct_update.theta{k,1}.alpha)+struct_update.theta{k,1}.tau*sum2);
                struct_update.theta{k,1}.W_sigma=inv(diag(struct_update.theta{k,1}.alpha)+struct_update.theta{k,1}.tau*struct_update.theta{k,1}.x'*struct_update.theta{k,1}.x);
                struct_update.theta{k,1}.W_sigma=(struct_update.theta{k,1}.W_sigma+struct_update.theta{k,1}.W_sigma')/2;
                struct_update.theta{k,1}.W_mean=struct_update.theta{k,1}.tau*(datak-kron(ones(lenk,1),struct_update.theta{k,1}.mu'))'*struct_update.theta{k,1}.x*struct_update.theta{k,1}.W_sigma;
%                 struct_update.theta{k,1}.W=reshape(mvnrnd(reshape(struct_update.theta{k,1}.W_mean,[q*v,1])',kron(struct_update.theta{k,1}.W_sigma,eye(q))),[q,v]);
%                 struct_update.theta{k,1}.W=struct_update.theta{k,1}.W_mean;
                for j=1:q
                    struct_update.theta{k,1}.W(j,:)=mvnrnd(struct_update.theta{k,1}.W_mean(j,:),struct_update.theta{k,1}.W_sigma);
                end
                struct_update.theta{k,1}.alpha_a=inithetaparameters.alpha_a+q/2;
                for j=1:real_v
                    struct_update.theta{k,1}.alpha_b(j,1)=inithetaparameters.alpha_b(j,1)+struct_update.theta{k,1}.W(:,j)'*struct_update.theta{k,1}.W(:,j)/2;
                    struct_update.theta{k,1}.alpha(j,1)=gamrnd(struct_update.theta{k,1}.alpha_a,1/(struct_update.theta{k,1}.alpha_b(j,1)));
%                     struct_update.theta{k,1}.alpha(j,1)=struct_update.theta{k,1}.alpha_a/struct_update.theta{k,1}.alpha_b(j,1);
                end
                struct_update.theta{k,1}.tau_c=inithetaparameters.tau_c+lenk*q/2;
                sum3=0;
                for j=1:lenk
                    sum3=sum3+datak(j,:)*datak(j,:)'+struct_update.theta{k,1}.mu'*struct_update.theta{k,1}.mu+trace(struct_update.theta{k,1}.W'*struct_update.theta{k,1}.W*struct_update.theta{k,1}.x(j,:)'*struct_update.theta{k,1}.x(j,:))+2*struct_update.theta{k,1}.mu'*struct_update.theta{k,1}.W*struct_update.theta{k,1}.x(j,:)'-2*datak(j,:)*struct_update.theta{k,1}.W*struct_update.theta{k,1}.x(j,:)'-2*datak(j,:)*struct_update.theta{k,1}.mu;
                end
                struct_update.theta{k,1}.tau_d=inithetaparameters.tau_d+0.5*sum3;
                struct_update.theta{k,1}.tau=gamrnd(struct_update.theta{k,1}.tau_c,1/(struct_update.theta{k,1}.tau_d));
%                 struct_update.theta{k,1}.tau=struct_update.theta{k,1}.tau_c/struct_update.theta{k,1}.tau_d;
%                 fprintf("k=%d,i=%d\n",k,i)
            end
            % prune W based on alpha
            occupied_col=struct_update.theta{k,1}.alpha_b>theta_alpha_threshold;
            struct_update.theta{k,1}.alpha_b=struct_update.theta{k,1}.alpha_b(occupied_col);
            struct_update.theta{k,1}.alpha=struct_update.theta{k,1}.alpha(occupied_col);
            struct_update.theta{k,1}.real_v=sum(occupied_col);
            struct_update.theta{k,1}.W_mean=struct_update.theta{k,1}.W_mean(:,occupied_col);
            struct_update.theta{k,1}.W_sigma=struct_update.theta{k,1}.W_sigma(occupied_col,occupied_col);
            struct_update.theta{k,1}.W=struct_update.theta{k,1}.W(:,occupied_col);
            struct_update.theta{k,1}.x_sigma=struct_update.theta{k,1}.x_sigma(occupied_col,occupied_col);
            struct_update.theta{k,1}.x_mean=struct_update.theta{k,1}.x_mean(:,occupied_col);
            struct_update.theta{k,1}.x=struct_update.theta{k,1}.x(:,occupied_col);
        elseif lenk==0
%                 struct_update.theta{k,1}.mu=mvnrnd(struct_update.theta{k,1}.mu_mean,struct_update.theta{k,1}.mu_sigma,1);
%                 struct_update.theta{k,1}.sigma=inv(wishrnd(inv(struct_update.theta{k,1}.IW_S),struct_update.theta{k,1}.IW_v));
              continue
        end

    end

    % update hyperparameters

    % sample \alpha
    r=zeros(L,1);
    s=zeros(L,1);
    m=auxiliary_variable.m;
    for i=1:hyperparameters_times
        if i==1
            alpha=struct_update.instantiation.alpha;
            for k=1:L
                r(k,1)=betarnd(alpha+1,sum(TFM(k,:))+rho(k,1),1);
                s(k,1)=rand()<((sum(TFM(k,:))+rho(k,1))/(sum(TFM(k,:))+alpha+rho(k,1)));
            end
            alpha=gamrnd(inihyperarameters.alpha_Gamma_a+sum(m,"all")-sum(s),1/(inihyperarameters.alpha_Gamma_b-sum(log(r))));
        else
            for k=1:L
                r(k,1)=betarnd(alpha+1,sum(TFM(k,:)),1);
                s(k,1)=rand()<(sum(TFM(k,:))/(sum(TFM(k,:))+alpha));
            end
            alpha=gamrnd(inihyperarameters.alpha_Gamma_a+sum(m,"all")-sum(s),1/(inihyperarameters.alpha_Gamma_b-sum(log(r))));
        end
    end
    struct_update.hyperparameters.alpha_Gamma_a=inihyperarameters.alpha_Gamma_a+sum(m,"all")-sum(s);
    struct_update.hyperparameters.alpha_Gamma_b=inihyperarameters.alpha_Gamma_b-sum(log(r));
    struct_update.instantiation.alpha=gamrnd(struct_update.hyperparameters.alpha_Gamma_a,1/(struct_update.hyperparameters.alpha_Gamma_b));

    % sample \gamma

    m=auxiliary_variable.m;
    for i=1:hyperparameters_times
        if i==1
            gamma=struct_update.instantiation.gamma;
            eta=betarnd(gamma+1,sum(m,"all"));
            zeta=rand()<(sum(m,"all")/(sum(m,"all")+gamma));
            gamma=gamrnd(inihyperarameters.gamma_Gamma_a+L-zeta,1/(inihyperarameters.gamma_Gamma_b-log(eta)));
        else
            eta=betarnd(gamma+1,sum(m,"all"));
            zeta=rand()<(sum(m,"all")/(sum(m,"all")+gamma));
            gamma=gamrnd(inihyperarameters.gamma_Gamma_a+L-zeta,1/(inihyperarameters.gamma_Gamma_b-log(eta)));
        end
    end
    struct_update.hyperparameters.gamma_Gamma_a=inihyperarameters.gamma_Gamma_a+L-zeta;
    struct_update.hyperparameters.gamma_Gamma_b=inihyperarameters.gamma_Gamma_b-log(eta);
    struct_update.instantiation.gamma=gamrnd(struct_update.hyperparameters.gamma_Gamma_a,1/(struct_update.hyperparameters.gamma_Gamma_b));

    % prune L
    occupied_state=state_frequency>0;
    real_L=sum(occupied_state>0);
    struct_update.theta=struct_update.theta(occupied_state,:);
    struct_update.instantiation.beta=struct_update.instantiation.beta(1,occupied_state);
    struct_update.instantiation.TPM=struct_update.instantiation.TPM(occupied_state,occupied_state);
    struct_update.instantiation.DUM=struct_update.instantiation.DUM(occupied_state,:);
    struct_update.instantiation.TPMraw=struct_update.instantiation.TPMraw(occupied_state,occupied_state);
    struct_update.Dirichlet.beta=struct_update.Dirichlet.beta(1,occupied_state);
    struct_update.Dirichlet.TPMraw=struct_update.Dirichlet.TPMraw(occupied_state,occupied_state);
    struct_update.Dirichlet.DUM=struct_update.Dirichlet.DUM(occupied_state,:);
    diminfo.L=real_L;
    

end

