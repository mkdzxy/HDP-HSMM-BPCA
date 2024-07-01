function auxiliary_variables = HDP_HSMM_BPCA_sample_auxiliary_variable(struct,diminfo,state_sequence)
    L=diminfo.L;
    T=diminfo.T;
%     data=diminfo.data;
    seed=diminfo.seed;
    rng(seed)
    TPMraw=struct.instantiation.TPMraw;
    TFM=zeros(L,L); % transition frequency matrix
    for t=1:T-1
        TFM(state_sequence(t,1),state_sequence(t+1,1))=TFM(state_sequence(t,1),state_sequence(t+1,1))+1;
    end
    for j=1:L
        TFM(j,j)=0;
    end
    % sample \rho
    TFM_jdot=sum(TFM,2);
    rho=zeros(L,1);
    for j=1:L
        rho(j,1)=sum(geornd(1-TPMraw(j,j),1,TFM_jdot(j)));
    end
    % sample m

    m=zeros(L,L);
    for j=1:L
        for k=1:L
            if j==k
                sum1=0;
                for i=1:rho(j,1)
                    r=rand();
                    x=(r<=(struct.instantiation.alpha*struct.instantiation.beta(1,k))/(struct.instantiation.alpha*struct.instantiation.beta(1,k)+i));
                    sum1=sum1+x;
                end
                m(j,k)=sum1;
            else
                sum1=0;
                for i=1:TFM(j,k)
                    r=rand();
                    x=(r<=(struct.instantiation.alpha*struct.instantiation.beta(1,k))/(struct.instantiation.alpha*struct.instantiation.beta(1,k)+i));
                    sum1=sum1+x;
                end
                m(j,k)=sum1;
            end
        end
    end
   
    auxiliary_variables.m=m;
    auxiliary_variables.rho=rho;
    auxiliary_variables.TFM=TFM;
end

