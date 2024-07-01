function state_message = backwardmessage_statesequence_HDP_HSMM_BPCA(struct,diminfo)

    TPM=struct.instantiation.TPM;
    DUM=struct.instantiation.DUM;
    T=diminfo.T;
    L=diminfo.L;
    q=diminfo.q;
    D=diminfo.D;
%     v=diminfo.v;
    data=diminfo.data;
    B=zeros(L,T);%t=1:T
    Bast=zeros(L,T);%t=0:T-1
    PPCApdf=zeros(L,T);
    for t=1:T
        for k=1:L
            real_v=struct.theta{k,1}.real_v;
            x_sigma=inv(eye(real_v)+struct.theta{k,1}.tau*struct.theta{k,1}.W'*struct.theta{k,1}.W);
            x_sigma=(x_sigma+x_sigma')/2;
            x_mean=struct.theta{k,1}.tau*x_sigma*struct.theta{k,1}.W'*(data(t,:)'-struct.theta{k,1}.mu);
%             x_sigma=struct.theta{k,1}.x_sigma;
%             [lenx,~]=size(struct.theta{k,1}.x_mean);
%             if lenx==1
%                 x_mean=struct.theta{k,1}.x_mean;
%             else
%                 x_mean=mean(struct.theta{k,1}.x_mean);
%             end
%             x_sample=mvnrnd(x_mean,x_sigma)';
%             x_sample=x_mean;
            PPCApdf(k,t)=numerical_overflow(mvnpdf(data(t,:),(struct.theta{k,1}.W*x_mean+struct.theta{k,1}.mu)',1/struct.theta{k,1}.tau*eye(q)+struct.theta{k,1}.W*x_sigma*struct.theta{k,1}.W'));
        end
        PPCApdf(:,t)=PPCApdf(:,t)./sum(PPCApdf(:,t));
    end

    for t=T:-1:0
        if t==T
            B(:,T)=repmat(1/L,L,1);
        elseif t<T&& t>0
            Dast=min(T-t,D);
            sum1=zeros(L,1);
            for j=1:L
                for d=1:Dast
                    sum1(j,1)=sum1(j,1)+B(j,t+d)*DUM(j,d)*prod(PPCApdf(j,t+1:t+d));
                end
            end
            if T-t<D
                sum2=zeros(L,1);
                for j=1:L
                    sum2(j,1)=(1-sum(DUM(j,1:T-t)))*prod(PPCApdf(j,t+1:T));
                    Bast(j,t+1)=numerical_overflow(sum1(j,1)+sum2(j,1));
                end
                Bast(:,t+1)=Bast(:,t+1)./sum(Bast(:,t+1));
                for j=1:L
                    for k=1:L
                        B(j,t)=B(j,t)+numerical_overflow(Bast(k,t+1)*TPM(j,k));
                    end
                end
                B(:,t)=B(:,t)./sum(B(:,t));
            else
                for j=1:L
                    Bast(j,t+1)=numerical_overflow(sum1(j,1));
                end
                Bast(:,t+1)=Bast(:,t+1)./sum(Bast(:,t+1));
                for j=1:L
                    for k=1:L
                        B(j,t)=B(j,t)+numerical_overflow(Bast(k,t+1)*TPM(j,k));
                    end
                end
                B(:,t)=B(:,t)./sum(B(:,t));
            end
            
            
        elseif t==0
            sum1=zeros(L,1);
            for j=1:L
                for d=1:D
                    sum1(j,1)=sum1(j,1)+B(j,t+d)*DUM(j,d)*prod(PPCApdf(j,t+1:t+d));
                end
                Bast(j,t+1)=numerical_overflow(sum1(j,1));
            end
            Bast(:,t+1)=Bast(:,t+1)./sum(Bast(:,t+1));

        end
    end

    statelen=0;
    sampletimes=1;
    state_sequence=zeros(T,1);
    while statelen< T
        if sampletimes==1
            pstate=Bast(:,1);
            r1=rand();
            a=find((cumsum(pstate)-r1)>=0);
            currentstate=a(1,1);
            for d=1:D
                pdu(d,1)=numerical_overflow(DUM(currentstate,d)*prod(PPCApdf(currentstate,1:d))*B(currentstate,d)/Bast(currentstate,1));
            end
            pdu=pdu/sum(pdu);
            r2=rand();
            b=find((cumsum(pdu)-r2)>=0);
            currentdu=b(1,1);
            state_sequence(1:currentdu,1)=currentstate;
            statelen=statelen+currentdu;
            sampletimes=sampletimes+1;
        else
            if T-statelen>=D
                for j=1:L
                    pstate(j,1)=numerical_overflow(Bast(j,1+statelen)*TPM(currentstate,j));
                end
                pstate=pstate/sum(pstate);
                r1=rand();
                a=find((cumsum(pstate)-r1)>=0);
                currentstate=a(1,1);
                for d=1:D
                    pdu(d,1)=numerical_overflow(DUM(currentstate,d)*prod(PPCApdf(currentstate,1+statelen:d+statelen))*B(currentstate,d+statelen)/Bast(currentstate,1+statelen));
                end
                pdu=pdu/sum(pdu);
                r2=rand();
                b=find((cumsum(pdu)-r2)>=0);
                currentdu=b(1,1);
                state_sequence(statelen+1:statelen+currentdu,1)=currentstate;
                statelen=statelen+currentdu;
                sampletimes=sampletimes+1;
            elseif T-statelen<D
                Dast=T-statelen;
                for j=1:L
                    pstate(j,1)=numerical_overflow(Bast(j,1+statelen)*TPM(currentstate,j));
                end
                pstate=pstate/sum(pstate);
                r1=rand();
                a=find((cumsum(pstate)-r1)>=0);
                currentstate=a(1,1);
                for d=1:Dast
                    pdu(d,1)=numerical_overflow(DUM(currentstate,d)*prod(PPCApdf(currentstate,1+statelen:d+statelen))*B(currentstate,d+statelen)/Bast(currentstate,1+statelen));
                end
                pdu=pdu/sum(pdu);
                r2=rand();
                b=find((cumsum(pdu)-r2)>=0);
                currentdu=b(1,1);
                state_sequence(statelen+1:statelen+currentdu,1)=currentstate;
                statelen=statelen+currentdu;
                sampletimes=sampletimes+1;
            end
        end
    end
    state_sequence=state_sequence(1:T,1);
    state_message.B=B;
    state_message.Bast=Bast;
    state_message.state_sequence=state_sequence;
    state_message.PPCApdf=PPCApdf;


















%     for t=T:-1:1
%         if t==T
%             backwardmessage(:,t)=repmat(1/L,L,1);
%         else
%             for k=1:L
%                 sum1=0;
%                 for j=1:L
%                     sum1=sum1+numerical_overflow(backwardmessage(j,t+1)*PPCApdf(j,t)*TPM(k,j));
%                 end
%                 backwardmessage(k,t)=sum1;
%             end
%             backwardmessage(:,t)=backwardmessage(:,t)./sum(backwardmessage(:,t));
%         end
%     end
% 
%     for t=1:T
%         p=zeros(L,1);
%         if t==1
%             for k=1:L
%                 p(k,1)=numerical_overflow(PPCApdf(k,t)*backwardmessage(k,t));
%             end
%             p=p./sum(p);
%             r=rand();
%             a=find((cumsum(p)-r)>=0);
%             state_sequence(t,1)=a(1,1);
%         else
%             for k=1:L
%                 p(k,1)=numerical_overflow(TPM(state_sequence(t-1,1),k)*PPCApdf(k,t)*backwardmessage(k,t));
%             end
%             p=p./sum(p);
%             r=rand();
%             a=find((cumsum(p)-r)>=0);
%             state_sequence(t,1)=a(1,1);
%         end
%     end
% 
%     state_message.backwardmessage=backwardmessage;
%     state_message.state_sequence=state_sequence;


end

