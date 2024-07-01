function re = HDP_HSMM_BPCA_data_simulate(P,meanvec,sigmacell,varargin)
% T: total length
% q: variable length
% n: state number
% P: transition probability matrix TPM
    vars=["T","q","n"];
    values=[1000,10,4];
    p=inputParser;
    for i=1:length(vars)
        addParameter(p,vars(i),values(i));
    end
    parse(p,varargin{:});
    T=p.Results.T;
    q=p.Results.q;
    n=p.Results.n;
    meancell=cell(n,1);
    for i=1:n
%         eqmean=length(meancell{i,1});
        eqsigma=length(sigmacell{i,1});
        meancell{i,1}=repmat(meanvec(1,i),1,q);
%         meancell{i,1}=[repmat(meanvec(1,i),1,eqsigma),zeros(1,q-eqsigma)];
        sigmacell{i,1}=diag([sigmacell{i,1},0.1*ones(1,q-eqsigma)]);
    end
    state_seqs=zeros(T,1);
    for i=1:T
        if i==1
            r=randi(n);
            state_seqs(i,1)=r;
        else
            runi=rand();
            a=find((cumsum(P(state_seqs(i-1,1),:))-runi)>=0);
            nextstate=a(1,1);
            state_seqs(i,1)=nextstate;
        end
    end
    data=zeros(T,q);
    for k=1:T
        data(k,:)=mvnrnd(meancell{state_seqs(k,1),1},sigmacell{state_seqs(k,1),1},1);
    end
    re.data=data;
    re.states=state_seqs;
    re.P=P;
    re.meancell=meancell;
    re.sigmacell=sigmacell;
    state_frequency=zeros(n,1);
    for i=1:n
        state_frequency(i,1)=sum(state_seqs==i)/T;
    end
    re.state_frequency=state_frequency;


end

