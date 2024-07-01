function re = simulation2datagenerate(P,seed,varargin)
    vars=["T","q","n"];
    values=[500,100,2];
    p=inputParser;
    for i=1:length(vars)
        addParameter(p,vars(i),values(i));
    end
    parse(p,varargin{:});
    T=p.Results.T;
    q=p.Results.q;
    n=p.Results.n;
    rng(seed)
    meanvec=[-1,1];
    for i=1:2
%         eqmean=length(meancell{i,1});
        meancell{i,1}=repmat(meanvec(1,i),1,q);
%         meancell{i,1}=[repmat(meanvec(1,i),1,eqsigma),zeros(1,q-eqsigma)];
    end

    M1=diag(rand(50,1));
    Z1=orth(rand(50,50));
    A1=Z1'*M1*Z1;
    A1(logical(eye(size(A1))))=diag(A1)+4;
%     diagA1=diag(A1);
%     A1(logical(eye(size(A1))))=diag(A1)+4;
%     A1=A1*10;
%     A1(logical(eye(size(A1))))=diagA1;




    M2=diag(rand(50,1));
    Z2=orth(rand(50,50));
    A2=Z2'*M2*Z2;
    A2(logical(eye(size(A2))))=diag(A2)+6;
%     A2(logical(eye(size(A2))))=diag(A2)+6;
%     A2(logical(eye(size(A2))))=rand(50,1);  



    sigmacell{1,1}=[A1,zeros(50);zeros(50),0.1*eye(50)];
    sigmacell{2,1}=[0.1*eye(50),zeros(50);zeros(50),A2];




%     sigmacell{1,1}=[0.01*diag(ones(30,1)),zeros(30,60);
%         zeros(30,30),3*ones(30,30),zeros(30,30);
%         zeros(30,60),3*ones(30,30)];
%     diag1=diag(sigmacell{1,1});
%     sigmacell{1,1}=0.3*sigmacell{1,1};
%     sigmacell{1,1}(logical(eye(size(sigmacell{1,1}))))=diag1;
%     
%     sigmacell{2,1}=[4*ones(30,30),zeros(30,60);
%         zeros(30,30),0.01*diag(ones(30,1)),zeros(30,30);
%         zeros(30,60),4*ones(30,30)];
%     diag2=diag(sigmacell{2,1});
%     sigmacell{2,1}=0.3*sigmacell{2,1};
%     sigmacell{2,1}(logical(eye(size(sigmacell{2,1}))))=diag2;
% 
%     sigmacell{3,1}=[5*diag(ones(30,1)),zeros(30,60);
%         zeros(30,30),5*ones(30,30),zeros(30,30);
%         zeros(30,30),zeros(30,30),0.01*diag(ones(30,1));];
%     diag3=diag(sigmacell{3,1});
%     sigmacell{3,1}=0.3*sigmacell{3,1};
%     sigmacell{3,1}(logical(eye(size(sigmacell{3,1}))))=diag3;

    
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

