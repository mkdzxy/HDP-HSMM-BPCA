clear
clc




sigmacell=cell(4,1);
sigmacell{1,1}=5*ones(1,3);
sigmacell{2,1}=4*ones(1,6);
sigmacell{3,1}=3*ones(1,9);
sigmacell{4,1}=2*ones(1,12);

meanvec=[-2,-1,1,2];

re = HDP_HSMM_BPCA_data_simulate(randmatrix_geo(4,0.45,0.55,300),meanvec,sigmacell,q=50,T=500,n=4);


[struct,diminfo,state_message,stateinfo,auxiliary_variable] = main_HDP_HSMM_BPCA_Bloked_Gibbs(re.data,10,500,50,100,400);







re = HDP_HSMM_BPCA_data_simulate(randmatrix_geo(2,0.45,0.55,300),meanvec,sigmacell,q=100,T=500,n=2);

[struct,diminfo,state_message,stateinfo,auxiliary_variable] = main_HDP_HSMM_BPCA_Bloked_Gibbs(re.data,10,150,50,100,300);



