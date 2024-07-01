function [struct,diminfo,state_message,stateinfo,auxiliary_variable] = main_HDP_HSMM_BPCA_Bloked_Gibbs(data,L,it,thetatimes,hypertimes,seed)
    for i=1:it
        if i==1
            [struct,diminfo] = ini_HDP_HSMM_BPCA(data,L,seed);
            T=diminfo.T;
            fprintf("Initialized, it=%d\n",i)
        else
            state_message = backwardmessage_statesequence_HDP_HSMM_BPCA(struct,diminfo);
            state_sequence=state_message.state_sequence;
            if i==2
                state_change=1;
                old_state_sequence=state_sequence;
                max_changestate=0;
            else
                state_change=sum(state_sequence~=old_state_sequence)/T;
                change_sequence=old_state_sequence(state_sequence~=old_state_sequence);
                [lenchange,~]=size(change_sequence);
                if lenchange>=1
                    change_sequence_tab=tabulate(change_sequence);
                    [~,maxindex]=max(change_sequence_tab(:,3));
                    max_changestate=change_sequence_tab(maxindex,1);
                else
                    max_changestate=0;
                end
                old_state_sequence=state_sequence;
            end
            auxiliary_variable = HDP_HSMM_BPCA_sample_auxiliary_variable(struct,diminfo,state_sequence);
            [struct,diminfo] = update_pi_beta_theta_hyper_HDP_HSMM_BPCA(struct,diminfo,state_sequence,auxiliary_variable,thetatimes,hypertimes);
            tabulate_stateseq=tabulate(state_sequence);
            fprintf("Sampling, it=%d, state_change=%d, max_changestate=%d, tabulate is following:\n",i,state_change,max_changestate)
            disp(tabulate_stateseq)
        end
    end
    T=diminfo.T;
    L=diminfo.L;
    state_frquency=zeros(L,1);
    for k=1:L
        state_frquency(k,1)=sum(state_sequence==k)/T;
    end
    stateinfo.state_sequence=state_sequence;
    stateinfo.state_frequencey=state_frquency;
    stateinfo.state_change=state_change;

end

