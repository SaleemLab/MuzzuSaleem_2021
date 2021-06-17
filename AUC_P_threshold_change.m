
% Functions to look at visual responses
if ~exist('ProjectData','var')
    [ProjectData AM_UnitResponses AM_Param AM_Speed AM_UOI SelectedResponses AM_UnitResponses_smooth] = LoadDataALL;
end

% first 7 animals
if  size(ProjectData,1)>10
    CTRL_exp = 0;
    Animal_1st_idx = [1 5 7 12 15 24 31];
    
    if ~exist('PertResp_units','var')
        % select only perturbation responsive units
        load('AUC_shuffled.mat')
        Sh_responses = AUC_shuffled(:,2:end);
        p_pert_th = prctile(Sh_responses(:),99);
        PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
        % select only pos. modulated perturbation responsive units
        load('DM_pert_shuffled.mat')
        DM = DM_sh(:,1);
        DM_sign_i(:,1) = DM>0;
        DM_sign_i(:,2) = DM<=0;
        % select only pos. modulated perturbation responsive units
        PertRespUnits_pos = PertResp_units & DM_sign_i(:,1);
        PertRespUnits_neg = PertResp_units & DM_sign_i(:,2);
    end
    
else
    CTRL_exp = 1;
    % naive animals
    Animal_1st_idx = [1 4 7];
    % select only perturbation responsive units
    load('AUC_shuffled_CTRL_1.mat')
    Sh_responses = AUC_shuffled(:,2:end);
    p_pert_th = prctile(Sh_responses(:),99);
    PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
    % select only pos. modulated perturbation responsive units
    load('DM_pert_shuffled_CTRL.mat')
    DM = DM_sh(:,1);
    DM_sign_i(:,1) = DM>0;
    DM_sign_i(:,2) = DM<=0;
    % select only pos. modulated perturbation responsive units
    PertRespUnits_pos = PertResp_units & DM_sign_i(:,1);
    PertRespUnits_neg = PertResp_units & DM_sign_i(:,2);
end

%% plot % of units resulting as pert. resp as function of p of AUC
% manual 
%PertRespUnits_pos = PertResp_units & DM_sign_i(:,1);
PertRespUnits_pos_m = DM_sign_i(:,1);
%PertRespUnits_pos_m(146) = false; % for naive animals
%PertRespUnits_pos_m(185) = false; % for naive animals
k = 1; clear PertResp_units_var PertResp_units_perc
for i = 90:0.1:99.9
    p_pert_th_var = prctile(Sh_responses(:),i);
    PertResp_units_var(:,k) = (AUC_shuffled(:,1)>p_pert_th_var);
    PertResp_units_perc(k) = sum(PertResp_units_var(:,k) & PertRespUnits_pos_m)/length(PertRespUnits_pos_m)*100;
    k = k+1;  
end

 
figure
plot((100-[90:0.1:99.9])/100,PertResp_units_perc,'-or')
hold on
plot((100-[90:0.1:99.9])/100,PertResp_units_perc_CTRL,'-ok')
xlabel('threshold for AUC performace (p value)')
ylabel('% pert. resp. units (MI>0)')
set(gca,'box','off','TickDir','out', 'XTick',[0.001 0.01 0.05 0.1],'XScale','log')
grid on
xlim([-0.01 0.11])
legend('Experienced','naive')

PertResp_units_var_CTRL = PertResp_units_var;
PertResp_units_perc_CTRL = PertResp_units_perc;
save('PercentPertUnits_CTRL.mat','PertResp_units_var_CTRL','PertResp_units_perc_CTRL')
