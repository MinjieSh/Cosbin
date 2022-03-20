norm_data_list <- list(totalcount_out,
                       cosbin_out$data,
                       tmm_out,
                       DESeq2_out,
                       tcc_out)


##################################################################

# CEG dislocation: 
CEG_dislocation <- cosine_dislocation(data_grnd, CEG_grnd, norm_data_list)

# DEG dislocation (exclude sDEG): 
DEG_dislocation <- cosine_dislocation(data_grnd, DEG_grnd, norm_data_list)

# sDEG dislocation: 
sDEG_dislocation <- groupwise_cosine_dislocation(data_grnd, norm_data_list, nGroup, pSMG)

# ROC (iCEG vs non-iCEG): 
ROC_iCEG <- ROC(norm_data_list, CEG_grnd, "iCEG", c("#AA000066", "#00AA0066", "#0000AA66",  "#AAAA0066", "#00AAAA66"))

# ROC (DEG (exclude SMG) vs non DEG (exclude SMG)): 
ROC_DEG <- ROC(norm_data_list, DEG_grnd, "DEG", c("#AA000066", "#00AA0066", "#0000AA66",  "#AAAA0066", "#00AAAA66"))

# average_of_pairwise_Log_Fold_Change_Mean_Squared_Error: 
mean_LFC_MSE <- average_of_pairwise_LFC_MSE(norm_data_list, CEG_grnd, nGroup)
