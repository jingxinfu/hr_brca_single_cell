scPipe --outdir data/workflow/CCG1112_ATAC_QC/Tumor_MO/ scenic-egrns-mo \
    --cistopic_object data/workflow/CCG1112_ATAC_QC/Tumor_MO/cistopic_w_model.pickle \
    --h5ad  data/workflow/GEX_CCG1112_LowMt/gex_qc.h5ad \
    --motifs_menr  data/workflow/CCG1112_ATAC_QC/Tumor_MO/motifs/menr.pkl \
    --variable Cellstate \
    --metadata data/result/cleaned_files/Cellstate_MO.csv \
    --tf_file data/external/TF_db/TF_names_v_1.01.txt \
    --n_cpu 90