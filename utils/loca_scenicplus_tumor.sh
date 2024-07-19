scPipe --outdir /tumor_mo/ scenic-egrns-mo \
    --cistopic_object /mnt/disks/brca_16_466_data/workflow/CCG1112_ATAC_QC/Tumor_MO/cistopic_w_model.pickle \
    --h5ad  /mnt/disks/brca_16_466_data/workflow/GEX_CCG1112_LowMt/gex_qc.h5ad \
    --motifs_menr  /mnt/disks/brca_16_466_data/workflow/CCG1112_ATAC_QC/Tumor_MO/motifs/menr.pkl \
    --variable Cellstate \
    --metadata /mnt/disks/brca_16_466_data/result/cleaned_files/Cellstate_MO.csv \
    --tf_file /mnt/disks/brca_16_466_data/external/TF_db/TF_names_v_1.01.txt \
    --n_cpu 80 
