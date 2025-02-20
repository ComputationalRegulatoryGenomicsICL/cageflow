// Creating markdown report from results

process CAGER_REPORT {
    label 'process_medium'
    stageInMode 'copy'

    input:
    path rmarkd_template
    tuple path(tss_hm_ta_plots), path(tss_hm_ta_data)
    path corr_m
    tuple path(cc_iqw_rc_plots), path(cc_txt), path(cc_iqw_rc_data)
    tuple path(tca_dn_n_plots), path(tca_dn_n_data)

    output:
    path "*.html"

    """
    #!/usr/bin/env Rscript
    library(rmarkdown)

    heatmap_cex_row <- ${heatmap_cex_row}
    rmarkdown::render('${rmarkd_template}')
    """
}