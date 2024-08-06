// 
// CAGEr analysis steps
// 

workflow CAGER {

    take:
        cager_rds
        ch_versions
    
    main:

    emit:
        ch_versions

}