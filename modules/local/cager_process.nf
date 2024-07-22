// 
// CAGEr analysis steps
// 

process CAGER_PROCESS {

    input:
    path(cagexp_object)

    output:
    path "versions.yml", emit: versions

    """
    """

}