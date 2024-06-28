process MAKE_INPUT_CSV {

    label 'prepare_sample_sheet'

    input:
    path(in_files)
    val(pairedness)
    
    output:
    path "samplesheet.csv"

    """
    // set SE or PE value
    def pairednessStr = "paired_end"
    def pairednessFlag = "False"

    if ($pairedness) {
        pairednessStr = "single_end"
        pairednessFlag = "True"
    }

    // write header
    // def outfile = new File("samplesheet.csv")
    // outfile.append("sample,fastq_1,fastq_2,$pairednessStr\n")

    // get sample names
    // def getR2SampleNames = { it.getSimpleName() }
    // def r2Files = listOfR2Files.collect(getR2SampleNames)

    // write content
    // for (r2File : r2Files) {
    //     def r1File = r2File.replaceAll(/_R2/, "_R1")
    //     def sampleName = r2File.split('_')[0]
    //     outfile.append("$sampleName,$r1File,$r2File,$pairednessFlag\n")
    // }

    println("CSV file created at: $outfile")

    return $in_files
    """

}

workflow {
    fasta_files = channel.fromFilePairs(params.folder + '*_R{1,2}*')
    pairedness_flag = channel.of('False')
    MAKE_INPUT_CSV(fasta_files, pairedness_flag) | view { it }
}
