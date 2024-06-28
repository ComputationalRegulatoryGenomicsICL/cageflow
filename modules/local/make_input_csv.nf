process MAKE_INPUT_CSV {

    label 'prepare_sample_sheet'

    input:
    path in_folder
    
    output:
    path "samplesheet.csv" outfile

    """
    // identify if SE or PE
    def pairedness = "paired_end"
    def pairednessFlag = "False"

    def listOfR2Files = file("$in_folder/**_R2*")

    if (listOfR2Files) {
        pairedness = "single_end"
        pairednessFlag = "True"
    }

    // write header
    def outfile = new File(options.o)
    outfile.append("sample,fastq_1,fastq_2,$pairedness\n")

    // get sample names
    def getR2SampleNames = { it.getSimpleName() }
    def r2Files = listOfR2Files.collect(getR2SampleNames)

    // write content
    for (r2File : r2Files) {
        def r1File = r2File.replaceAll(/_R2/, "_R1")
        def sampleName = r2File.split('_')[0]
        // check if R1 truly exists
        if (r1Files*.exists()) {
            throw new IOException("Missing read 1 for $sampleName")
        }
        outfile.append("$sampleName,$r1File,$r2File,$pairednessFlag\n")
    }

    println("CSV file created at: $outfile")
    """

}

workflow {
    MAKE_INPUT_CSV | view { it }
}
