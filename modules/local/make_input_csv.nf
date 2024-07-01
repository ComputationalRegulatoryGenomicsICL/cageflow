process MAKE_INPUT_CSV {

    label 'prepare_sample_sheet'

    input:
    val(inFile)
    
    output:
    stdout
    // path "results/samplesheet.csv"

    exec:
    outpath = "results/samplesheet.csv"
    outfile = file(outpath)
    if (! outfile.exists()){
            header = "sample,fastq_1,fastq_2,single_end\n"
            outfile.text = header
        }
    if (params.paired) {
        println("hello")
        inFile1 = inFile[1][0]
        inFile2 = inFile[1][1]
        inFile1File = file(inFile1)
        baseName = inFile1File.getSimpleName()
        sampleName = baseName.minus(~/_R\d.*$/)
        outfile.append(sampleName + ',' + inFile1 + ',' + inFile2 +',False\n')
    }else {
        println("hi")
        inFileFile = file(inFile)
        baseName = inFileFile.getSimpleName()
        sampleName = baseName.minus(~/_R1.*$/)
        outfile.append(sampleName + ',' + inFile + ',,True\n')
    }
    
    return outpath

}

workflow {
    if (params.paired) {
        myPath = channel.fromFilePairs("$params.folder/**_R{1,2}*")
        MAKE_INPUT_CSV(myPath)
    }else{
        myPath = channel.fromPath("$params.folder/**_R1*")
        MAKE_INPUT_CSV(myPath)
    }
        
}
