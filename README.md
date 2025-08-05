# &beta;

## To-do for version 2

### Features to implement

1. Include plotting motifs around TSSs on both strands separately to check if a pyrimidine-purine (initiator-like) motif is present on both strands. This lets a user check if TSSs are shifted (are not a pyrimidine-purine pair) and/or initiator motifs are different on the two strands (neither should happen).

2. Track generation for the genome browser (normalized counts).

3. Investigate and ideally resolve the issue with `CAGEr` using only one thread when reading samples and working within the pipeline. Get in touch with Charles Plessy after a reasonable investigation. (Damir discovered that CAGEr uses the number of thread equal to the number of read input files, independently of the number of threads set to it; but it is still unclear why CAGEr uses only one thread for multiple input samples when run within the pipeline.)

### Finishing up

4. Cite in `CITATIONS.md` all the tools that we used.

5. Make it possible to run the pipeline by providing the GitHub repository name (and, possibly, a version name / commit hash), instead of making the user clone the repository first.

- [ReadMe](docs/README.md)
  - The actual documentation
