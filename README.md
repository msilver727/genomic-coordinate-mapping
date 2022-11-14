# genomic-coordinate-mapping
This project enables mapping between transcript and genomic coordinates.

Two input files are required for the translation.
Input file 1 (transcripts input) contains a set of transcripts, each one consisting of a name, chromosome, start position, and CIGAR string. Can optionally include a fifth column for orientation, which must be provided as either "forward" or "reverse". If the fifth column is not provided, all transcripts will be assumed as "forward" orientation (5' to 3').
Input file 2 (queries input) contains a set of queries of specific coordinates to be mapped to the opposite coordinate system.

The function produces an output file of each query provided in input file 2 along with the corresponding position in the opposite coordinate system.

The default behavior is to convert from transcript position to genomic position, but the opposite translation can be enabled by setting `start_from_genomic` as True.
