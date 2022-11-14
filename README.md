# genomic-coordinate-mapping
This project enables mapping between transcript and genomic coordinates.

Two input files are required for the translation.
Input file 1 contains a set of transcripts, each one consisting of a name, chromosome, start position, and CIGAR string.
Input file 2 contains a set of queries of specific coordinates to be mapped to the opposite coordinate system.

The function produces an output file of each query provided in input file 2 along with the corresponding position in the opposite coordinate system.

The default behavior is to convert from transcript position to genomic position, but the opposite translate can be enabled by setting `start_from_genomic` as True.

There is an optional argument for using reverse orientation instead of forward orientation for the transcript coordinates.
The limitation for this orientation is that the start position for each transcript only corresponds to the forward orientation,
so there is no way of knowing what the starting transcript position would be for the reverse orientation.
