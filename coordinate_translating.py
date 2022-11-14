import argparse
import re
import os


# Break the CIGAR string into each component part
# Compile single time for improved efficiency
cigar_pttn = re.compile(r'(\d+)([MXDIN])')


def get_converted_position(transcript_info, requested_pos, reverse_orientation=False, start_from_genomic=False):
    """ Convert from requested position in one coordinate system to corresponding position in the opposite coordinate system.
    By default, convert from transcript position to a corresponding genomic position, but pass
    start_from_genomic=True to convert from genomic position to a corresponding transcript position.
    Can pass reverse_orientation=True to use 3' to 5' instead of the default of 5' to 3'
    """

    # Each segment consisting of a segment length and cigar type
    parsed_cigar = cigar_pttn.findall(transcript_info['cigar'])

    # If reverse orientation is requested, just reverse the cigar segments
    # NOTE: does not account for potentially different transcript starting position
    if reverse_orientation:
        parsed_cigar = parsed_cigar[::-1]

    # Use the starting point of the transcript as the initial genomic position
    genomic_pos = transcript_info['start']

    # Transcript position will always start at 0
    transcript_pos = 0

    # Iterate through the cigar string to keep incrementing the genomic and transcript positions until the requested position is reached
    for segment_length, cigar_type in parsed_cigar:
        segment_length = int(segment_length)

        # Use the position from the starting coordinate system as the "current position of interest"
        current_pos_of_interest = genomic_pos if start_from_genomic else transcript_pos

        # Break from loop once the current position of interest is equal to or beyond the requested position in the starting coordinate system
        if current_pos_of_interest >= requested_pos:
            break

        # For either matches ('M') or mismatches ('X'), increase both the genomic position and transcript position by
        # the remaining number of bases left in the cigar segment until the requested position is reached for the starting coordinate system
        if cigar_type in ('M', 'X'):
            bases_into_cigar = min(segment_length, requested_pos - current_pos_of_interest)
            genomic_pos += bases_into_cigar
            transcript_pos += bases_into_cigar

        # For deletions ('D') or gaps ('N'), increase the genomic position by full segment length without altering the transcript position
        elif cigar_type in ('D', 'N'):
            genomic_pos += segment_length

        # For insertions ('I'), increase the transcript position by full segment length without altering the genomic position
        elif cigar_type == 'I':
            transcript_pos += segment_length

        # Uncomment following print statement if needed for testing purposes
        # print('cigar: {}{}, requested: {}, current: {}, transcript: {}, genomic: {}, bases_into_cigar: {}'.format(
        #     segment_length, cigar_type, requested_pos, current_pos_of_interest, transcript_pos, genomic_pos, bases_into_cigar))

    # Return the final position in the opposite coordinate system from the starting coordinate system
    final_pos_of_interest = transcript_pos if start_from_genomic else genomic_pos
    return final_pos_of_interest


def translate_coordinates(transcripts_input_p, queries_input_p, output_file_p, reverse_orientation, start_from_genomic):
    """ Given two input files, one containing information about various transcripts and another containing information
    regarding a requested coordinate within a transcript, output the set of translated coordinates."""

    # Validate that the input file paths point to actual files
    if not os.path.exists(transcripts_input_p):
        raise Exception('Invalid transcripts input path provided, {} does not point to a valid file.'.format(transcripts_input_p))
    if not os.path.exists(queries_input_p):
        raise Exception('Invalid queries input path provided, {} does not point to a valid file.'.format(queries_input_p))

    # Collect info from transcripts input into a dictionary of each transcript
    transcripts = {}
    with open(transcripts_input_p) as transcripts_input:
        for row in transcripts_input:
            cols = row.strip().split('\t')

            # Validate that the right number of columns were provided
            if len(cols) != 4:
                raise Exception('Invalid transcripts input provided. Expecting 4 columns, received {} columns.'.format(len(cols)))

            # Grab data from the columns
            transcript_name, chrom, start, cigar = cols

            # Validate that the start is an integer
            try:
                start = int(start)
            except ValueError:
                raise ValueError('Invalid start position of `{}` provided, must be an integer'.format(start))

            # Validate transcript uniqueness
            if transcript_name in transcripts:
                raise Exception('Duplicate transcripts provided, each transcript should be unique')

            transcripts[transcript_name] = {'chrom': chrom, 'start': start, 'cigar': cigar}

    # Read through queries input to translate each requested position into the opposite coordinate system
    with open(queries_input_p) as queries_input, open(output_file_p, 'w') as output_file:
        for row in queries_input:
            cols = row.strip().split('\t')

            # Validate that the right number of columns were provided
            if len(cols) != 2:
                raise Exception('Invalid queries input provided. Expecting 2 columns, received {} columns.'.format(len(cols)))

            # Grab data from the columns
            transcript_name, requested_pos = cols

            # Validate that the requested position is an integer
            try:
                requested_pos = int(requested_pos)
            except ValueError:
                raise ValueError('Invalid requested position of `{}` provided, must be an integer'.format(requested_pos))

            # Grab relevant transcript info and validate that the transcript name was provided in the transcripts input file
            transcript_info = transcripts.get(transcript_name)
            if not transcript_info:
                raise Exception('Unknown transcript provided in queries input: `{}`'.format(transcript_name))

            # Convert the requested position from its original coordinate system to the opposite coordinate system
            # (default is to go from transcript coordinate to genomic coordinate, but vice versa is also possible)
            converted_pos = get_converted_position(transcript_info, requested_pos,
                                                   reverse_orientation=reverse_orientation, start_from_genomic=start_from_genomic)

            # Combine all of the outputs together and write to the output file
            output = cols + [transcript_info['chrom'], str(converted_pos)]
            output_file.write('\t'.join(output) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert between transcript and genomic coordinates')
    parser.add_argument('-t', '--transcripts_input_p', default='inputs/transcripts_inputs.tsv', help='Path of the transcripts input file.')
    parser.add_argument('-q', '--queries_input_p', default='inputs/queries_inputs.tsv', help='Path of the queries input file.')
    parser.add_argument('-o', '--output_file_p', default='outputs/output.tsv', help='Path of the output file.')
    parser.add_argument('-r', '--reverse_orientation', action='store_true',
                        help="Optionally choose reverse orientation to indicate 3' to 5' instead of the default of 5' to 3'")
    parser.add_argument('-g', '--start_from_genomic', action='store_true',
                        help='Optionally convert from genomic coordinate to transcript coordinate instead of the default conversion of transcript to genomic')
    args = parser.parse_args()
    translate_coordinates(args.transcripts_input_p, args.queries_input_p, args.output_file_p, args.reverse_orientation, args.start_from_genomic)
