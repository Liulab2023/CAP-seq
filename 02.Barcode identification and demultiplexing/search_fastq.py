import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import multiprocessing as mp


# set parameter parser
parser = argparse.ArgumentParser(description='Find a specific sequence in FASTQ files')
parser.add_argument('-i', '--fastq_folder', required=True, type=str, help='Path to the folder containing FASTQ files')
parser.add_argument('-o', '--output_folder', required=True, help='Output fastq file for matched sequences')
parser.add_argument('-s', '--specific_seq', required=True, type=str, help='The specific sequence to be found in the FASTQ files')
parser.add_argument('-t', '--threshold', type=int, default=20, help='The threshold score for matching sequences (default=20)')
args = parser.parse_args()

if not os.path.exists(args.output_folder):
    os.makedirs(args.output_folder)

# def function: reverse complement
def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


# def function：Smith-Waterman align
def smith_waterman(seq1, seq2):
    # initialize aligner
    aligner = PairwiseAligner()
    aligner.mode = 'local'  # local alignment
    aligner.match_score = 2  # score for match
    aligner.mismatch_score = -1  # score for mismatch
    aligner.open_gap_score = -1  # score for open gap
    aligner.extend_gap_score = -1  # score for extend gap

    # align
    alignments = aligner.align(seq1, seq2)

    # default for not align
    if not alignments:
        return [None, None, 0, 0, None, None]

    # sort alignments
    sorted_alignments = sorted(alignments, key=lambda x: -x.score)

    # get the best alignments
    top_alignment = sorted_alignments[0]
    top_score = top_alignment.score
    top_aln_start = top_alignment.aligned[0][0][0]  # start position
    top_aln_end = top_alignment.aligned[0][-1][-1]  # end position

    # second best alignment
    try:
        second_top_score = sorted_alignments[1].score
    except IndexError:
        second_top_score = 0

    # get aligned seq
    seq1_match = str(top_alignment.target)
    seq2_match = str(top_alignment.query)

    return [seq1_match, seq2_match, int(top_score), int(second_top_score), int(top_aln_start), int(top_aln_end)]

def process_fastq(file_path, specific_seq, threshold, process_id, output_folder):
    # creat temp folder for process
    temp_folder = os.path.join(output_folder, f"process_{process_id}")
    os.makedirs(temp_folder, exist_ok=True)
    output_path_match1 = os.path.join(temp_folder, "output_match1.fastq")
    output_path_match2 = os.path.join(temp_folder, "output_match2.fastq")

    total_seq = 0
    matched_seq = 0
    re_matched_seq = 0
    complementary_matched_seq = 0
    multiple_matched_seq = 0
    fail_matched_seq = 0
    results = []

    for record in SeqIO.parse(file_path, "fastq"):
        seq = str(record.seq)

        # match the forward sequence
        top_aln = smith_waterman(seq, specific_seq)
        # match the reverse complement sequence
        rc_top_aln = smith_waterman(seq, reverse_complement(specific_seq))

        # Determine whether the match was successful
        matched = 0 #1 for forward match, 2 for reverse match, 3 for complementary match, 4 for multiple match, 5 for fail match
        if top_aln[2] >= threshold and rc_top_aln[2] >= threshold:
            complementary_matched_seq += 1
        elif top_aln[3] >= threshold or rc_top_aln[3] >= threshold:
            multiple_matched_seq += 1
        elif top_aln[2] < threshold and rc_top_aln[2] < threshold:
            fail_matched_seq += 1
        elif top_aln[2] >= threshold:
            matched_seq += 1
            matched = 1
            sub_seq = top_aln[0][top_aln[4]-29:top_aln[4]]
            results.append((record, matched, sub_seq))

        elif rc_top_aln[2] >= threshold:
            re_matched_seq += 1
            matched = 2
            sub_seq = rc_top_aln[0][rc_top_aln[5]:rc_top_aln[5]+29]
            results.append((record, matched, sub_seq))
        
        total_seq += 1
    
    with open(output_path_match1, "a") as out_file1, open(output_path_match2, "a") as out_file2:
        for record, match, sub_seq in results:
            if match == 1:
                write_output(record, sub_seq, out_file1)
            elif match == 2:
                write_output(record, sub_seq, out_file2)

    return [total_seq, matched_seq, re_matched_seq,complementary_matched_seq, multiple_matched_seq, fail_matched_seq, temp_folder]

def write_output(record, sub_seq, out_file):
    seq = str(record.seq)
    out_file.write('@' + record.id + '_'+str(sub_seq)+'\n')
    out_file.write(seq + '\n')
    out_file.write('+' + '\n')
    quality_string = ''.join([chr(score + 33) for score in record.letter_annotations['phred_quality']])
    out_file.write(quality_string + '\n')
      
def main():
    fastq_files = [os.path.join(args.fastq_folder, file_name) for file_name in os.listdir(args.fastq_folder) if file_name.endswith('.fastq')]

    pool = mp.Pool(mp.cpu_count())

    # When using starmap, pass all the necessary parameters
    args_for_starmap = [(file_path, args.specific_seq, args.threshold, process_id, args.output_folder) for process_id, file_path in enumerate(fastq_files)]
    results = pool.starmap(process_fastq, args_for_starmap)

    # Merge temporary files and clean up
    #temp_folders = [os.path.join(args.output_folder, f"process_{pid}") for pid in range(len(fastq_files))]
    #final_output_paths = [os.path.join(args.output_folder, "output_match1.fastq"), os.path.join(args.output_folder, "output_match2.fastq")]
    #merge_and_clean(temp_folders, final_output_paths)

    # Output statistical results
    total_seq = sum(res[0] for res in results)
    matched_seq = sum(res[1] for res in results)
    re_matched_seq = sum(res[2] for res in results)
    complementary_matched_seq = sum(res[3] for res in results)
    multiple_matched_seq = sum(res[4] for res in results)
    fail_matched_seq = sum(res[5] for res in results)
    parent_folder = os.path.dirname(args.output_folder)
    output_file_path = os.path.join(parent_folder, "statistical.txt")
    with open(output_file_path, "w") as f:
        f.write('Total reads: ' + str(total_seq) + '\n')
        f.write('Percent for single forward match: ' + str(matched_seq / total_seq) + '\n')
        f.write('Percent for single reverse match: ' + str(re_matched_seq / total_seq) + '\n')
        f.write('Percent for complementary match: ' + str(complementary_matched_seq / total_seq) + '\n')
        f.write('Percent for multiple match: 多' + str(multiple_matched_seq / total_seq) + '\n')
        f.write('Percent for failed match: ' + str(fail_matched_seq / total_seq) + '\n')

if __name__ == "__main__":
    main()


