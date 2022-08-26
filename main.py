import numpy as np
from Bio.SeqIO import parse
from pickle import dump, load
from methods import collect_variant_counts, is_superset, is_subset, local_filter_seqs, adjust_letter, extend_template, generate_reference_freqs



from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument('-motifs', help='motifs file', type=str, required=True)
parser.add_argument('-ref', help='reference fasta', type=str, required=True)

parser.add_argument('-out', help='output file', type=str, default='output.txt')
parser.add_argument('-max_motifs', help='the maximum expected number of motifs extracted', default=20, type=int)
parser.add_argument('-min_conf', help='the minimal confidence value. Default is 1000', type=int, default=1000)
parser.add_argument('-ref_index', help='precalculatd dumped reference motifs file (will be created if not specified)', type=str, default=None)


args = parser.parse_args()

print()

print('PARAMETERS:')
for arg in vars(args):
     print(arg, ':\t', getattr(args, arg))
print()

f = parse(args.motifs, format='fasta')
seqs = [str(rec.seq) for rec in f]


ref_file = parse(args.ref, format='fasta')
for rec in ref_file:
    reference = str(rec.seq)
    break


N_REF = len(set(
    [reference[i:i+11] for i in range(len(reference) - 11)]
))

lengths = [4,5,6]



if args.ref_index is not None:
    with open(args.ref_index, 'rb') as fin:
        ref_motifs_counter = load(fin)

else:
    ref_motifs_counter, N_REF = generate_reference_freqs(reference, 11, lengths=lengths)
    with open('ref_counter.dump', 'wb') as fout:
        dump(ref_motifs_counter, fout)


new_seqs = seqs.copy()
seq_array = np.array([list(s) for s in new_seqs])


MOTIFS_SET = []
DETAILED_MOTIF_SET = []

print('ITERATION 1 ({} unexplained 11mers):'.format(len(seq_array)))


variants_counter_list = collect_variant_counts(seq_array, ref_motifs_counter, N_REF, lengths=lengths)


ITERATION = 2
while variants_counter_list[0][0] > args.min_conf:

    
    
    top_variant = variants_counter_list[0]


    extended_top_variant = extend_template(top_variant, maxlength=11)

    positions_to_adjust = []

    for i, pos in enumerate(extended_top_variant[2]):
        if extended_top_variant[1][i] == '.':
            positions_to_adjust.append((pos, i))
    
    print(extended_top_variant, positions_to_adjust)
    
    modifiable_extended_top_variant = [
        extended_top_variant[0],
        list(extended_top_variant[1]),
        list(extended_top_variant[2])
    ]

    
    for pos in positions_to_adjust:

        adjusted_pos_letter = adjust_letter(seq_array, extended_top_variant, pos[0], reference)
        modifiable_extended_top_variant[1][pos[1]] = adjusted_pos_letter

    extended_top_variant = (
        extended_top_variant[0],
        tuple(modifiable_extended_top_variant[1]),
        tuple(modifiable_extended_top_variant[2]),
    )

    print(extended_top_variant)

    is_superset_check = False
    is_subset_check = False

    for i, motif in enumerate(MOTIFS_SET):
        is_superset_check = is_superset(motif, ''.join(extended_top_variant[1]))
        is_subset_check = is_subset(motif, ''.join(extended_top_variant[1]))

        if is_subset_check:
            break

        if is_superset_check:
            break

    
    if is_superset_check == False:
        MOTIFS_SET.append(''.join(extended_top_variant[1]))
        DETAILED_MOTIF_SET.append(extended_top_variant)
    
    
    if len(MOTIFS_SET) == args.max_motifs:
        break

    else:
        print('{} already has a supermotif!'.format(extended_top_variant[1]))

    print(MOTIFS_SET)

    new_seqs = local_filter_seqs(new_seqs, extended_top_variant[2], extended_top_variant[1])


    seq_array = np.array([list(s) for s in new_seqs])
    
    print('ITERATION {} ({} unexplained 11mers):'.format(ITERATION, len(seq_array)))
    ITERATION += 1
    
    variants_counter_list = collect_variant_counts(seq_array, ref_motifs_counter, N_REF, lengths=lengths)

    
with open(args.out, 'w') as f:
    
    f.write('confidence\tmotif\tpos in 11mer(CORRECT!!!))\n')
    for rec in DETAILED_MOTIF_SET:
        f.write(
            '{}\t{}\t{}\n'.format(round(rec[0],2), rec[1], rec[2])
            )

print('Done!')