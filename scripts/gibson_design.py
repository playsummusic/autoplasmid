import primer3
import pandas as pd
from .utils import gibson_design_args
OVS = {
    'pppt4':{
    'FWD_OV': "TCAAAAAACAACTAATTATTGAAAGAATTAAAACGATG".upper(),
    'RV_OV': "gtggatgtgaccaagcagagcggcc".upper(),
    'RV_OV_COMP': "ggccgctctgcttggtcacatccac".upper()
    },
    'pppt4a':{
    'FWD_OV': "gctaaggaagagggtgtctctctcgagaagagagaggccgaagct".upper(),
    'RV_OV': "gtggatgtgaccaagcagagcggcc".upper(),
    "RV_OV_COMP": "ggccgctctgcttggtcacatccac".upper()
    }
}


def design_binding_primers(insert_seq, seq_id=0, mode='pppt4', forward_only=False):

    insert_seq = insert_seq[gibson_design_args[mode]['ATG_cut']:]
    global_args = {
        'PRIMER_MIN_SIZE': 10,
        'PRIMER_MAX_SIZE': 36,
        'PRIMER_OPT_TM': 61.0,
        'PRIMER_MIN_TM': 56.0,
        'PRIMER_MAX_TM': 65.0,
        'PRIMER_PRODUCT_SIZE_RANGE': [[len(insert_seq), len(insert_seq)]],
        'PRIMER_MAX_HAIRPIN_TH': 100.0,
        'PRIMER_MAX_SELF_ANY_TH': 100.0,
        'PRIMER_MAX_SELF_END_TH': 100.0,
        'PRIMER_MIN_GC': 15,
        'PRIMER_MAX_GC': 85,
        'PRIMER_MAX_POLY_X': 8
    }
    seq_args = {
        'SEQUENCE_ID': 'insert',
        'SEQUENCE_TEMPLATE': insert_seq,
        'SEQUENCE_FORCE_LEFT_START': 0,
        'SEQUENCE_FORCE_RIGHT_START': len(insert_seq)-1
    }
    result = primer3.bindings.design_primers(seq_args, global_args)

    # Extract primer data - try up to 5 primer pairs
    primers_data = []
    
    # Try PRIMER_LEFT_0 through PRIMER_LEFT_4
    for PRIMER_LFNR in range(5):
        try:
            left_pos, left_len = result[f'PRIMER_LEFT_{PRIMER_LFNR}']
            primers_data.append({
                'binding_seq': result[f'PRIMER_LEFT_{PRIMER_LFNR}_SEQUENCE'],
                'dir': "fwd",
                'gc': round(result[f'PRIMER_LEFT_{PRIMER_LFNR}_GC_PERCENT'], 1),
                'mt': round(result[f'PRIMER_LEFT_{PRIMER_LFNR}_TM'], 1),
                'binding_length': len(result[f'PRIMER_LEFT_{PRIMER_LFNR}_SEQUENCE']),
                'binding_start': left_pos,  # Position im Insert
                'binding_end': left_pos + left_len
            })
            if forward_only is False:
                # Reverse primer (PRIMER_RIGHT_N)
                right_pos, right_len = result[f'PRIMER_RIGHT_{PRIMER_LFNR}']
                right_start = right_pos - right_len + 1
                primers_data.append({
                    'binding_seq': result[f'PRIMER_RIGHT_{PRIMER_LFNR}_SEQUENCE'],
                    'dir': "rv",
                    'gc': round(result[f'PRIMER_RIGHT_{PRIMER_LFNR}_GC_PERCENT'], 1),
                    'mt': round(result[f'PRIMER_RIGHT_{PRIMER_LFNR}_TM'], 1),
                    'binding_length': len(result[f'PRIMER_RIGHT_{PRIMER_LFNR}_SEQUENCE']),
                    'binding_start': right_start,  # Position im Insert (0-based)
                    'binding_end': right_pos + 1
                })
            
            # Success - return primers_data
            return primers_data
            
        except KeyError:
            # This primer pair doesn't exist, try next one
            continue
    # If we get here, none of the 5 primer pairs worked
    print("-"*50)
    print(f"Couldnt design primers for sequence ID {seq_id}")
    # Get error explanation from primer3 result
    if 'PRIMER_ERROR' in result:
        print(f"  Error: {result['PRIMER_ERROR']}")
    if 'PRIMER_PAIR_EXPLAIN' in result:
        print(f"  Pair Explain: {result['PRIMER_PAIR_EXPLAIN']}")
    if 'PRIMER_LEFT_EXPLAIN' in result:
        print(f"  Left Explain: {result['PRIMER_LEFT_EXPLAIN']}")
    if 'PRIMER_RIGHT_EXPLAIN' in result:
        print(f"  Right Explain: {result['PRIMER_RIGHT_EXPLAIN']}")
    print("-"*50)
    primers_fallback = [
        {
            'binding_seq': None,
            'dir': "fwd",
            'gc': None,
            'mt': None,
            'binding_length': None,
            'binding_start': None,  # Position im Insert (0-based)
            'binding_end': None
        }
    ]
    if not forward_only:
        primers_fallback.append(
            {
                'binding_seq': None,
                'dir': "rv",
                'gc': None,
                'mt': None,
                'binding_length': None,
                'binding_start': None,  # Position im Insert (0-based)
                'binding_end': None
            }
        )
    return primers_fallback

def design_gibson_overhangs(primer, mode='pppt4'):
    # Return early if no binding sequence was found
    if primer['binding_seq'] is None:
        primer['full_seq'] = None
        primer['ov_length'] = None
        primer['full_length'] = None
        return primer
    
    if primer['dir'] == "fwd":
        primer['full_seq'] = f"{OVS[mode]['FWD_OV']}{primer['binding_seq']}"
        primer['ov_length'] = len(OVS[mode]['FWD_OV'])

    elif primer['dir'] == "rv":
        primer['full_seq'] = f"{OVS[mode]['RV_OV']}{primer['binding_seq']}"
        primer['ov_length'] = len(OVS[mode]['RV_OV'])
    
    
    if len(primer['full_seq']) > 60:
        diff = len(primer['full_seq']) - 60
        ov_len = primer['ov_length'] - diff
        if ov_len >= 36 and mode == 'pppt4a':
            primer['full_seq'] = primer['full_seq'][diff:]
            primer['ov_length'] -= diff
        elif mode == 'pppt4':
            primer['full_seq'] = primer['full_seq'][diff:]
            primer['ov_length'] -= diff 
    primer['full_length'] = len(primer['full_seq'])
    return primer

# Example usage
if __name__ == '__main__':
    inserts = [
        'ATGCATCATCACCATCACCATGCTCACATCCACGACCTGGCGCCGGAAGTAAGCAACTACTCTTCTGGTCGCCTGACCCCGCCGACCCCAGTTAGGTTCCCGCGCACCCCAGTGTTCGCATCTATGAACAAACCGTGCCGCTTCGAAGGTGACGTTTTCGACCTGGAAGTTTCTGGTGCTATCCCGCCGGACATCGACGGTACCTTCTTCCGCGTTCAGCCGGACCACCGCTTCCCGCCGCTGTTCGAAGACGACATCCACTTCAACGGTGACGGTTCTGTTACCGCTATCCGCATCTCTGGTGGTCACGCTGACCTGCGCCAGCGCTACGTTCGCACCGAACGCTACCTGCTGGAAACCCGCGCTCGCCGCTCTCTGTTCGGTCGCTACCGCAACCCGTGGACCGACAACGAATCTGTTCGCGGTGTTATCCGCACCGCTTCTAACACCAACGTTGTTTTCTGGCGCGGTGCTCTGCTGGCTATGAAAGAAGACGGTCCGCCGTTCGCTATGGACCCGGTTACCCTGGAAACCCTGGGTCGCTACGACTTCGAAGGTCAGATCCTGTCTCCGACCTTCACCGCTCACCCGAAAATCGACCCGGACACCGGTGAAATGGTTTGCTTCGCTTACGAAGCTGGTGGTGACGGTTCTGACTGCTCTGTTGACGTTGCTGTTTGGACCGTTGACGCTGACGGTAAAAAAGTTGAAGAATGCTGGTACAAAGCTCCGTTCGCTGGTATGATCCACGACTGCGGTATCACCAAAAACTGGGTTGTTCTGCCGCTGACCCCGATCAAAATGGACCTGGAACGCATGAAACGCGGTGGTAACAAATTCGCTTGGGACCCGTCTGAAGACCAGTGGTACGGTGTTGTTCCGCGCCGCGGTGCTAAATCTGACGACATCATCTGGTTCCGCGCTGACAACGGCTTCCACGGTCACGTTGCTGGTTGCTACGAACTGCCGTCTGGTGAAATCGTTTTCGACCTGACAGTTGCGGACGGCAACGTCTTCTTCTTCTTCCCGCCGGACGACAACATCACCCCGCCGGCTGACGGTGTTGCTAAACGCAACCGCCTGTCTTCTCCGACCGTTCGCTGGATCTTCGACCCGAAAGCTAAAAAATCTGCTATCCGCACCGAAGCTGCTGGTGACGCTGACATCTGGGTTGCTGACGAACGCGTTAAACCGGCTCTGACCTGGCCGACCAACGGTGAATTCTCTCGCATCGACGACCGCTACGTTACCAAACCGTACCGCCACTTCTGGCAGGCTGTTGTTGACCCGACCCGCCCGTACGACTTCGAAAAATGCGGTCCGCCGGCTGGTGGTCTGTTCAACTGCCTGGGTCACTACACCTGGTCTGACCAGAACTACCACCACGGTCACAACACCGGTGACCCGTCTGGTGACGGTCGCTCTAACGGTTCTGCTGAAGAAGCTACCGCTGGTAAATTCGGCCTGCAGGACGTATACTTCGCCGGTCCGACCATGACCTTCCAGGAACCGACCTTCATCCCGCGCCAGGGTGCTGCTGAAGGTGAAGGTTACCTGATCGCTCTGCTGAACCACCTGGACGAACTGCGCAACGACGTTGTTATCTTCGAAGCTCGCAACCTGGGTAAAGGTCCGCTGGCTGTTATCCACCTGCCGCTGAAACTGAAACTGGGTCTGCACGGTAACTGGGTTGACTCTCGCGAAATCGAAGCTTGGCGCCGCCGCCGCGCTGAAAACGGTGACGTTGGTCCGCTGCGCGTTGCTAAAGAACCGCTGCCGTGGCAGAAAAAATTCGCTGCTGCTGCTCAGAACGGTTCTAACGGTGTT',
    ]
    
    all_primers = []
    
    for i, insert_seq in enumerate(inserts, 1):
        primers = design_binding_primers(insert_seq)
        for primer in primers:
            primer['Nr'] = i
            primer = design_gibson_overhangs(primer=primer)
            all_primers.append(primer)
    #print(all_primers)
    #primers_df = pd.DataFrame(all_primers)
