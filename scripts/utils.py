import re
import os
import datetime
EXT_PER_BP = 2000
PRIMER_START_NUM = 14
PLASMID_START_NUM = 11
REACTION_START_NUM = 0
EXPORT_FOLDER_NAME = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
MIN_SEQ_LEN_FOR_PRIMER = 125
os.makedirs(f"output/{EXPORT_FOLDER_NAME}", exist_ok=True)

FWD_RV_TRANSLATE = {
    'fwd': "5'-3'",
    'rv': "3'-5'"
}

gibson_design_args = {
        'pppt4': {
            'vector_path': 'input/pppt4_bb.gb',
            'delete_start': 1102,
            'delete_end': 1121,
            'output_dir': f'output/{EXPORT_FOLDER_NAME}/plasmids',
            'prefix': 'pPpT4',
            'excel_output': f'output/{EXPORT_FOLDER_NAME}/gibson_primers_pppt4.xlsx',
            'ATG_add': 'ATG',
            'ATG_cut': 3,
        },
        'pppt4a': {
            'vector_path': 'input/pppt4alpha_bb.gb',
            'delete_start': 1351,
            'delete_end': 1371,
            'output_dir': f'output/{EXPORT_FOLDER_NAME}/plasmids',
            'prefix': 'pPpT4a',
            'excel_output': f'output/{EXPORT_FOLDER_NAME}/gibson_primers_pppt4a.xlsx',
            'ATG_add': '',
            'ATG_cut': 0
        }
    }

def safe_string(value):
    value = str(value)
    value = re.sub(r'[\\/:*?"<>|]+', '_', value)
    value = re.sub(r'\s+', '_', value)
    value = re.sub(r'_+', '_', value)
    return value.strip('._')

