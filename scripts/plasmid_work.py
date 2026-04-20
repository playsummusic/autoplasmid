from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from .gibson_design import OVS

def _resolve_primer_location(record_seq, insert_start, primer):
    full_seq = primer.get('full_seq', primer['binding_seq'])
    if primer['binding_length'] is None:
        return None, None, None
    binding_len = int(primer['binding_length'])

    overhang_len = len(full_seq) - binding_len

    binding_start_abs = insert_start + int(primer['binding_start'])
    binding_end_abs = insert_start + int(primer['binding_end'])

    if primer['dir'] == 'fwd':
        expected_start = max(0, binding_start_abs - overhang_len)
        expected_end = binding_end_abs
        search_seq = full_seq
        strand = 1
    else:
        expected_start = binding_start_abs
        expected_end = min(len(record_seq), binding_end_abs + overhang_len)
        search_seq = str(Seq(full_seq).reverse_complement())
        strand = -1

    found_at = str(record_seq).find(search_seq)
    if found_at != -1:
        return found_at, found_at + len(search_seq), strand

    return expected_start, expected_end, strand

def add_primers_to_construct(record, insert_start, primers_list):
    """
    Fügt Primer-Features ins Plasmid ein.
    
    Args:
        record: SeqRecord des Plasmids
        insert_start: Position, wo der Insert im Plasmid anfängt
        primers_list: Liste von Primer-Dicts aus gibson_design.design_binding_primers()
    """
    for primer in primers_list:
        plasmid_start, plasmid_end, strand = _resolve_primer_location(
            record.seq,
            insert_start,
            primer
        )
        if plasmid_start is None:
            continue
        
        primer_feature = SeqFeature(
            FeatureLocation(plasmid_start, plasmid_end, strand=strand),
            type="primer_bind",
            qualifiers={
                "label": f"{primer['primer_name']}",
                "sequence": primer['full_seq'],
                "tm": [str(round(primer['mt'], 1))],
                "gc": [str(round(primer['gc'], 1))]
            }
        )
        record.features.append(primer_feature)
    
    return record

def insert_into_vector(
    vector_gb_path,
    insert_seq,
    delete_start,
    delete_end,
    output_path,
    construct_name,
    insert_name="INSERT",
    primers=None,
    features_to_draw=None,
    mode='pppt4'
):
    # 1️⃣ Vector laden
    record = SeqIO.read(vector_gb_path, "genbank")

    # Locus / Accession anpassen
    locus_name = construct_name[:16] if len(construct_name) > 16 else construct_name
    record.id = locus_name
    record.name = locus_name
    record.description = "."
    record.annotations["accessions"] = [construct_name]

    original_seq = record.seq
    insert_seq = Seq(insert_seq)

    # 2️⃣ Bereich löschen und Insert einfügen
    insert_len = len(insert_seq)
    deleted_len = delete_end - delete_start

    new_seq = original_seq[:delete_start] + insert_seq + original_seq[delete_end:]
    record.seq = new_seq

    # 3️⃣ Features anpassen
    shift = insert_len - deleted_len
    new_features = []

    for feature in record.features:
        start = int(feature.location.start)
        end = int(feature.location.end)

        # Feature komplett hinter dem ersetzten Bereich → verschieben
        if start >= delete_end:
            new_location = FeatureLocation(
                start + shift,
                end + shift,
                strand=feature.location.strand
            )
            feature.location = new_location

        # Feature komplett im gelöschten Bereich → überspringen (löschen)
        elif start >= delete_start and end <= delete_end:
            continue

        # Feature überlappt den Bereich → nicht automatisch handhaben
        elif (start < delete_end and end > delete_start):
            continue

        new_features.append(feature)

    record.features = new_features

    # 4️⃣ Insert-Feature hinzufügen
    insert_feature = SeqFeature(
        FeatureLocation(delete_start, delete_start + insert_len),
        type="CDS",
        qualifiers={"label": insert_name}
    )

    record.features.append(insert_feature)

    # 5️⃣ Primer hinzufügen (optional)
    if primers:
        add_primers_to_construct(record, delete_start, primers)

    # 6️⃣ Zusätzliche Features hinzufügen (optional)
    if features_to_draw:
        for feature_def in features_to_draw:
            rel_start = int(feature_def.get("start", 0))
            rel_end = int(feature_def.get("end", 0))
            if rel_end <= rel_start:
                continue

            if feature_def.get("relative_to_insert", True):
                start = delete_start + rel_start - len(OVS[mode]['FWD_OV'])
                end = delete_start + rel_end + len(OVS[mode]['RV_OV'])
                if mode == 'pppt4':
                    start += 3
            else:
                start = rel_start
                end = rel_end

            strand = feature_def.get("strand")
            feature_type = feature_def.get("type", "misc_feature")
            qualifiers = dict(feature_def.get("qualifiers", {}))
            label = feature_def.get("label")
            if label:
                qualifiers.setdefault("label", label)
            record.features.append(
                SeqFeature(
                    FeatureLocation(start, end, strand=strand),
                    type=feature_type,
                    qualifiers=qualifiers,
                )
            )

    # 7️⃣ Speichern
    SeqIO.write(record, output_path, "genbank")

