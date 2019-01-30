import argparse
import functools
import logging
from typing import NamedTuple

import gffutils

from general_gene_model import retrieve_segment_annotation

CACHE_MAX_SIZE = 2048

logger = logging.getLogger(__name__)


class Variant(NamedTuple):
    chrom: str
    start: int
    end: int
    tx_id: str
    classification: str


def find_overlapped_tx(chrom, start, end, db):
    # Search the best transcript id
    overlapped_txs = db.region(seqid=chrom, start=start, end=end, featuretype='transcript')
    overlapped_txs = [tx for tx in overlapped_txs if 'tag' in tx.attributes and 'basic' in tx.attributes['tag']]
    overlapped_txs = sorted(overlapped_txs, key=lambda tx: db.children_bp(tx), reverse=True)
    if overlapped_txs:
        return overlapped_txs[0]
    else:
        return None


def read_variant(pth, db):
    with open(pth) as f:
        columns = next(f)
        for i, line in enumerate(f, 1):
            record = line[:-1].split('\t')
            tx_id = record[33]
            chrom, start, end, classification = record[:4]
            if chrom == 'NA':
                # Use ICGC (PCAWG) part
                chrom, start, end, classification =  record[111:115]
                # Ad-hoc transcript annotation
                tx = find_overlapped_tx('chr' + chrom, int(start), int(end), db)
                if tx is None:
                    logger.warning(
                        f'Variant(chrom={chrom}, start={start}, end={end}) '
                        f'does not have any overlapped annotation'
                    )
                    yield None
                    continue
                tx_id = tx.id
            yield Variant(chrom, int(start), int(end), tx_id, classification)


def locate(variant, tx, segments_annotation):
    """Locate a variant in the general gene model."""
    if tx.strand == '-':
        segments_annotation = list(reversed(segments_annotation))

    try:
        segment, annotation = next(
            (segment, annotation) for segment, annotation in segments_annotation
            if variant.start >= segment.start and variant.end <= segment.end
        )
        if tx.strand == '+':
            segment_pos = (variant.start - segment.start + 1) / (segment.end - segment.start + 1)
        else:
            segment_pos = (segment.end - variant.end + 1) / (segment.end - segment.start + 1)
    except StopIteration:
        return None, None, None
    return segment, segment_pos, annotation


def calc_proximity(variant, tx, segment, annotation, segments_annotation):
    """Calculate the proximity to the exon boundary, exon1 and exon3 only."""
    if annotation == 'exon1':
        if tx.strand == '+':
            proximity = variant.start - segment.start
        else:
            proximity = segment.end - variant.end
    elif annotation == "5'UTR":
        exon1 = next(seg for seg, anno in segments_annotation if anno == 'exon1')
        if tx.strand == '+':
            proximity = variant.end - exon1.start
        else:
            proximity = exon1.end - variant.start
    elif annotation == 'exon3':
        if tx.strand == '+':
            proximity = segment.end - variant.end
        else:
            proximity = variant.start - segment.start
    elif annotation == "3'UTR":
        try:
            exon = next(seg for seg, anno in segments_annotation if anno == 'exon3')
        except StopIteration:
            exon = next(seg for seg, anno in segments_annotation if anno == 'exon1')
        if tx.strand == '+':
            proximity = exon.end - variant.start
        else:
            proximity = variant.end - exon.start
    else:
        proximity = 'NA'
    return proximity


def main(db_pth, variant_pth):
    # Load the GENCODE database
    logger.info(f'Gene model annotation from {db_pth}')
    db = gffutils.FeatureDB(db_pth)
    db.execute('PRAGMA cache_size=-16000000')
    db.execute('PRAGMA temp_store=MEMORY')

    # Set up the variant reader
    logger.info(f'Reading variants from {variant_pth}')
    variant_reader = read_variant(variant_pth, db)

    # Create a transcript ID resolver
    tx_id_to_full = {}
    for row in db.execute("SELECT id FROM features WHERE featuretype='transcript'"):
        full_tx_id = row[0]
        tx_id = full_tx_id.rsplit('.', 1)[0]
        tx_id_to_full[tx_id] = full_tx_id
    logger.info(f'Detected {len(tx_id_to_full):,d} transcript IDs')

    # A cached annotation retriever
    cached_annotate = functools.lru_cache(maxsize=CACHE_MAX_SIZE)(
        functools.partial(retrieve_segment_annotation, db=db)
    )

    HEADER = ('Annotation', 'Normalized position', 'Proximity', 'Full Transcript ID')
    EMPTY_OUTPUT = ('NA', 'NA', 'NA')

    print(*HEADER, sep='\t')
    for i, variant in enumerate(variant_reader, 1):
        if i % 5000 == 0:
            logger.info(f'Processed {i:,d} variants')

        # Skip unannotated variants
        if variant is None:
            print(*EMPTY_OUTPUT, 'NA', sep='\t')
            continue

        if '.' not in variant.tx_id:
            full_tx_id = tx_id_to_full[variant.tx_id]
        else:
            full_tx_id = variant.tx_id
        tx = db[full_tx_id]
        segments_annotation = cached_annotate(tx)
        segment, normalized_pos, annotation = locate(variant, tx, segments_annotation)
        if annotation is None:
            print(*EMPTY_OUTPUT, full_tx_id, sep='\t')
            continue
        proximity = calc_proximity(variant, tx, segment, annotation, segments_annotation)
        print(annotation, normalized_pos, proximity, full_tx_id, sep='\t')

    logger.info(f'Successfully annotated {i:,d} variants')
    logger.info(f'Annotation cache usage: {cached_annotate.cache_info()}')


def setup_cli():
    # Setup console logging
    console = logging.StreamHandler()
    all_loggers = logging.getLogger()
    all_loggers.setLevel(logging.INFO)
    all_loggers.addHandler(console)
    log_fmt = '[%(asctime)s][%(levelname)-7s] %(message)s'
    log_formatter = logging.Formatter(log_fmt, '%Y-%m-%d %H:%M:%S')
    console.setFormatter(log_formatter)

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        'db_pth', help="Path to the GENCODE database made by gffutils"
    )
    parser.add_argument(
        'variant_pth', help="Path to the variant TSV file."
    )
    return parser


if __name__ == '__main__':
    parser = setup_cli()
    args = parser.parse_args()

    main(args.db_pth, args.variant_pth)

