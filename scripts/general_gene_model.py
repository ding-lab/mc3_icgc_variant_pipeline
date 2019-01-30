import logging

logger = logging.getLogger(__name__)


def get_gene_model_segments(tx, db):
    if tx.strand == '+':
        gene_model_segments = list(
            db.children(tx, order_by='start', featuretype=['UTR', 'CDS']))
    else:
        gene_model_segments = list(
            db.children(tx, order_by='end', reverse=True, featuretype=['UTR', 'CDS'])
        )
    return gene_model_segments


def annotate(tx, gene_model_segments):
    """Anntatote the current gene model using the generalized gene model"""
    segments_annotation = []
    segment_iter = iter(gene_model_segments)
    cur_segment = next(segment_iter)

    try:
        while cur_segment.featuretype != 'CDS':
            segments_annotation.append("5'UTR")
            cur_segment = next(segment_iter)

        # There will be no exon2 if total number of CDS <= 2
        num_cds = sum(1 for seg in gene_model_segments if seg.featuretype == 'CDS')
        if num_cds >= 2:
            # Annotate all CDS except for the last one with exon2
            segments_annotation.append('exon1')
            cur_segment = next(segment_iter)
            for _ in range(num_cds - 2):
                segments_annotation.append('exon2')
                cur_segment = next(segment_iter)
            segments_annotation.append('exon3')
            cur_segment = next(segment_iter)
        elif num_cds == 1:
            segments_annotation.append('exon1')
            cur_segment = next(segment_iter)
        else:
            logger.error(f'Transcript {tx.id} has no CDS.')
            raise ValueError(f'Transcript {tx.id} has no CDS.')

        while cur_segment.featuretype == 'UTR':
            segments_annotation.append("3'UTR")
            cur_segment = next(segment_iter)
    except StopIteration:
        pass

    return segments_annotation


def retrieve_segment_annotation(tx_id, db):
    """Retrieve the annotated gene model segments using general gene model."""
    tx = db[tx_id]
    gene_model_segments = get_gene_model_segments(tx, db)

    if gene_model_segments:
        segments_annotation = annotate(tx, gene_model_segments)
    else:
        # logger.warning(f'Transcript {tx.id} does not have any gene model segment')
        segments_annotation = []

    return list(zip(gene_model_segments, segments_annotation))
