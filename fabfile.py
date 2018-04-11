from datetime import datetime
from pathlib import Path
from fabric.api import local, task, lcd, env
from fabric.contrib.console import confirm
from fabric.utils import abort, puts

FABFILE_DIR = str(Path(env.real_fabfile).parent)
TARGET_SNAKEFILE = Path(
    '/gscmnt/gc7210/dinglab/medseq/TCGA/MC3/mbailey/ICGC_MC3/Snakefile'
)
GENERATED_OUTPUT_DIRS = [
    Path(d) for d in [
        'annotations',
        'processed_data',
    ]
]


@task
def update_snakefile():
    """Copy Bailey's Snakefile to local directory"""
    if not TARGET_SNAKEFILE.exists():
        abort('Cannot find the target Snakefile at %s' % TARGET_SNAKEFILE)
    old_infix = '{:%Y-%m-%d_%H:%M}'.format(datetime.now())
    with lcd(FABFILE_DIR):
        if Path('Snakefile.copied').exists():
            local(
                'cp Snakefile.copied Snakefile.copied.{}'.format(old_infix)
            )
        local(
            'cp {} Snakefile.copied'.format(TARGET_SNAKEFILE)
        )


@task
def gen_sample_names():
    """Generate all exome and genome samples names"""
    with lcd(FABFILE_DIR):
        local(
            r"awk -F'\t' '{print $12}' "
            r"processed_data/exome.broadbed.gaf4bed.exon.maf | "
            "sort -u > exome_samples.list"
        )
        local(
            r"awk -F'\t' '{print $46}' "
            r"processed_data/genome.broadbed.gaf4bed.exon.maf | "
            "sort -u > genome_samples.list"
        )


@task
def clean():
    """Remove all generated outputs."""
    dirs_to_remove = [
        '%s/' % folder for folder in GENERATED_OUTPUT_DIRS
    ]
    puts(
        'Fabric is going to remove all generated outputs '
        'under folders: {!s}'.format(', '.join(dirs_to_remove))
    )
    confirm('Proceed to remove all outputs in these folders?', False)

    with lcd(FABFILE_DIR):
        for folder in GENERATED_OUTPUT_DIRS:
            local('rm -rf %s' % folder)
