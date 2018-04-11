.mode csv
.separator \t
-- Read in exome headers
.import processed_data/exome.maf.for_db.header exome
.import processed_data/exome.pipe exome

-- Read in genome
.import processed_data/genome.maf.header genome
.import processed_data/genome.pipe genome

-- Read in genome aliquot mapping
.import /diskmnt/Projects/ICGC_MC3/ID_Mapping/SQL_Mapping.txt genome_aliquot
