
DROP TABLE full;
CREATE TABLE full AS
SELECT e.*,g.*,a.mc3_exome_barcode FROM exome AS e
INNER JOIN genome_aliquot AS a ON e.Tumor_Sample_Barcode = a.mc3_exome_barcode
LEFT JOIN genome AS g ON a.pcawg_donor_id = g.Donor_ID
  AND e.Chromosome = g.Chromosome
  AND e.End_Position >= g.Start_position
  AND e.Start_Position <= g.End_position

UNION ALL

SELECT e.*,g.*,a.mc3_exome_barcode FROM genome AS g
INNER JOIN genome_aliquot AS a ON g.Donor_ID = a.pcawg_donor_id
LEFT JOIN exome AS e ON e.Tumor_Sample_Barcode = a.mc3_exome_barcode
  AND e.Chromosome = g.Chromosome
  AND e.End_Position >= g.Start_position
  AND e.Start_Position <= g.End_position
WHERE e.Tumor_Sample_Barcode IS NULL;

.headers on
.separator "\t"
.output full.tsv
.nullvalue "NA"
SELECT * FROM full
EOF;
