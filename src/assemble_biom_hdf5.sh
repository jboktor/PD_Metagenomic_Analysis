#!/bin/bash



### Merged Data
# Species
biom convert -i data/Nestedness/biom/Merged/species.tsv \
-o data/Nestedness/biom/Merged/species_hdf5.biom \
--table-type="OTU table" --to-hdf5
biom add-metadata -i data/Nestedness/biom/Merged/species_hdf5.biom \
-o data/Nestedness/biom/Merged/species_hdf5_meta.biom \
--sample-metadata-fp data/Nestedness/biom/Merged/metadata.tsv \
--sample-header donor_id,donor_group,PD
# Genera
biom convert -i data/Nestedness/biom/Merged/genera.tsv \
-o data/Nestedness/biom/Merged/genera_hdf5.biom \
--table-type="OTU table" --to-hdf5
biom add-metadata -i data/Nestedness/biom/Merged/genera_hdf5.biom \
-o data/Nestedness/biom/Merged/genera_hdf5_meta.biom \
--sample-metadata-fp data/Nestedness/biom/Merged/metadata.tsv \
--sample-header donor_id,donor_group,PD
# Family
biom convert -i data/Nestedness/biom/Merged/family.tsv \
-o data/Nestedness/biom/Merged/family_hdf5.biom \
--table-type="OTU table" --to-hdf5
biom add-metadata -i data/Nestedness/biom/Merged/family_hdf5.biom \
-o data/Nestedness/biom/Merged/family_hdf5_meta.biom \
--sample-metadata-fp data/Nestedness/biom/Merged/metadata.tsv \
--sample-header donor_id,donor_group,PD
# Order
biom convert -i data/Nestedness/biom/Merged/order.tsv \
-o data/Nestedness/biom/Merged/order_hdf5.biom \
--table-type="OTU table" --to-hdf5
biom add-metadata -i data/Nestedness/biom/Merged/order_hdf5.biom \
-o data/Nestedness/biom/Merged/order_hdf5_meta.biom \
--sample-metadata-fp data/Nestedness/biom/Merged/metadata.tsv \
--sample-header donor_id,donor_group,PD
# Class
biom convert -i data/Nestedness/biom/Merged/class.tsv \
-o data/Nestedness/biom/Merged/class_hdf5.biom \
--table-type="OTU table" --to-hdf5
biom add-metadata -i data/Nestedness/biom/Merged/class_hdf5.biom \
-o data/Nestedness/biom/Merged/class_hdf5_meta.biom \
--sample-metadata-fp data/Nestedness/biom/Merged/metadata.tsv \
--sample-header donor_id,donor_group,PD
# Phylum
biom convert -i data/Nestedness/biom/Merged/phylum.tsv \
-o data/Nestedness/biom/Merged/phylum_hdf5.biom \
--table-type="OTU table" --to-hdf5
biom add-metadata -i data/Nestedness/biom/Merged/phylum_hdf5.biom \
-o data/Nestedness/biom/Merged/phylum_hdf5_meta.biom \
--sample-metadata-fp data/Nestedness/biom/Merged/metadata.tsv \
--sample-header donor_id,donor_group,PD

#--------------------------------------------------------

### TBC Data
# Species
biom convert -i data/Nestedness/biom/TBC/species.tsv \
-o data/Nestedness/biom/TBC/species_hdf5.biom \
--table-type="OTU table" --to-hdf5
biom add-metadata -i data/Nestedness/biom/TBC/species_hdf5.biom \
-o data/Nestedness/biom/TBC/species_hdf5_meta.biom \
--sample-metadata-fp data/Nestedness/biom/TBC/metadata.tsv \
--sample-header donor_id,donor_group,PD
# Genera
biom convert -i data/Nestedness/biom/TBC/genera.tsv \
-o data/Nestedness/biom/TBC/genera_hdf5.biom \
--table-type="OTU table" --to-hdf5
biom add-metadata -i data/Nestedness/biom/TBC/genera_hdf5.biom \
-o data/Nestedness/biom/TBC/genera_hdf5_meta.biom \
--sample-metadata-fp data/Nestedness/biom/TBC/metadata.tsv \
--sample-header donor_id,donor_group,PD
# Family
biom convert -i data/Nestedness/biom/TBC/family.tsv \
-o data/Nestedness/biom/TBC/family_hdf5.biom \
--table-type="OTU table" --to-hdf5
biom add-metadata -i data/Nestedness/biom/TBC/family_hdf5.biom \
-o data/Nestedness/biom/TBC/family_hdf5_meta.biom \
--sample-metadata-fp data/Nestedness/biom/TBC/metadata.tsv \
--sample-header donor_id,donor_group,PD
# Order
biom convert -i data/Nestedness/biom/TBC/order.tsv \
-o data/Nestedness/biom/TBC/order_hdf5.biom \
--table-type="OTU table" --to-hdf5
biom add-metadata -i data/Nestedness/biom/TBC/order_hdf5.biom \
-o data/Nestedness/biom/TBC/order_hdf5_meta.biom \
--sample-metadata-fp data/Nestedness/biom/TBC/metadata.tsv \
--sample-header donor_id,donor_group,PD
# Class
biom convert -i data/Nestedness/biom/TBC/class.tsv \
-o data/Nestedness/biom/TBC/class_hdf5.biom \
--table-type="OTU table" --to-hdf5
biom add-metadata -i data/Nestedness/biom/TBC/class_hdf5.biom \
-o data/Nestedness/biom/TBC/class_hdf5_meta.biom \
--sample-metadata-fp data/Nestedness/biom/TBC/metadata.tsv \
--sample-header donor_id,donor_group,PD
# Phylum
biom convert -i data/Nestedness/biom/TBC/phylum.tsv \
-o data/Nestedness/biom/TBC/phylum_hdf5.biom \
--table-type="OTU table" --to-hdf5
biom add-metadata -i data/Nestedness/biom/TBC/phylum_hdf5.biom \
-o data/Nestedness/biom/TBC/phylum_hdf5_meta.biom \
--sample-metadata-fp data/Nestedness/biom/TBC/metadata.tsv \
--sample-header donor_id,donor_group,PD

#--------------------------------------------------------

### RUSH Data
# Species
biom convert -i data/Nestedness/biom/RUSH/species.tsv \
-o data/Nestedness/biom/RUSH/species_hdf5.biom \
--table-type="OTU table" --to-hdf5
biom add-metadata -i data/Nestedness/biom/RUSH/species_hdf5.biom \
-o data/Nestedness/biom/RUSH/species_hdf5_meta.biom \
--sample-metadata-fp data/Nestedness/biom/RUSH/metadata.tsv \
--sample-header donor_id,donor_group,PD
# Genera
biom convert -i data/Nestedness/biom/RUSH/genera.tsv \
-o data/Nestedness/biom/RUSH/genera_hdf5.biom \
--table-type="OTU table" --to-hdf5
biom add-metadata -i data/Nestedness/biom/RUSH/genera_hdf5.biom \
-o data/Nestedness/biom/RUSH/genera_hdf5_meta.biom \
--sample-metadata-fp data/Nestedness/biom/RUSH/metadata.tsv \
--sample-header donor_id,donor_group,PD
# Family
biom convert -i data/Nestedness/biom/RUSH/family.tsv \
-o data/Nestedness/biom/RUSH/family_hdf5.biom \
--table-type="OTU table" --to-hdf5
biom add-metadata -i data/Nestedness/biom/RUSH/family_hdf5.biom \
-o data/Nestedness/biom/RUSH/family_hdf5_meta.biom \
--sample-metadata-fp data/Nestedness/biom/RUSH/metadata.tsv \
--sample-header donor_id,donor_group,PD
# Order
biom convert -i data/Nestedness/biom/RUSH/order.tsv \
-o data/Nestedness/biom/RUSH/order_hdf5.biom \
--table-type="OTU table" --to-hdf5
biom add-metadata -i data/Nestedness/biom/RUSH/order_hdf5.biom \
-o data/Nestedness/biom/RUSH/order_hdf5_meta.biom \
--sample-metadata-fp data/Nestedness/biom/RUSH/metadata.tsv \
--sample-header donor_id,donor_group,PD
# Class
biom convert -i data/Nestedness/biom/RUSH/class.tsv \
-o data/Nestedness/biom/RUSH/class_hdf5.biom \
--table-type="OTU table" --to-hdf5
biom add-metadata -i data/Nestedness/biom/RUSH/class_hdf5.biom \
-o data/Nestedness/biom/RUSH/class_hdf5_meta.biom \
--sample-metadata-fp data/Nestedness/biom/RUSH/metadata.tsv \
--sample-header donor_id,donor_group,PD
# Phylum
biom convert -i data/Nestedness/biom/RUSH/phylum.tsv \
-o data/Nestedness/biom/RUSH/phylum_hdf5.biom \
--table-type="OTU table" --to-hdf5
biom add-metadata -i data/Nestedness/biom/RUSH/phylum_hdf5.biom \
-o data/Nestedness/biom/RUSH/phylum_hdf5_meta.biom \
--sample-metadata-fp data/Nestedness/biom/RUSH/metadata.tsv \
--sample-header donor_id,donor_group,PD

