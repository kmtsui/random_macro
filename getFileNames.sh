gfal-ls root://hepgrid11.ph.liv.ac.uk//t2k.org/t2kliverpooldisk/production6T/mc/neut/flattree/run8a/root >filenames.txt
awk '{ printf "root://hepgrid11.ph.liv.ac.uk//t2k.org/t2kliverpooldisk/production6T/mc/neut/flattree/run8a/root/"; print }' filenames.txt > fileNames_prod6T_run8a.txt
rm filenames.txt
split -l 20 -d fileNames_prod6T_run8a.txt fileNames_prod6T_run8a.txt_

gfal-ls root://hepgrid11.ph.liv.ac.uk//t2k.org/t2kliverpooldisk/production6T/mc/neut/flattree/run8w/root >filenames.txt
awk '{ printf "root://hepgrid11.ph.liv.ac.uk//t2k.org/t2kliverpooldisk/production6T/mc/neut/flattree/run8w/root/"; print }' filenames.txt > fileNames_prod6T_run8w.txt
rm filenames.txt
split -l 20 -d fileNames_prod6T_run8w.txt fileNames_prod6T_run8w.txt_

gfal-ls root://hepgrid11.ph.liv.ac.uk//t2k.org/t2kliverpooldisk/production6T/mc/neut/flattree/run7/root >filenames.txt
awk '{ printf "root://hepgrid11.ph.liv.ac.uk//t2k.org/t2kliverpooldisk/production6T/mc/neut/flattree/run7/root/"; print }' filenames.txt > fileNames_prod6T_run7.txt
rm filenames.txt
split -l 20 -d fileNames_prod6T_run7.txt fileNames_prod6T_run7.txt_
