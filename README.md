KNLab Scripts
=============

Miscellaneous scripts used in the [Kris Niyogi lab][1] at UC Berkeley.
Some might be well-documented, some might be useful, but nothing is guaranteed!
Each lab member has their own setup, so check the README in their subfolder for instructions.

| script | description |
|:-------|:------------|
| simple-growth-curves | Plots growth curves given a spreadsheet of OD measurements and a list of treatments
| led-array-programs | Describe a pattern of light levels in Python, and this writes a program for our custom growth chamber LED arrays
| [merge-metadata](jeff/qrlabels/) | Merges multiple plate maps into one big tidy metadata table
| tidy-tecan | Convert Excel spreadsheets generated by Tecan plate readers into tidy data for plotting in R
| [merge-tecan](jeff/mergetecan/) | Merges the outputs of `merge-metadata` and `tidy-tecan`
| tecan-growth-curves | High-throughput growth curves using the output of `merge-tecan`
| tecan-absorbance-scans | Plots absorbance scans from Tecan data
| [pcc7942-cpf1-primers](jeff/pcc7942-cpf1-primers/) | Generates initial primers to order for editing the PCC 7942 genome with CRISPR/Cpf1
| [qrcode-labels](jeff/qrlabels/) | Generates sheets of random QR codes to print on labels

[1]: http://niyogilab.berkeley.edu
