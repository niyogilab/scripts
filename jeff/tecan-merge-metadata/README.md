Merge Tecan readings with Metadata
==================================

This is the final step in my standard Tecan workflow: after making one big
metadata table with [tecan-merge-metadata](../tecan-merge-metadata/), and one
or more tables of Tecan readings with
[tecan-extract-table](../tecan-extract-table/), use this to merge them. You get
Tecan readings annotated with well-level and plate-level metadata so they're
easy to plot in R (or Python, probably).

Note that this drops readings which don't match a well with metadata! If you
want to make sure the right ones are being kept/dropped, turn up the verbosity
with `-v`.
