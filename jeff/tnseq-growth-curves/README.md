TnSeq Growth Curves
===================

Takes a table of OD750 measurements and a table of treatment colors, and plots
growth curves. Takes dilutions into account to add up total generations grown,
which is useful for TnSeq experiments.

To run the example:

    nix-build
    ./result/bin/tnseq-growth-curves -c examples/colors.tsv -o examples/od750.tsv -p examples/plots2.png
