Tecan Growth Curves
===================

Overall protocol
----------------

To prepare data from a Tecan timecourse, you should:

#. Write a table of plate metadata describing the conditions each plate is
   grown under, which antibiotics they use etc. I usually name this
   `<start date>_<experiment name>_plates.csv`.

#. For each field that doesn't fit in one cell because it varies across
   a plate, write a separate plate map and reference it from the main table.
   These files need to have names that include the *column* they appear in
   (`antibiotic`, `light`, `strain` etc) and also the identifier used in
   their cell (`p0210` if you referred to `@p0210`). I usually name them
   `<start date>_<experiment name>_<plate id>_<column name>.csv`. It sounds
   brittle, but works pretty well in practice. `mergemeta` will abort if
   a reference is ambiguous (but the error might be less than helpful).

#. Run `mergemeta` to consolidate the above files into a master metadata
   table. You only need to do this once at the beginning of the experiment,
   unless you decide to change conditions later.

#. As the timecourse progresses and you accumulate `.xlsx` files, run
   `tidytecan` on each of them. If you're measuring multiple things (I often
   do OD~750~ to track growth and also a full absorbance scan from 400 to
   800nm to track pigments), you'll want to set their labels in iControl
   or Magellan so you can extract them separately with the `--label` option.

#. At any point during or after the timecourse, you can use `mergetecan` to
   make a big table of annotated Tecan readings from your metadata table +
   readings. This is what you'll want when graphing with `ggplot` or doing
   other analyses.
