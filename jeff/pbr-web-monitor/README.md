ePBR Web Monitor
================

We have some photobioreactors, which are really cool and run by a program called AlgalCommand.
Unfortunately they're also 1) fairly unreliable and 2) way down in the basement.
I put together this remote monitoring script so we can intervene quickly when something goes wrong.
For example if a temperature sensor breaks and starts reading 1000 degrees, the heaters will shut off to compensate and we need to fix that before the cells get too cold.

AlgalCommand continuously updates a bunch of separate data tables.
This merges them together, plots selected variables, and `rsync`s the plots to my bench computer.
It also checks that the rolling average of recent measurements is within a defined range and sends an email alert if not.
