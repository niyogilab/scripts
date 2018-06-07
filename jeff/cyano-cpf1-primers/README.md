qrlabels
========

Generates sheets of QR codes that can be printed on labels. I use the
[example commands](https://github.com/jefdaj/qrlabels/tree/master/examples) with
["Tough Spot" labels](https://www.amazon.com/Diversified-Tough-Spots-SPOT-1000-Polyvinyl-Diameter/dp/B00APK3EJE)
to track laboratory samples, but the program should work for anything that
comes as a grid of stickers on regular-size paper. Be sure to use the highest
quality print settings for readability of small QR codes.

```
Usage:
  qrlabels [-v] -p PREFIX -n NCHAR -m MARGINS -d DIMENSIONS -s SIDE PDFPATH
  qrlabels --help

Options:
  --help         Show this text
  -v             Print the text of each QR code [defualt: False]
  -p PREFIX      Text to start QR codes with, for example "http://my.site/qrcode/"
  -n NCHAR       Number of characters in the random portion of each QR code
  -d DIMENSIONS  Dimensions of the QR codes table in rows and columns,
                 for example "4x6" or "12x16". No spaces!
  -s SIDE        Length of each side of the QR codes. Allowed units: mm, cm, in, px
  -m MARGINS     Comma-separated margins. Specify one, two, or all four:
                 * one  means all margins are equal
                 * two  means top/bottom, left/right
                 * four means top, bottom, left, right
                 Allowed units are mm, cm, in, px. No spaces!
                 Examples: "10mm,20mm" "1in" "25px,25px,20px,20px"
     PDFPATH     Where to save the output PDF
```

Features I plan to implement eventually (but would love help with!):

* "Test mode" that prints rulers or something on the paper
* Accept a list of codes from a file rather than generating them
* Groups of repeated barcodes for when you want more than one of each
* Keep a history file of previous codes to prevent short "UUIDs" from colliding
* Print multiple pages at once
* Add optional text labels below/around the QR codes
* Support paper sizes besides letter
