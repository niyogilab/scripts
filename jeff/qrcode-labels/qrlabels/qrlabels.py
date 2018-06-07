#!/usr/bin/env python

# TODO figure out lingering margin issues:
#      - when size margins made really big, bottom increases too
#      - nrow = 1 or 20 expands to two pages

# TODO remove settings once printed on paper:
#      - any qr size .30-.35 works great with spot-2000 on white boxes
#      - 0.18 is best for spot-1000 on the sides of 200uL tubes (trickiest location)
#      - 0.22 is good for spot-1000 otherwise
#      - 0.45 margins for all

'''
Generates sheets of QR codes that can be printed on labels.

Usage:
  qrlabels [-v] -p PREFIX -n NCHAR -m MARGINS -d DIMENSIONS -s SIDE PDFPATH
  qrlabels (-h | --help)

Options:
  -h, --help     Show this text
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
'''

from docopt                        import docopt
from reportlab.graphics.barcode.qr import QrCodeWidget
from reportlab.graphics.shapes     import Drawing
from reportlab.lib.colors          import lightgrey
from reportlab.lib.pagesizes       import letter
from reportlab.lib.styles          import getSampleStyleSheet
from reportlab.lib.units           import mm, cm, inch
from reportlab.platypus            import SimpleDocTemplate, TableStyle, Paragraph
from reportlab.platypus.tables     import Table
from shortuuid                     import uuid
from sys                           import argv

### generate qrcodes ###

def qrcode(args):
  text = args['prefix'] + uuid()[:args['nchar']]
  if args['verbose']: print text
  return QrCodeWidget(text)

def qrcodes(args):
  rows = []
  for r in range(1, args['nrow']+1):
    row = []
    for c in range(1, args['ncol']+1):
      widget = qrcode(args)
      bounds = widget.getBounds()
      label  = Drawing(
        args['side'], args['side'],
        transform=[args['side']/bounds[2], 0, 0,
                   args['side']/bounds[3], 0, 0]
      )
      label.add(widget)
      row.append(label)
    rows.append(row)
  return rows

### generate pdf ###

def table(doc, style, args):
  data = qrcodes(args)
  widths  = [ doc.width                   / args['ncol']] * args['ncol']
  heights = [(doc.height - style.leading) / args['nrow']] * args['nrow']
  tbl = Table(
    data,
    colWidths  = widths,
    rowHeights = heights,
    hAlign     = 'CENTER',
    vAlign     = 'MIDDLE'
  )
  tbl.setStyle(TableStyle([
    ('ALIGN'    , (0,0), (-1,-1), 'CENTER'),
    ('VALIGN'   , (0,0), (-1,-1), 'MIDDLE'),
    ('INNERGRID', (0,0), (-1,-1), 1.0, lightgrey),
    ('BOX'      , (0,0), (-1,-1), 1.0, lightgrey)
  ]))
  return tbl

def kludge(canvas, doc):
  frame = doc.pageTemplates[0].frames[0]
  frame.leftPadding   = 0
  frame.rightPadding  = 0
  frame.topPadding    = 0
  frame.bottomPadding = 0
  canvas.saveState()

def pdf(args):
  style = getSampleStyleSheet()['Normal']
  cmd = Paragraph(' '.join(['qrlabels']+argv[1:]), style)
  doc = SimpleDocTemplate(
    args['pdfpath'],
    pagesize     = letter,
    rightMargin  = args['right'],
    leftMargin   = args['left'],
    topMargin    = args['top'] - style.leading,
    bottomMargin = args['bottom']
  )
  t = table(doc, style, args)
  return doc.build([cmd, t], onFirstPage=kludge, onLaterPages=kludge)

### parse arguments ###

def pixels(desc):
  lengths = {'mm': mm, 'cm': cm, 'in': inch, 'px': 1}
  for unit in lengths:
    index = desc.rfind(unit)
    if index != -1:
      return float(desc[:index]) * lengths[unit]
  raise Exception('Invalid measurement: %s' % desc)

def margins(desc):
  descs = desc.split(',')
  if len(descs) == 1:
    top = bottom = left = right = pixels(descs[0])
  elif len(descs) == 2:
    top  = bottom = pixels(descs[0])
    left = right  = pixels(descs[1])
  elif len(descs) == 4:
    top    = pixels(descs[0])
    bottom = pixels(descs[1])
    left   = pixels(descs[2])
    right  = pixels(descs[3])
  else:
    raise Exception('Invalid margins: %s' % desc)
  return (top, bottom, left, right)

def dimensions(desc):
  dims = desc.split('x')
  return int(dims[0]), int(dims[1])

def parse(args):
  top, bottom, left, right = margins(args['-m'])
  nrow, ncol = dimensions(args['-d'])
  return {
    'verbose' : args['-v'],
    'prefix'  : args['-p'],
    'top'     : top,
    'bottom'  : bottom,
    'left'    : left,
    'right'   : right,
    'side'    : pixels(args['-s']),
    'nchar'   : int(args['-n']),
    'nrow'    : nrow,
    'ncol'    : ncol,
    'pdfpath' : args['PDFPATH']
  }

def main():
  pdf(parse(docopt(__doc__, version='qrlabels 0.1')))
