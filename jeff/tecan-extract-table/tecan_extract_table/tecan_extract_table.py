from csv        import reader, writer
from docopt     import docopt
from glob       import glob
from logging    import Logger, StreamHandler
from os.path    import join, splitext, basename
from pkgutil    import get_data
from shutil     import rmtree
from subprocess import check_call
from sys        import stdout, stderr
from tempfile   import mkdtemp
from time       import strptime, strftime

### logging ###

log = Logger('tecan-extract-table')
log.addHandler(StreamHandler(stderr))

### parse each type of tecan reading ###

def parse_readings_single(rows):
    cols = []
    readings = []
    for row in rows:
        if len(cols) > 0:
            if len(row[0]) == 0 or row[0] == 'End Time:':
                # reached end of table
                return readings
            else:
                for n in range(len(cols)):
                    well = row[0] + cols[n]
                    readings.append({'well': well, 'reading': row[n+1]})
        elif row[0] == '<>':
            # reached table header
            cols = [c for c in row[1:] if len(c) > 0]

# TODO: unify with parse_readings_scan?
def parse_readings_multi(rows):
    header   = []
    readings = []
    for row in rows:
        row = [c for c in row if len(c) > 0]
        if header:
            if len(row[0]) == 0 or row[0] == 'End Time:':
                # reached end of table
                return readings
            for n in range(len(header)):
                (xpos, ypos) = header[n].split(';')
                reading = \
                    { 'well'    : row[0]
                    , 'xpos'    : xpos
                    , 'ypos'    : ypos
                    , 'reading' : row[n+3]
                    }
                readings.append(reading)
        elif row[0] == 'Well':
            # reached table header
            header = row[3:]

def parse_readings_scan(rows):
    wlens    = []
    readings = []
    for row in rows:
        if wlens:
            if len(row) == 0 or row[0] == 'End Time:':
                # reached end of table
                return readings
            # one row represents the complete scan of one well
            for n in range(len(wlens)):
                reading = \
                    { 'well'       : row[0]
                    , 'wavelength' : wlens[n]
                    , 'reading'    : row[n+1]
                    }
                readings.append(reading)
        elif row[0] == 'Wavel.':
            # reached table header
            wlens = [c for c in row[1:] if len(c) > 0]

def corner_matches(rows, label):
    for row in rows:
        if row[0] == label:
            return True
    return False

def parse_readings(rows):
    'dispatches to one of the section parsers above'
    parser = None
    parsers = \
        { parse_readings_multi  : 'Well'
        , parse_readings_scan   : 'Wavel.'
        , parse_readings_single : '<>'
        }
    for fn in parsers:
        corner = parsers[fn]
        if corner_matches(rows, corner):
            parser = fn
            break
    if parser is None:
        raise Exception('No parser written for this type of data yet :(')
    log.info('    dispatching to %s' % parser.__name__)
    return parser(rows)

### parse whole tecan files ###

def xlsx2csv(xlsxfile, tmpd):
    cmd = ['xlsx2csv', '--all', xlsxfile, tmpd]
    check_call(cmd)
    csv = glob(join(tmpd, '*.csv'))
    return csv

def parse_csv(infile):
    with open(infile, 'rb') as f:
        try:
            rdr = reader(f)
            return list(r for r in rdr)
        except:
            msg = "Could not load csv file '%s'" % infile
            raise Exception(msg)

def parse_time(cells):
    # time rows look like: Start Time:,5/3/2015 6:52:20 PM,,,,
    for line in cells:
        if line[0] == 'Start Time:':
            time = strptime(line[1], '%m/%d/%Y %I:%M:%S %p')
            time = strftime('%Y-%m-%d %H:%M:%S', time)
            return time
    print(cells)
    raise Exception('failed to parse time')

def parse_sheet(rows, label):
    'extracts sections (lists of rows) with the given label'
    sections = []
    current  = []
    for row in rows:
        if row[0].startswith('Label:') and len(current) > 0:
            sections.append(current)
            current = []
        current.append(row)
    if len(current) > 0:
        sections.append(current)
    sections = [s for s in sections if s[0][0].find(label) != -1]
    return sections

def parse_tecan(tecanpath, csvpath, label):
    meta = {'plate': splitext(basename(csvpath))[0]}
    cells = parse_csv(csvpath)
    if cells == []:
        return []
    # if len(cells) == 0:
    #     msg = 'empty CSV file: %s' % csvpath
    #     raise Exception(msg)
    log.info('  parsing %s' % csvpath)
    sections = parse_sheet(cells, label)
    data = []
    # if len(sections) == 0:
    #     msg = "No sections matching '%s' found in %s"
    #     raise Exception(msg % (label, csvpath))
    for section in sections:
        log.debug('cells in section: %s' % section)
        data.append(parse_readings(section))
    meta['measure'] = parse_time(cells)
    rows = []
    for section in data:
        for reading in section:
            reading.update(meta)
            rows.append(reading)
    return rows

def parse_tecans(tecanpaths, label):
    allrows = []
    for tecanpath in tecanpaths:
        log.info('parsing %s readings from %s' % (label, tecanpath))
        try:
            tmpd = mkdtemp()
            csvpaths = xlsx2csv(tecanpath, tmpd)
            for csvpath in csvpaths:
                newrows = parse_tecan(tecanpath, csvpath, label)
                allrows += newrows
        except:
            raise
        finally:
            rmtree(tmpd)
    return allrows

### write tidy csv ###

def header(rows):
    klst = [r.keys() for r in rows]
    keys = set()
    for k in klst:
        new = set(k)
        keys = keys.union(new)
    return sorted(keys)

def write_rows(rows, outcsv):
    if outcsv:
        handle = open(outcsv, 'wb')
    else:
        handle = stdout
    try:
        wtr = writer(handle)
        hdr = header(rows)
        wtr.writerow(hdr)
        for point in rows:
            wtr.writerow([point[k] for k in hdr])
    except:
        raise
    finally:
        if outcsv:
            handle.close()

### command line interface ###

def main():
    args = docopt(get_data('tecan_extract_table', 'usage.txt'))
    log.setLevel(40 - 10 * args['-v'])
    log.info('read command line args:\n%s' % args)
    rows = parse_tecans(args['<tecan>'], args['--label'])
    write_rows(rows, args['--out'])
