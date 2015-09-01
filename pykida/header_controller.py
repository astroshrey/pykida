import re

def read_old_vals(filename):
    """Reads all the values from a specified header.f90 file, and returns a
    dictionary of each value with its key ('NS' -> 489, 'NRTOT -> 6900', etc.)"""
    f = open(filename)
    all_old_vals = {}
    for i, line in enumerate(f):
        if i == 12:
            j = 0
            equals =  [(m.start(0), m.end(0)) for m in re.finditer("=", line)]
            commas =  [(m.start(0), m.end(0)) for m in re.finditer("[,)]", line)]
            while j < 5:
                if j == 0:
                    all_old_vals["NS"] = int(line[equals[j][1]:commas[j][0]])
                if j == 1:
                    all_old_vals["NRTOT"] = int(line[equals[j][1]:commas[j][0]])
                if j == 2:
                    all_old_vals["NTIME"] = int(line[equals[j][1]:commas[j][0]])
                if j == 3:
                    all_old_vals["NX"] = int(line[equals[j][1]:commas[j][0]])
                if j == 4:
                    all_old_vals["NELEM"] = int(line[equals[j][1]:commas[j][0]])
                j += 1
    f.close()
    return all_old_vals

def check_vals(all_data):
    """Given a dictionary of all the values and keys in a header.f90 file,
    this method prints out a summary of those values and keys."""
    print "Number of Species: " + str(all_data['NS'])
    print "Number of Reactions in KIDA: " + str(all_data['NRTOT'])
    print "Number of Output Times in plot.dat: " + str(all_data['NTIME'])
    print "Number of Spatial Points: " + str(all_data['NX'])
    print "Number of Elements + 1 (for charges): " + str(all_data['NELEM'])
    
def change_val(param_to_change, old_vals, new_val):
    """Given a parameter to change and a dictionary of values and keys from
    a header.f90 file, this method replaces the appropriate parameter's value
    with new_val."""
    old_vals[param_to_change] = new_val

def rewrite(all_data):
    """Given a dictionary of values and keys, this method creates
    the appropriate new line for replacement in the header.f90 file."""
    tags = ['NS', 'NRTOT', 'NTIME', 'NX', 'NELEM']
    new_line = "       parameter ("
    new_line += (tags[0] + "=" + str(all_data[tags[0]]) + ", ")
    new_line += (tags[1] + "=" + str(all_data[tags[1]]) + ", ")
    new_line += (tags[2] + "=" + str(all_data[tags[2]]) + ",")
    new_line += (tags[3] + "=" + str(all_data[tags[3]]) + ",")
    new_line += (tags[4] + "=" + str(all_data[tags[4]]) + ")\t\n")
    return new_line

def write_new_vals(filename, all_data):
    """Given a dictionary of values and keys, this method writes a
    header.f90 file while preserving the formatting necessary for
    Nahoon to process the data."""
    initial_comment = "!       NS = NUMBER OF SPECIES\n!       NRTOT = NUMBER OF REACTIONS IN THE NETWORK (THE REAL ONE!!!)\n!       NX = NUMBER OF SPATIAL POINTS\n!       NTIME = NUMBER OF OUTPUT TIMES\n!       NELEM = NUMBER OF ELEMENTS + 1 (FOR CHARGES)\n!       XCO_0 and XH2_0 are the abundances of CO and H2 at the edge of \n!       the cloud to compute shelf sheilding of CO and H2. ONE NEEDS TO PAY\n!       ATTENTION TO THESE PARAMETERS SINCE THEY CAN STRONGLY AFFECT THE \n!       CHEMISTRY\n!       \n"
    line1 = "       integer NS, NRTOT, NX, NTIME, NELEM \n"
    line2 = "       double precision XCO_0,XH2_0\n"
    line3 = rewrite(all_data)
    line4 = "       parameter (XCO_0=1.D-4, XH2_0=5.D-1)\n"
    rewrite(all_data)
    g = open(filename, 'w')
    g.write(initial_comment)
    g.write(line1)
    g.write(line2)
    g.write(line3)
    g.write(line4)
    g.close()