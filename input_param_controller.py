import re

def read_old_vals(filename):
    """Reads all the data from the input parameter file into an array, and
    returns the array. The array also contains all comments - this is done to
    preserve formatting."""
    vals = []
    comments = []
    f = open(filename, 'r')
    for i, line in enumerate(f):
	if i > 8 and i < 18: #these are the line numbers not in the initial comment
            data = re.split("\t\t\t?", line)
            vals.append(data[0])
            comments.append(data[1])
    f.close()
    return [vals, comments]

def check_vals(all_data):
    """Given the output array from a read_old_vals call, this method
    prints all the relevant input parameter data from the file."""
    print "H density: " +  str(all_data[0][1])
    print "Temperature: " +  str(all_data[0][2])
    print "Visual Extinction: " + str(all_data[0][3])
    print "Cosmic Ray Ionization Rate: " + str(all_data[0][5])

def fortran_format(new_val):
    """Given a value, this method formats that value in such a way that
    it is processable by Nahoon."""
    return str('%.3E' % (new_val)).replace('E', 'D')

def change_val(to_be_changed, old_vals, new_val):
    """Given a parameter to be changed, this method goes into the array of
    old values (obtained from a call to read_old_vals) and replaces the
    appropriate value with new_val."""
    if to_be_changed == 'TEMPERATURE':
        old_vals[0][2] = fortran_format(float(new_val))
    elif to_be_changed == 'DENSITY':
        old_vals[0][1] = fortran_format(float(new_val))
    elif to_be_changed == 'EXTINCTION':
        old_vals[0][3] = fortran_format(float(new_val))
    elif to_be_changed == 'CRIR':
        old_vals[0][5] = fortran_format(float(new_val))
    else:
        return "ERROR: You didn't specify a valid parameter to be changed!"

def write_new_vals(filename, params):
    """Given a filename to write to and an array of values (in the same format
    as the matrix obtained with a call to read_old_vals), this method writes
    that array back to the file while preserving the original input_parameter.dat
    formatting."""
    initcomment = "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc\nc	input_parameter.dat\nc	May 2006 Wakelam Valentine\nc	input file of the parameters\nc	Please respect the format\nc	If you are running the uncertainties, the gas temperature and density \nc	will not be read from this file\ncccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc\ncccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc\n"
    vals = params[0]
    comments = params[1]
    g = open(filename, 'w')
    g.write(initcomment)
    for i in range(len(vals)):
        spacing = "\t\t"
        if i == 0:
            spacing += "\t"
        g.write(str(vals[i]) + spacing + str(comments[i]))
    g.close()