class initial_abundance:
    def __init__(self, data):
        """Initializes an initial abundance object using all the data from a
        given line of the cond_initial.dat file and that line's number."""
        self.data = data
        self.first_char = 0
        for i in range(len(data)):
            if data[i] != '':
                self.first_char = i
                break
        self.number = int(data[self.first_char])
        self.species = data[self.first_char + 1]
        self.abundance = float(data[-1][0:14].replace('D','E'))
    
    def get_fortran_abundance(self):
        """Formats all the data for an initial abundance in the standard KIDA format.
        Returns a string containing all the formatted data for the reaction."""
        return str('%.8E' % (self.abundance)).replace('E', 'D') + "\n"
    
    def rewrite(self):
        """After one of the instance variables in the initial_abundance is changed,
        this method rewrites the object to reflect the changed variable."""
        self.data[self.first_char] = str(self.number)
        self.data[self.first_char + 1] = self.species
        self.data[-1] = self.get_fortran_abundance()

def read_old_vals(filename):
    """This function takes in a cond_initial.dat file and creates an array of
    initial_abundance objects.""" 
    f = open(filename)
    all_old_vals = []
    for i, line in enumerate(f):
        if i > 0:
            data = line.split(" ")
            all_old_vals.append(initial_abundance(data))
    f.close()
    return all_old_vals

def check_vals(all_data):
    """Given the output array from a read_old_vals call, this method
    prints all the relevant initial condition data from the file."""
    for val in all_data:
        if val.abundance != 0:
            print str(val.number) + "\t" + val.species + "\t" + val.get_fortran_abundance()[:-1]
    
def change_val(species_to_be_changed, old_vals, new_val):
    """Given a value, this method formats that value in such a way that
    it is processable by Nahoon."""
    for val in old_vals:
        if val.species == species_to_be_changed:
            val.abundance = new_val

def write_new_vals(filename, all_data):
    """Given a filename to write to and an array of values (in the same format
    as the matrix obtained with a call to read_old_vals), this method writes
    that array back to the file while preserving the original cond_initial.dat
    formatting."""
    JSPACE_val = 0
    g = open(filename, 'w')
    g.write("JSPACE = " + str(JSPACE_val) + "\n")
    for val in all_data:
        val.rewrite()
        data = val.data
        for i in range(len(data)):
            if i == len(data) - 1:
                g.write(data[i])
            else:
                g.write(data[i] + " ")
    g.close()