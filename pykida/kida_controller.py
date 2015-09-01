class reaction:
    def __init__(self, data, index):
        """Initializes a reaction using all the data from a given line of the
        kida.uva.xxxx file and that line's number."""
        self.data = data
        components = self.data.split()
        self.index = index
        self.reactants = components[0:2]
        self.products = components[2:len(components) - 13]
        self.rate_constant_coefficients = components[-13:-10]
        self.uncertainty_coefficients = components[-10:-8]
        self.uncertainty_type = components[-8]
        self.type_of_reaction = components[-7]
        self.valid_temperature_range = components[-6:-4]
        self.formula = components[-4]
        self.number = components[-3]
        self.number_of_same_reaction = components[-2]
        self.recommendation = components[-1]
    
    def to_kida_format(self):
        """Formats all the data for a reaction in the standard KIDA format.
        Returns a string containing all the formatted data for the reaction."""
        kida = self.data.rsplit("               ", 1)
        
        coeffs = self.rate_constant_coefficients[0]
        if float(self.rate_constant_coefficients[1]) < 0:
            coeffs += (" " + self.rate_constant_coefficients[1])
        else:
            coeffs += ("  " + self.rate_constant_coefficients[1])
            
        if float(self.rate_constant_coefficients[2]) < 0:
            coeffs += (" " + self.rate_constant_coefficients[2])
        else:
            coeffs += ("  " + self.rate_constant_coefficients[2])
        
        form = str(self.formula)
        if(len(form) == 1):
            string = kida[0] + "               " + coeffs + kida[1][len(coeffs):73] + form + kida[1][74:]
        elif(len(form) == 2):
            string = kida[0] + "               " + coeffs + kida[1][len(coeffs):72] + form + kida[1][74:]
        return string
            
    def visualize(self):
        """Provides a quick way to visualize a reaction as:
        reactant1 + reactant2 ----> product1 + product2 + ..."""
        visual = ""
        for i, reactant in enumerate(self.reactants):
            if i == len(self.reactants) - 1:
                visual += reactant + " "
            else:
                visual += (reactant + " + ")
        visual += "----> "
        for i, product in enumerate(self.products):
            if i == len(self.products) - 1:
                visual += product
            else:
                visual += (product + " + ")
        print visual

def initialize_network(filename):
    """This function takes in a kida.uva.xxxx file and creates an array of
    reaction objects. This array is quickly iterable, allowing for faster 
    traversal of the network than currently exists on the KIDA webpage.""" 
    network = []
    network_file = open(filename)
    for i, line in enumerate(network_file):
        network.append(reaction(line, i))
    return network

def is_in_KIDA(reactants, products, network):
    """This function takes in an array or tuple of reactants and products and
    a network (an array of reaction objects). It checks to see if a reaction
    with those reactants and products exists within the network. It returns
    (True, reaction) if the reaction is found in the network, and (False, -1)
    if not."""
    for reaction in network:
        if sorted(reaction.reactants) == sorted(reactants):
            if sorted(reaction.products) == sorted(products):
                return (True, network[reaction.index])
    return (False, -1)

def fortran_format(new_vals):
    """This function takes an array of values and returns those values
    formatted for use in the kida.uva.xxxx file (such that Nahoon
    doesn't complain about the formatting)."""
    return_vals = []
    for val in new_vals:
        return_vals.append(str('%.3E' % (float(val))).replace('E', 'e'))
    return return_vals

def change_val(reaction, network, new_rate_coef_constants):
    """This function takes a reaction, the network in which that reaction
    appears, and an array of rate coefficients constants. It changes the
    reaction in the network to use the new rate coefficient constants."""
    network[reaction.index].rate_constant_coefficients = fortran_format(new_rate_coef_constants)
    
def change_form(reaction, network, new_formula):
    """This function takes a reaction, the network in which that reaction
    appears, and a new formula number to calculate the rate coefficient.
    It changes the reaction in the network to use the new formula number."""
    network[reaction.index].formula = str(new_formula)

def write_new_vals(network, filename):
    """This takes a network of reactions and a filename to write those reactions
    to. It formats those reactions in a KIDA-appropriate fashion and writes
    them to the file specified."""
    g = open(filename, 'w')
    for val in network:
        g.write(val.to_kida_format())
    g.close()