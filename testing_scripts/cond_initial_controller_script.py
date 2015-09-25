import os
import sys
current_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(current_path + "/src/")
from cond_initial_controller import *
current_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
os.chdir(current_path + "/Nahoon_kida.uva.2014/")

def is_in_KIDA(old_vals, species):
    """Given a species, checks if the species is one of the species
    in the input parameter file."""
    for val in old_vals:
        if val.species == species:
            return True
    return False

###########GETTING THE INITIAL CONDITION FILE#######################
f = 'cond_initial_kida.uva.2014.dat'
good_file = True

#good_file = False
#while good_file == False:
#    try:
#        f = raw_input("Specify the initial condition file (CTRL+C to quit): ")
#        good_file = True
#    except IOError as e:
#        print "I/O error({0}): {1}".format(e.errno, e.strerror)
#        good_file = False

################CHECKING IF THE FILE MUST BE CHANGED##########

if good_file == True:
    old_vals = read_old_vals(f)
    check_vals(old_vals)
    good_yes_no = False
    while good_yes_no == False:
        to_continue =  raw_input("Are these the initial abundances that you want for Nahoon? (YES or NO): ")   
        if to_continue == "YES" or to_continue == "yes" or to_continue == "Y" or to_continue == "y":
            print "You're good to go!"
            good_yes_no = True
        elif to_continue == "NO" or to_continue == "no" or to_continue == "N" or to_continue == "n":
            good_to_go = False
            while good_to_go == False:
                good_input_1 = False
                while good_input_1 == False:
                    ##########CHECKING WHAT SPECIES MUST BE CHANGED##############
                    to_change = raw_input("What should be changed? (separate all species to be changed with ', '): ")
                    params = to_change.split(', ')
                    for to_be_changed in params:
                        if is_in_KIDA(old_vals, to_be_changed) == True:
                            good_input_1 = True
                            good_input_2 = False
                            while good_input_2 == False:        
                                try:  
                                    ##########GETTING A NEW VALUE FOR THE SPECIES TO BE CHANGED##############      
                                    new_val = raw_input("What is the new value for the " + to_be_changed + " initial abundance? : ")
                                    print
                                    change_val(to_be_changed, old_vals, float(new_val))
                                    good_input_2 = True
                                except ValueError:
                                    print "That wasn't a valid number!!!!"
                        else:
                            print "Oops! " + to_be_changed + " is not in the initial abundance file you specified!"
                            good_input_1 = False
                            break
                    write_new_vals(f, old_vals)
                print "Done! New initial abundances: "
                print
                check_vals(old_vals)
                ##########RE-CHECKING IF THE PARAMETERS ARE CORRECT##############
                are_you_good =  raw_input("Are these the initial abundances that you want for Nahoon? (YES or NO): ")
                if are_you_good == "YES" or are_you_good == "yes" or are_you_good == "Y" or are_you_good == "y":
                    print "You're good to go!"
                    good_to_go = True
                    good_yes_no = True
        else:
            print "Whoops! That was a YES or NO question. Try again and don't mess up."

os.chdir(current_path + "/testing_scripts")