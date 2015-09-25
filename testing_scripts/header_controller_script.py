import os
import sys
current_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(current_path + "/src/")
from header_controller import *
current_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
os.chdir(current_path + "/Nahoon_kida.uva.2014")

def is_valid_param(param):
    """Given a parameter string, this method checks whether that parameter
    is a modifiable parameter in header.f90."""
    return param == "NS" or param == "NRTOT" or param == "NX" or param == "NTIME" or param == "NELEM"

##################GETTING THE HEADER.F90 FILE#####################
f = 'header.f90'
good_file = True

#good_file = False
#while good_file == False:
#    try:
#        f = raw_input("Specify the header file (CTRL+C to quit): ")
#        good_file = True
#    except IOError as e:
#        print "I/O error({0}): {1}".format(e.errno, e.strerror)
#        good_file = False

##################CHECKING IF THE INFORMATION IN THE HEADER.F90 FILE IS CORRECT#####################
if good_file == True:
    old_vals = read_old_vals(f)
    check_vals(old_vals)
    good_yes_no = False
    while good_yes_no == False:
        to_continue =  raw_input("Are these the header parameters that you want for Nahoon? (YES or NO): ")
        if to_continue == "YES" or to_continue == "yes" or to_continue == "Y" or to_continue == "y":
            print "You're good to go!"
            good_yes_no = True
        elif to_continue == "NO" or to_continue == "no" or to_continue == "N" or to_continue == "n":
            good_to_go = False
            while good_to_go == False:
                good_input_1 = False
                ##################GETTING THE PARAMETERS TO BE CHANGED#####################
                while good_input_1 == False:
                    print
                    print "To change the number of species, type NS"
                    print "To change the number of reactions, type NRTOT"
                    print "To change the number of times writted to plot.dat, type NTIME"
                    print "To change the number of spatial points, type NX"
                    print "To change the number of elements (+1 for charges), type NELEM"
                    print "NOTE: if you want to change the shielding parameters, you can't do it through this program. You'll have to do it yourself!"
                    to_change = raw_input("What should be changed? (separate all parameters to be changed with ', '): ")
                    params = to_change.split(', ')
                    for to_be_changed in params:
                        if is_valid_param(to_be_changed) == True:
                            good_input_1 = True
                            good_input_2 = False
                            while good_input_2 == False:        
                                try:
                                    ##################GETTING NEW VALS FOR THE PARAMETER TO BE CHANGED#####################
                                    new_val = raw_input("What is the new value for " + to_be_changed + "? : ")
                                    change_val(to_be_changed, old_vals, int(new_val))
                                    good_input_2 = True
                                except ValueError:
                                    print "That wasn't a valid number!!!!"
                        else:
                            print "Oops! " + to_be_changed + " isn't a valid header parameter to be changed!"
                            good_input_1 = False
                            break
                    write_new_vals(f, old_vals)
                print "Done! New parameters: "
                print
                check_vals(old_vals)
                ##################CHECKING IF THE UPDATED PARAMETERS ARE CORRECT#####################
                are_you_good =  raw_input("Are these the header parameters that you want for Nahoon? (YES or NO): ")
                if are_you_good == "YES" or are_you_good == "yes" or are_you_good == "Y" or are_you_good == "y":
                    print "You're good to go!"
                    good_to_go = True
                    good_yes_no = True
        else:
            print "Whoops! That was a YES or NO question. Try again and don't mess up."
os.chdir(current_path + "/testing_scripts")