import os
current_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
os.chdir(current_path + "/pykida")
from kida_controller import *
current_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
os.chdir(current_path + "/Nahoon_kida.uva.2014")

#f = 'kida.uva.2014'
####################GETTING THE KIDA FILE###########
good_file = False
while good_file == False:
    try:
        f = raw_input("Specify the initial condition file (CTRL+C to quit): ")
        old_vals = initialize_network(f)
        good_file = True
    except IOError as e:
        print "I/O error({0}): {1}".format(e.errno, e.strerror)
        good_file = False

if good_file == True: 
    good_yes_no = False
    while good_yes_no == False:
        ################CHECKING IF THE FILE MUST BE CHANGED#####################
        to_continue =  raw_input("Do you want to leave the " + f + " network unchanged? (YES or NO): ")   
        if to_continue == "YES" or to_continue == "yes" or to_continue == "Y" or to_continue == "y":
            print "You're good to go!"
            good_yes_no = True
        elif to_continue == "NO" or to_continue == "no" or to_continue == "N" or to_continue == "n":
            good_to_go = False
            while good_to_go == False:
                how_many_input = False
                while how_many_input == False:
                    try:
                        ##########CHECKING HOW MANY REACTIONS TO CHANGE##############
                        how_many = raw_input("How many reactions would you like to change?: ")
                        how_many_input = True
                        for i in range(int(how_many)):
                            good_input_1 = False
                            while good_input_1 == False:
                                ##############GETTING THE REACTIONS TO CHANGE#################
                                to_change = raw_input("What reaction should be changed? (enter reactants followed by a '; ' followed by the products, and separate the species within the reactants and products with ', '): ")
                                reactants_and_products = to_change.split('; ')
                                kida_tuple = is_in_KIDA(reactants_and_products[0].split(', '), reactants_and_products[1].split(', '), old_vals)
                                if kida_tuple[0] == True:
                                    reaction = kida_tuple[1]
                                    good_input_1 = True
                                    good_input_2 = False
                                    while good_input_2 == False:        
                                        try:
                                            ##########GETTING THE NEW COEFFICIENTS FOR THE REACTION################        
                                            new_vals_raw = raw_input("What are the new Arrhenius coefficients for the " + reaction.visualize() + " reaction? (separate coefficients with ', '): ")
                                            print
                                            new_vals = new_vals_raw.split(", ")
                                            fortran_formatted_vals = fortran_format(new_vals)
                                            change_val(reaction, old_vals, fortran_formatted_vals)
                                            good_input_2 = True
                                        except ValueError:
                                            print "That wasn't a valid number!!!!"
                                else:
                                    print "Oops! That reaction is not in KIDA!"
                                    good_input_1 = False
                                    
                    except ValueError:
                        print "That wasn't a valid number!!!!"
                        how_many_input = False
                        
                write_new_vals(f, old_vals)
                print "Done! The changes have been saved in " + f
                ##########ASKING FOR ANY ADDITIONAL CHANGES##############
                are_you_good =  raw_input("Want to make any more changes to the " + f + " network? (YES or NO): ")
                if are_you_good == "NO" or are_you_good == "no" or are_you_good == "N" or are_you_good == "n":
                    print
                    print "You're good to go!"
                    good_to_go = True
                    good_yes_no = True
        else:
            print "Whoops! That was a YES or NO question. Try again and don't mess up."