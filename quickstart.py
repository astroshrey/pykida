import os
import sys
import shutil
import numpy as np
current_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(current_path + "/src/")
import input_param_controller
import kida_controller
os.chdir(current_path + "/testing_scripts/")

def fitfunc(T, A, N, c_1, c_2, c_3, c_4, T_1, T_2, T_3, T_4):
    """Here, we define a non-standard fitting function. Currently, these functions can't be accepted into KIDA.
        This script allows the inclusion of a non-standard function."""
    return A*((300/T)**N) + T**(-3.0/2)*c_1*np.exp(-T_1/T) + T**(-3.0/2)*c_2*np.exp(-T_2/T) + T**(-3.0/2)*c_3*np.exp(-T_3/T) + T**(-3.0/2)*c_4*np.exp(-T_4/T)

#setting constants and array of temps
command1 = 'sh build.sh'
command2 = 'sh run.sh'
input_file = 'input_parameter.dat'
kida_file = 'kida.uva.2014'
output_directory = "output"
rate_coefficient_parameters = [7.39E-10, 1.46E-1, 6.32E-8, -3.00E-6, -1.17E-5, 5.76E-4, 7.47E1, 5.40E2, 1.83E3, 1.90E4]
temps = np.logspace(np.log10(10), np.log10(400), 100)

#executing some initial checks
execfile('cond_initial_controller_script.py')
execfile('header_controller_script.py')

os.chdir(current_path + "/Nahoon_kida.uva.2014/")

#creating the directory into which plot.dat files will be placed
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

for i, temp in enumerate(temps):
    #changing the temperature in input_parameter.dat
    old_vals = input_param_controller.read_old_vals(input_file)
    input_param_controller.change_val('TEMPERATURE', old_vals, temp)
    input_param_controller.write_new_vals(input_file, old_vals)
    input_param_controller.check_vals(old_vals)
    
    #changing the rate coefficient value in KIDA to the correct value (based on the previously defined "fitfunc")
    rate_coefficient = fitfunc(temp, *rate_coefficient_parameters)
    old = kida_controller.initialize_network(kida_file)
    reaction = kida_controller.is_in_KIDA(['O', 'H3+'], ['H', 'H2O+'], old)[1]
    kida_controller.change_val(reaction, old, [rate_coefficient,0,0])
    kida_controller.write_new_vals(old, kida_file)
    
    #running KIDA
    os.system(command1)
    os.system(command2)
    
    #copying and shutiling the plot.dat and verif.dat files to the new directory
    f = open('plot.dat')
    g = open(str(temp), 'w')
    shutil.copyfileobj(f, g)
    f.close()
    g.close()
    shutil.move(str(temp), './' + output_directory)

shutil.move(output_directory, current_path)