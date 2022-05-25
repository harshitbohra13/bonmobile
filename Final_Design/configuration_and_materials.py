#Structures data


"""""
This is the configuration, you can initialize parameters that will be iterated

"""""


#Class is basically a datatype like integer or string
class configure:
    
    """""
    Initialize parameters here: This will give default values to any test run
    """""
    Rotor_number = 12
    Rotor_number = 12
    fuselage_rotor_clearance = 0
    safety_factor = 1.5
    def __init__(self, name):
        self.name = name


"""""
How to use config class?
Simple. 
Add all your parameters in the class above then follow the code below:

To import config into your file: 

from configuration_and_materials import config

To initialize the class:
bonmobile = config()

To use the properties, which you can change in your code: 
This will print the number of rotors

print(bonmobile.Rotor_number)

Why do we use this?
You can pass class in functions like variables!!! Basically

def fuselage_mass(bonmobile):
    Rotormax = bonmobile.Rotor_number * 69 * 69 
    return(bonmobile)

This function basically sends all the parameters to your function and back, not great for computations but ezz for iterations
basically basically
"""""

class materials:
    """
    Material class
    Feel free to add more properties like I did for E modulus
    """
    def __init__(self, name, E):
        self.name = name
        self.E = E



"""""
This is an example of how to initialize a material. 
You can always add more material properties and they need not be initialized for all the materials

"""""
Aluminiumbadboy = materials("Aluminium", E = 201301) 