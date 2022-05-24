from sqlalchemy import null
from Battery_mass import Max_Total_thrust
from configuration_and_materials import Rotor_number

# Total_Supporting_structure_mass = 0


#Thrust per rotor
Max_thrust_per_rotor = Max_Total_thrust/Rotor_number

def sup_struct_mass(Rotor_number, ):
    return null