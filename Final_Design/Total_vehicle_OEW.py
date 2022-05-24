#importing different subsytems' mass
from Motor_mass import Total_Motor_mass
from Propeller_mass import Total_Propeller_mass
from Battery_mass import Total_Battery_mass
from Sup_struct_mass import sup_struct_mass
from Skid_mass import Total_Skid_mass
from Fuselage_mass import Total_Fuselage_mass
from Electronics_mass import Total_Electronics_mass

#Estimation of vehicle OEW

class Finaldesign: 
    mass_estimations = 0 
    sup_structure_mass = 0 
    def __init__(self, name):
        self.name = name 
    def Motor(self):
        self.mass_estimations += Total_Motor_mass
    def Propeller(self):
        self.mass_estimations += Total_Propeller_mass

Bonmobile = Finaldesign("bon")

Bonmobile.sup_structure_mass = sup_struct_mass()