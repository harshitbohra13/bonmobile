#importing different subsytems' mass
from Motor_mass import Total_Motor_mass
from Propeller_mass import Total_Propeller_mass
from Battery_mass import Total_Battery_mass
from Sup_struct_mass import Total_Supporting_structure_mass
from Skid_mass import Total_Skid_mass
from Fuselage_mass import Total_Fuselage_mass
from Electronics_mass import Total_Electronics_mass

#Estimation of vehicle OEW
Total_vehicle_OEW = Total_Motor_mass + Total_Propeller_mass + Total_Battery_mass + Total_Supporting_structure_mass + Total_Skid_mass + Total_Fuselage_mass + Total_Electronics_mass
print("Total vehicle OEW: ",Total_vehicle_OEW)

