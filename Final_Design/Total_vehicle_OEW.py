#importing different subsytems' mass
from Fuselage_mass import frame
from configuration_and_materials import configure

#Estimation of vehicle OEW
Bonmobile = configure("name")
print(frame(Bonmobile).Rotor_number)
