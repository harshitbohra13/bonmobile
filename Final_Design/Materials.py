#fibers
# carbon_fibre = {"strength":(1730+7060)/2, "E":(155+860)/2, "density":(1740+2190)/2 }
carbon_fibre_hexcel= {"strength":4930, "E":243, "density":1790 }
flax_fibre = {"strength":(343+2000)/2, "E":(27.6+103)/2, "density":(1400+1500)/2 }
bamboo_fibre = {"strength":(140+800)/2, "E":(11+32)/2, "density":(600+1100)/2 }
Ramie_fibre = {"strength":(400+1000)/2, "E":(24.5+128)/2, "density":(1000+1550)/2 }

#resin
Vertimer_Mallinda = {"strength":135 , "E":3.42 , "density":1050} #Flexural strength 135, Flexural Modulus 3.42
PLA = {"strength":(48+60)/2 , "E":(3.45+3.83)/2, "density":(1210+1250)/2}
Epoxy = {"strength":(36+71)/2 , "E":(2.35+3.08)/2, "density":(1110+1140)/2}
#Percentage of fiber content
Fibre_perc_cont = 0.6 # 47% carbon fiber mallinda

#pick fibre and resin
fibre = carbon_fibre_hexcel
resin = Epoxy


#Calcualting composite
comp_strength = Fibre_perc_cont * fibre["strength"] + (1-Fibre_perc_cont) * resin["strength"]
comp_E = Fibre_perc_cont * fibre["E"] + (1-Fibre_perc_cont) * resin["E"]
comp_density = Fibre_perc_cont * fibre["density"] + (1-Fibre_perc_cont) * resin["density"] 

#Printing values
print("Composite stength: ", comp_strength, "MPa")
print("Composite E modulus: ",comp_E, "GPa")
print("Composite density: ",comp_density, "kg/m^3")

#mallinda composite properties T130:
# Tensile strnegth: 1722.5 MPa
# compressive strnegth: 895.7 MPa
# E modulus: 113.685 GPa

#How recyling changes the values
#directionality