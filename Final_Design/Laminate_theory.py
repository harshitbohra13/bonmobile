from composites import laminated_plate
import numpy as np
laminaprop = (255, 3, 0.15, 2, 2, 2) #?
plyt = 0.1 #?
stack = [0, 90, 90, 0 ] #?
plate = laminated_plate(stack, plyt=plyt, laminaprop=laminaprop)
# print(plate.ABD)


LOADS = np.array([[1],[1],[1],[1],[1],[1]]) #?
# print(LOADS)
strains  = np.dot(np.linalg.inv(plate.ABD), LOADS)
print(strains)
