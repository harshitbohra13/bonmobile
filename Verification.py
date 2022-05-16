from re import X
import unittest
import drag_calc as drag
import Noise_Calculations as noise


class testclass(unittest.TestCase):
    def test_list_drag(self):
        """
        This is a test check if the function rotordrag passes the correct result
        """
        data1 = 0
        data2 = 0
        data3 = 0
        result = drag.rotor_drag(data1, data2, data3)
        self.assertEqual(result, 0.5)
    
    def test_list_drag90(self):
        """
        This is a test check if the function rotordrag90 passes the correct result
        """
        data1 = 0
        data2 = 0
        result = drag.rotor_drag90(data1, data2)
        self.assertEqual(result, 10)

    def test_list_noise(self):
        """
        This is a test check if the function rotordrag90 passes the correct result
        """
        data1 = 10
        data2 = 100
        data3 = 10
        result = noise.rotor_solidity(data1, data2, data3)
        self.assertEqual(result, 10)
    
    def test_list_noise(self):
        """
        This is a test check if the function rotordrag90 passes the correct result
        """
        data1 = 10
        data2 = 100
        data3 = 10
        data4 = 10
        data5 = 10
        data6 = 10
        result = noise.tip_velocity(data1,data2,data3,data4,data5,data6)
        self.assertEqual(result, 10)
    

        
if __name__ == '__main__':
    unittest.main()
    

