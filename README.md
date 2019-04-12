# Crack-_Mesh_DAT_PVD
This project shows how to convert a .dat file into .csv, .vtk, .vtp and .pvd file compatible for time varying visualization in Para View. A sample input .dat file format is included for reference. 

There are two .py files crackGenPlane.py and crackGen3D.py.

crackGenPlane.py actually reads the .dat input file and creates the corresponding .csv, .vtk, .vtp and .pvd files for planar view of crack as present in the input file.

crackGen3D.py reads the output .csv file from crackGenPlane.py to create a 3D view of the crack based on vertex normals.

Parameters to run crackGenPlane.py:

    def main():
        no_of_stages = 98
        createfolder()
        convertToCSV(no_of_stages)
        createPVDFile(no_of_stages)
    
Parameters to run crackGen3D.py:

    def main():
        no_of_stages = 98
        mode = "UP"
        crackGen3d(no_of_stages, mode)
        createPVDFile(no_of_stages, mode)
        mode = "DOWN"
        crackGen3d(no_of_stages, mode)
        createPVDFile(no_of_stages, mode)
    
crackGen3D.py reads the csv file output from crackGenPlane.py so you will need to run crackGenPlane.py first and then run crackGen3D.py.
