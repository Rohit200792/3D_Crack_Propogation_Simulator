# 3D_Crack_Propogation_Simulator
A 3D Visualization for Crack Propagation using Boundary Element Method.

## Overview
This source code shows how to convert a 2D planar data for crack boundary into a corresponding 3D view. It uses two python scripts to convert a .dat file obtained from a crack simulationn into corresponding .vtk, .vtp and .pvd file compatible for time varying visualization in Para View. A sample input .dat file format is included for reference.

## Links
[1. Project Report](./Docs/Project_Report.pdf)

[2. Project Presentation](./Docs/Project_Presentation.pdf)

## API Description
There are two .py files crackGenPlane.py and crackGen3D.py.

crackGenPlane.py actually reads the .dat input file and creates the corresponding .csv, .vtk, .vtp and .pvd files for planar view of crack as present in the input file.

crackGen3D.py reads the output .csv file from crackGenPlane.py to create a 3D view of the crack based on vertex normals.

### Parameters to run crackGenPlane.py:

    def main():
        no_of_stages = 98 #total number of stages of crack opening
        createfolder() #create folder structure
        data_formatted=readData("sp2_lc2_test_crack2.dat") #input file
        colormap=createColorMap(data_formatted,no_of_stages) #generate a colormap
        convertToCSV(data_formatted, no_of_stages) #convert to csv(optional)
        for i in range(1,no_of_stages+1):
            convertToVTK(data_formatted, colormap, i) #generate a .vtk and .vtp file for each stage
        createPVDFile(no_of_stages) #create a .pvd file linking all stages)
    
### Parameters to run crackGen3D.py:

    def main():
        no_of_stages = 98 #total number of stages of crack opening
        filepath='Crack Data\crack_data.csv' #input .csv file for planar crack data
        mode = "UP" #3D upper plane
        data_formatted=readData(filepath)
        colormap=createColorMap(data_formatted,no_of_stages)

        data_final=crackGen3d(filepath,no_of_stages, mode) #generate coordinates for upper plane in 3D view
        for i in range(1,no_of_stages+1):
            convertToVTK(data_final, colormap, i, mode) #generate a .vtk and .vtp file for each stage
        createPVDFile(no_of_stages, mode) #create a .pvd file linking all stages

        mode = "DOWN" #3D lowerplane
        data_final=crackGen3d(filepath,no_of_stages, mode) #generate coordinates for lower plane in 3D view
        for i in range(1,no_of_stages+1):
            convertToVTK(data_final, colormap, i, mode) #generate a .vtk and .vtp file for each stage
        createPVDFile(no_of_stages, mode) #create a .pvd file linking all stages
    
crackGen3D.py reads the csv file output from crackGenPlane.py so you will need to run crackGenPlane.py first and then run crackGen3D.py.

## Sample Output on Paraview:

### Planar Crack from Simulation Data:

![Sample Output](./Docs/SS1.PNG)

### Reconstructed 3D View for  Initial State of Crack:

![Sample Output](./Docs/SS2.PNG)

### Reconstructed 3D View for  Final State of Crack Within Material:

![Sample Output](./Docs/SS3.PNG)

### Crack Opening in Final Stage:

![Sample Output](./Docs/SS4.PNG)


