#!/usr/bin/env python3
#Author: Rohit Singh
#Date: 04/01/2019
#convert a .dat file to .csv, .vtk, .vtp and .pvd file for planar view of crack
#use this planar data in crackGen3D.py to calculate 3D view of crack based on crack opening value

import re
import os, sys
import vtk
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from collections import  defaultdict
from collections import OrderedDict

#creates required folder structure to save the files
def createfolder():
    try:
        if os.path.exists("./Crack Data"):
            return
        os.makedirs("./Crack Data")
        os.makedirs("./Crack Data/Plane")
        os.makedirs("./Crack Data/Plane/VTK")
        os.makedirs("./Crack Data/Plane/VTP")

        os.makedirs("./Crack Data/3D_UP")
        os.makedirs("./Crack Data/3D_UP/VTK")
        os.makedirs("./Crack Data/3D_UP/VTP")

        os.makedirs("./Crack Data/3D_DOWN")
        os.makedirs("./Crack Data/3D_DOWN/VTK")
        os.makedirs("./Crack Data/3D_DOWN/VTP")
    except:
        print("Could not create directory for storing output files")


#reads the data from the .dat input file.
#read is dependent on the data formatting in the input file so refer to the sample input file
def readData(file_path):
    with open(file_path, "r") as file:
        data=file.read()
    data_row= data.split("\n")
    data_formatted=""
    for row in data_row:
        row=row.strip(" ")
        row=re.sub("        ", ",", row)
        row = re.sub("    ", ",", row)
        data_formatted+=row+"\n"
    return data_formatted

#saves the input data in .csv format
def convertToCSV(data_formatted, no_of_stages):
    with open('Crack Data\crack_data.csv', mode='w') as file:
        file.write(data_formatted)

#parses thorugh the input data and generates the .vtk format for the input data
def convertToVTK(data_formatted, colormap, stage_to_extract):
    data_list=data_formatted.split("\n")
    header="# vtk DataFile Version 2.0 \n" \
            "Crack Data \n" \
            "ASCII \n" \
            "DATASET POLYDATA \n" \
            "POINTS"

    dataVTK = ""
    dataVTK+=header
    count=0
    stage=0
    skip=0
    no_points=0
    no_polygons=0
    no_triangles=0

    adjacent_faces=defaultdict(set)

    #specifying vertex positions-(x,y,z)
    for row in data_list:
        if row[:4] == "ZONE":
            stage += 1
        if(stage <stage_to_extract):
            skip+=1
        if(stage==stage_to_extract):
            if row[:4]=="ZONE":
                i=row.find("N=")
                no_points=int(row[i+2:i+row[i:].find(",")])+1
                j=row.find("E=")
                b=row[j:].find("\n")
                no_polygons = int(row[j + 2:])
                dataVTK+=" "+str(no_points)+" float \n" + \
                    "0.0000,0.0000,0.0000 \n"
                count=1
            else:
                if(count<no_points):
                    dataVTK+=row[:len(row)-row[::-1] .find(",")-1]+" \n"
                    count+=1
                else:
                    if(count == no_points):
                        dataVTK+="POLYGONS"
                        no_vertices=row.count(",")
                        temp = row.split(",")
                        if no_vertices==2:
                            no_triangles+=1
                            #create a dictionary with all adjacent vertices for the vertex
                            adjacent_faces[int(temp[0].strip())].add(int(temp[1].strip()))
                            adjacent_faces[int(temp[0].strip())].add(int(temp[2].strip()))
                            adjacent_faces[int(temp[1].strip())].add(int(temp[2].strip()))
                            adjacent_faces[int(temp[1].strip())].add(int(temp[0].strip()))
                            adjacent_faces[int(temp[2].strip())].add(int(temp[0].strip()))
                            adjacent_faces[int(temp[2].strip())].add(int(temp[1].strip()))
                        else:
                            adjacent_faces[int(temp[0].strip())].add(int(temp[3].strip()))
                            adjacent_faces[int(temp[0].strip())].add(int(temp[1].strip()))
                            adjacent_faces[int(temp[1].strip())].add(int(temp[0].strip()))
                            adjacent_faces[int(temp[1].strip())].add(int(temp[2].strip()))
                            adjacent_faces[int(temp[2].strip())].add(int(temp[1].strip()))
                            adjacent_faces[int(temp[2].strip())].add(int(temp[3].strip()))
                            adjacent_faces[int(temp[3].strip())].add(int(temp[2].strip()))
                            adjacent_faces[int(temp[3].strip())].add(int(temp[0].strip()))
                        dataVTK+=str(no_vertices+1)+","+row+" \n"
                        count+=1
                    else:
                        no_vertices = row.count(",")
                        temp = row.split(",")
                        if no_vertices==2:
                            no_triangles+=1
                            adjacent_faces[int(temp[0].strip())].add(int(temp[1].strip()))
                            adjacent_faces[int(temp[0].strip())].add(int(temp[2].strip()))
                            adjacent_faces[int(temp[1].strip())].add(int(temp[2].strip()))
                            adjacent_faces[int(temp[1].strip())].add(int(temp[0].strip()))
                            adjacent_faces[int(temp[2].strip())].add(int(temp[0].strip()))
                            adjacent_faces[int(temp[2].strip())].add(int(temp[1].strip()))
                        else:
                            adjacent_faces[int(temp[0].strip())].add(int(temp[3].strip()))
                            adjacent_faces[int(temp[0].strip())].add(int(temp[1].strip()))
                            adjacent_faces[int(temp[1].strip())].add(int(temp[0].strip()))
                            adjacent_faces[int(temp[1].strip())].add(int(temp[2].strip()))
                            adjacent_faces[int(temp[2].strip())].add(int(temp[1].strip()))
                            adjacent_faces[int(temp[2].strip())].add(int(temp[3].strip()))
                            adjacent_faces[int(temp[3].strip())].add(int(temp[2].strip()))
                            adjacent_faces[int(temp[3].strip())].add(int(temp[0].strip()))
                        dataVTK += str(no_vertices+1)+","+ row+" \n"
                        count += 1
            if count==no_points+no_polygons:
                break
    count=0

    #specifying vertex ordering for polygons
    index = dataVTK.find("POLYGONS")
    dataVTK=dataVTK[:(index+8)]+" "+str(no_polygons)+" "+str(no_polygons*5-no_triangles)+"\n" + dataVTK[index+8:]
    dataAtt ="POINT_DATA "+str(no_points)+" \n" \
            "SCALARS point_scalars float 1 \n" \
             "LOOKUP_TABLE my_table"+"\n"
    dataVTK += dataAtt
    dataVTK = re.sub(",", " ", dataVTK)
    scalar_val=[0.0]
    max_val=0

    #specifying scalar values at vertices
    for row in data_list:
        if(count==0):
            dataVTK +="0"+"\n"
            count+=1
        elif(count>skip and count-skip<no_points):
            n=len(row)
            val=float(row[n-row[::-1].find(",")+1:])
            scalar_val.append(val)
            if(val>max_val):
                max_val=val
            dataVTK+=str(val)+ " \n"
            count += 1
        else:
            count+=1

    #perform linear interpolation to determine scalar value at a vertex
    interpolated_scalars=set()
    last_color = 0
    for i in range(0, len(scalar_val)):
        neighbors=adjacent_faces[i]
        sum_scalars=scalar_val[i]
        count=1
        if(len(neighbors)!=0):
            for j in neighbors:
                if(j!=i):
                    sum_scalars+=scalar_val[j]
                    count+=1
            last_color=sum_scalars/(count)
        else:
            pass
        interpolated_scalars.add(last_color)

    #specifying color map for vertices
    dataVTK+="LOOKUP_TABLE my_table"+" "+str(len(interpolated_scalars))+"\n"

    interpolated_scalars=sorted(list(interpolated_scalars))
    max_val=max(interpolated_scalars)
    min_val=min(interpolated_scalars)
    for i in interpolated_scalars:
        color_val=(i-min_val)/(max_val-min_val)
        color=colormap(color_val)
        rgb=str(color[0])+" "+str(color[1])+" "+str(color[2])+" 1.0 \n"
        dataVTK+=rgb

    #store the formatted data in .vtk and .vtp format
    fileno=""
    if(stage_to_extract<10):
        fileno="0"+str(stage_to_extract)
    else:
        fileno=str(stage_to_extract)
    with open('Crack Data\Plane\VTK\crack_data_vtk'+ \
              fileno+'.vtk', mode='w') as file:
        file.write(dataVTK)

    vtkf='Crack Data\Plane\VTK\crack_data_vtk'+fileno+'.vtk'
    vtpf='Crack Data\Plane\VTP\crack_data_vtp'+fileno+'.vtp'
    vtk2vtp(vtkf, vtpf, binary=False)

#source: https://gist.github.com/thomasballinger/1281457
def vtk2vtp(invtkfile, outvtpfile, binary=False):
    """What it says on the label"""
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(invtkfile)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(outvtpfile)
    if binary:
        writer.SetFileTypeToBinary()
    writer.SetInputConnection(reader.GetOutputPort())
    writer.Update()


#link all.vtp files to produce a time varying data set
def createPVDFile(no_of_stages):
    header="<?xml version=\"1.0\"?> \n" \
                "<VTKFile type=\"Collection\" version=\"0.1\">\n" \
                    "  <Collection>\n"

    tail="  </Collection>\n" \
            "</VTKFile>\n"

    data=header

    for i in range(1,no_of_stages+1):
        if i<10:
            data += "    <DataSet timestep=\"" + str(i - 1) + "\" file=\"crack_data_vtp0" + str(i) + ".vtp\"/>\n"
        else:
            data += "    <DataSet timestep=\"" + str(i - 1) + "\" file=\"crack_data_vtp" + str(i) + ".vtp\"/>\n"

    data+=tail

    with open('D:\CS 6635 Vis or Data Sceince\Final Project\Crack Data\Plane\VTP\crack_data.pvd', mode='w') as file:
        file.write(data)


#create a cool to warm color map for the input data set based on max. crack opening
def createColorMap(data_formatted, no_of_stages):
    data_list=data_formatted.split("\n")
    scalar_val=[0]
    skip=0
    stage=0
    count=0
    for row in data_list:
        if row[:4] == "ZONE":
            stage += 1
        if(stage <no_of_stages):
            skip+=1
        if stage==no_of_stages:
            if row[:4] == "ZONE":
                i = row.find("N=")
                no_points = int(row[i + 2:i + row[i:].find(",")]) + 1
                count=1
            else:
                if (count < no_points):
                    n = len(row)
                    val = float(row[n - row[::-1].find(",") + 1:])
                    scalar_val.append(val)
                    count+=1
    max_val=max(scalar_val)
    min_val=min(scalar_val)
    range_val= len(set(scalar_val))
    for i in range(0, range_val):
        scalar_val[i]=(scalar_val[i]-min_val)/(max_val-min_val)
    color=cm.get_cmap('coolwarm', range_val)
    newcolors=color(np.linspace(0, 1, range_val))
    newcmp = ListedColormap(newcolors, name='CoolWarm')
    '''
    fg, ax = plt.subplots()
    psm = ax.pcolormesh(np.array([sorted(scalar_val)]), cmap=newcmp, rasterized=True, vmin=0, vmax=1)
    fg.colorbar(psm, ax=ax)
    plt.show()
    '''
    return newcmp


def main():
    no_of_stages = 98 #total number of stages of crack opening
    createfolder() #create folder structure
    data_formatted=readData("sp2_lc2_test_crack2.dat") #input file
    colormap=createColorMap(data_formatted,no_of_stages) #generate a colormap
    convertToCSV(data_formatted, no_of_stages) #convert to csv(optional)
    for i in range(1,no_of_stages+1):
        convertToVTK(data_formatted, colormap, i) #generate a .vtk and .vtp file for each stage
    createPVDFile(no_of_stages) #create a .pvd file linking all stages


if __name__ == '__main__':
    main()
