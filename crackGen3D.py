#!/usr/bin/env python3
#!/usr/bin/env python
#Author: Rohit Singh
#Date: 04/01/2019
#use the planar crack data from crackGenPlane.py to calculate 3D view of crack based on crack opening value
#save the data for up and down plane in .csv, .vtk, .vtp and .pvd files for 3D view of crack

import re
import os, sys
import vtk
import numpy as np
from collections import defaultdict
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from collections import  defaultdict
from collections import OrderedDict

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

#calculates up and down plane coordinates from planar crack coordinates using vertex normals
#stores the calculated coordinates in .csv file
#mode=UP for Upper plane and mode=DOWN for lower plane
def crackGen3d(filepath,no_of_stages, mode):
    data=""
    data_final=""
    data_list=[]
    data_list_final=[]
    count=0
    stage_number=0
    next_stage=0
    index=0
    no_points=0
    no_polygons=0
    no_triangles=0

    with open(filepath, 'r') as file:
        data=file.read()
    data_list=data.split("\n")


    for row in data_list:
        if(index==next_stage):
            if row[:4] == "ZONE":
                data_list_final.append(row)

                stage_number+=1
                normals_dict = defaultdict(list) #dictionary with normals at the verices
                points_dict = defaultdict(list)
                calculateNormals(data_list, normals_dict, points_dict, stage_number)
                for key, values in normals_dict.items():
                    sum = np.array([0, 0, 0])
                    for i in values:
                        sum = np.add(sum, i)
                    sum = sum / len(values)
                    normals_dict[key] = sum

                i = row.find("N=")
                no_points = int(row[i + 2:i + row[i:].find(",")])
                j = row.find("E=")
                b = row[j:].find("\n")
                no_polygons = int(row[j + 2:])
                count = 0
                next_stage+=1
            else:
                if (count < no_points):
                    count += 1
                    row_val=row.split(",")

                    for i in range(len(row_val)):
                        row_val[i]=float(row_val[i].strip())
                    #calculate the UP coordinates from normal new_a^=a^+crack_val.n^/2
                    if(mode.upper()=="UP"):
                        if(normals_dict[count]!=[]):
                            new_val=np.add(row_val[:-1], np.multiply \
                                (normals_dict[count], row_val[3]/2))
                            row_val[0]=new_val[0]
                            row_val[1]=new_val[1]
                            row_val[2]=new_val[2]
                    # calculate the DOWN coordinates from normal new_a^=a^-crack_val.n^/2
                    elif(mode.upper()=="DOWN"):
                        if (normals_dict[count] != []):
                            new_val = np.subtract(row_val[:-1], np.multiply \
                                (normals_dict[count], row_val[3]/2))
                            row_val[0] = new_val[0]
                            row_val[1] = new_val[1]
                            row_val[2] = new_val[2]
                    data_list_final.append((",").join(str(i) for i in row_val))
                    next_stage += 1
            if count == no_points:
                next_stage+=no_polygons
        elif(index<next_stage):
            data_list_final.append(row)
        index+=1

    data_final=("\n").join(data_list_final[:-1])
    with open('Crack Data\crack_data3D_'+mode+'.csv', mode='w') as file:
        file.write(data_final)

    return data_final


#method to calculate normal at all vertex based on Right-hand rule for vector cross product
def calculateNormals(data_list, normals_dict, points_dict, stage_number):
    count=0
    index=0
    stage=0
    no_points=0
    no_polygons=0
    no_triangles=0

    for row in data_list:
        if row[:4] == "ZONE":
            i = row.find("N=")
            no_points = int(row[i + 2:i + row[i:].find(",")])
            j = row.find("E=")
            b = row[j:].find("\n")
            no_polygons = int(row[j + 2:])
            count = 0
            stage+=1

        if(stage==stage_number):
            if 0<count<=no_points:
                point_coords=row.split(",")
                for i in point_coords:
                    points_dict[count].append(float(i.strip()))

            elif count>no_points:
                row_val = row.split(",")
                for i in row_val:
                    if i.strip().isnumeric()==False:
                        row_val.remove(i)
                for i in range(-1, len(row_val)-1):
                    point_min1=row_val[i-1]
                    point=row_val[i]
                    point_plus1=row_val[i+1]

                    point_min1=int(point_min1.strip())
                    point=int(point.strip())
                    point_plus1=int(point_plus1.strip())

                    cross_prod=np.cross(np.subtract(points_dict[point][:-1], \
                                                                    points_dict[point_min1][:-1]),\
                                                        np.subtract(points_dict[point_plus1][:-1],\
                                                                    points_dict[point][:-1]))
                    if(np.linalg.norm(cross_prod)!=0):
                        normals_dict[point].append(cross_prod/np.linalg.norm(cross_prod))
        count += 1


#parses thorugh the formatted 3D data and generates the .vtk format for UP and DOWN plane
def convertToVTK(data_formatted, colormap, stage_to_extract, mode):
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
    with open('Crack Data\\3D_'+mode+'\VTK\crack_data_3D_vtk'+ \
              fileno+'.vtk', mode='w') as file:
        file.write(dataVTK)

    vtkf='Crack Data\\3D_'+mode+'\VTK\crack_data_3D_vtk'+fileno+'.vtk'
    vtpf='Crack Data\\3D_'+mode+'\VTP\crack_data_3D_vtp'+fileno+'.vtp'
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
def createPVDFile(no_of_stages, mode):
    header="<?xml version=\"1.0\"?> \n" \
                "<VTKFile type=\"Collection\" version=\"0.1\">\n" \
                    "  <Collection>\n"

    tail="  </Collection>\n" \
            "</VTKFile>\n"

    data=header

    for i in range(1,no_of_stages+1):
        if i<10:
            data += "    <DataSet timestep=\"" + str(i - 1) + "\" file=\"crack_data_3D_vtp0" + str(i) + ".vtp\"/>\n"
        else:
            data += "    <DataSet timestep=\"" + str(i - 1) + "\" file=\"crack_data_3D_vtp" + str(i) + ".vtp\"/>\n"

    data+=tail

    with open('Crack Data\\3D_'+mode+'\VTP\crack_data_3D_'+mode+'.pvd', mode='w') as file:
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

if __name__ == '__main__':
    main()
