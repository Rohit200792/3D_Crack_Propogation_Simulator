#!/usr/bin/env python
import re
import os, sys
import vtk
import numpy as np
from collections import defaultdict


def crackGen3d(no_of_stages, mode):
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

    with open('Crack Data\crack_data.csv', 'r') as file:
        data=file.read()
    data_list=data.split("\n")



    for row in data_list:
        if(index==next_stage):
            if row[:4] == "ZONE":
                data_list_final.append(row)

                stage_number+=1
                normals_dict = defaultdict(list)
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

                    if(mode.upper()=="UP"):
                        if(normals_dict[count]!=[]):
                            new_val=np.add(row_val[:-1], np.multiply \
                                (normals_dict[count], row_val[3]/2))
                            row_val[0]=new_val[0]
                            row_val[1]=new_val[1]
                            row_val[2]=new_val[2]
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
    with open('Crack Data\crack_data3D_'+folder_type+'.csv', mode='w') as file:
        file.write(data_final)

    for i in range(1,no_of_stages+1):
        convertToVTK(data_final, i, folder_type)

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
                #print(row_val)
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



def convertToVTK(data_formatted, stage_to_extract, folder_type):
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
                        if no_vertices==2:
                            no_triangles+=1
                        dataVTK+=str(no_vertices+1)+","+row+" \n"
                        count+=1
                    else:
                        no_vertices = row.count(",")
                        if no_vertices==2:
                            no_triangles+=1
                        dataVTK += str(no_vertices+1)+","+ row+" \n"
                        count += 1
            if count==no_points+no_polygons:
                break
    count=0
    index = dataVTK.find("POLYGONS")
    dataVTK=dataVTK[:(index+8)]+" "+str(no_polygons)+" "+str(no_polygons*5-no_triangles)+"\n" + dataVTK[index+8:]

    dataAtt ="POINT_DATA "+str(no_points)+" \n" \
            "SCALARS point_scalars float 1 \n" \
             "LOOKUP_TABLE my_table"+"\n"
    dataVTK += dataAtt
    dataVTK = re.sub(",", " ", dataVTK)
    max_val=0
    min_val=0
    for row in data_list:
        if(count==0):
            dataVTK +="0"+"\n"
            count+=1
        elif(count>skip and count-skip<no_points):
            n=len(row)
            val=float(row[n-row[::-1].find(",")+1:])
            if(val>max_val):
                max_val=val
            dataVTK+=str(val)+ " \n"
            count += 1
        else:
            count+=1

    dataVTK+="LOOKUP_TABLE my_table"+" "+str(no_points)+"\n"
    scale=1/max_val
    count = 0
    for row in data_list:
        if count==0:
            dataVTK+="0.0 0.0 0.0 1.0 \n"
        elif(count>skip and count-skip<no_points):
            print(row)
            row=float((row.split(",")[3]).strip())
            row=row*scale
            rgb=str(row)+" "+"0.0"+" "+str(row)+" 1.0 \n"
            dataVTK+=rgb
        count+=1

    fileno=""
    if(stage_to_extract<10):
        fileno="0"+str(stage_to_extract)
    else:
        fileno=str(stage_to_extract)
    with open('Crack Data\\3D_'+folder_type+'\VTK\crack_data_3D_vtk'+ \
              fileno+'.vtk', mode='w') as file:
        file.write(dataVTK)

    vtkf='Crack Data\\3D_'+folder_type+'\VTK\crack_data_3D_vtk'+fileno+'.vtk'
    vtpf='Crack Data\\3D_'+folder_type+'\VTP\crack_data_3D_vtp'+fileno+'.vtp'
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


def createPVDFile(no_of_stages, folder_type):
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

    with open('Crack Data\\3D_'+folder_type+'\VTP\crack_data_3D_'+folder_type+'.pvd', mode='w') as file:
        file.write(data)

no_of_stages=98
folder_type="UP"
crackGen3d(no_of_stages,folder_type)
createPVDFile(no_of_stages, folder_type)

folder_type="DOWN"
crackGen3d(no_of_stages,folder_type)
createPVDFile(no_of_stages, folder_type)