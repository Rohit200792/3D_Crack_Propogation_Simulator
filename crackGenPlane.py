#!/usr/bin/env python
import re
import os, sys
import vtk

def convertToCSV(no_of_stages):
    with open("sp2_lc2_test_crack2.dat", "r") as file:
        data=file.read()

    data_row= data.split("\n")
    data_formatted=""

    for row in data_row:
        row=row.strip(" ")
        row=re.sub("        ", ",", row)
        row = re.sub("    ", ",", row)
        data_formatted+=row+"\n"

    with open('Crack Data\crack_data.csv', mode='w') as file:
        file.write(data_formatted)

    for i in range(1,no_of_stages+1):
        convertToVTK(data_formatted, i)

def convertToVTK(data_formatted, stage_to_extract):
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
            row=float((row.split()[3]).strip())
            row=row*scale
            rgb=str(row)+" "+"0.0"+" "+str(row)+" 1.0 \n"
            dataVTK+=rgb
        count+=1

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



no_of_stages=98
convertToCSV(no_of_stages)
createPVDFile(no_of_stages)