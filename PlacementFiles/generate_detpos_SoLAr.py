xlen = 100
ylen = 100
distance = 1


File_SiPM = open(r"SiPM", "w")
File_Pixel = open(r"Pixel", "w")


def place_array(center, num, FilePix, FileSiPM):
    # Get the center coordinate of the overall array (5 Pixel, 1 SiPM)
    # Adjust the SiPM and the 5 pixels with the correct spacing
    # Drawing of File in Docs
    numPix, numSiPM = num
    x = center[0]
    y = center[1]
    SiPMCenter = [x+0.15, y+0.15]
    Pi1Center = [x-0.35, y-0.35]
    Pi2Center = [x-0.35, y-0.017]
    Pi3Center = [x-0.35, y+0.316]
    Pi4Center = [x-0.017, y-0.35]
    Pi5Center = [x+0.316, y-0.35]
    # Write the SiPM and the 5 pixels to the corresponding files using the orientation along z and type 1
    FileSiPM.write(str(numSiPM)+ " " + str(SiPMCenter[0]) + " " + str(SiPMCenter[1]) + " " + str(0) + " " + str(1) + " " + str(3) +"\n")
    numSiPM += 1
    FilePix.write(str(numPix)+ " " + str(Pi1Center[0]) + " " + str(Pi1Center[1]) + " " + str(0) + " " + str(1) + " " + str(3) +"\n")
    numPix += 1
    FilePix.write(str(numPix)+ " " + str(Pi2Center[0]) + " " + str(Pi2Center[1]) + " " + str(0) + " " + str(1) + " " + str(3) +"\n")
    numPix += 1
    FilePix.write(str(numPix)+ " " + str(Pi3Center[0]) + " " + str(Pi3Center[1]) + " " + str(0) + " " + str(1) + " " + str(3) +"\n")
    numPix += 1
    FilePix.write(str(numPix)+ " " + str(Pi4Center[0]) + " " + str(Pi4Center[1]) + " " + str(0) + " " + str(1) + " " + str(3) +"\n")
    numPix += 1
    FilePix.write(str(numPix)+ " " + str(Pi5Center[0]) + " " + str(Pi5Center[1]) + " " + str(0) + " " + str(1) + " " + str(3) +"\n")
    numPix += 1
    num = [numPix, numSiPM]
    return num


# Determine the number of placements in x and y direction
x_num = xlen / distance
y_num = ylen / distance

current_number = [0, 0]
i = j = 0;
while( i<x_num ):
    j = 0
    while( j<y_num ):
        current_number = place_array([distance/2 + i*distance,distance/2 + j * distance], current_number, File_Pixel, File_SiPM)
        j += 1
    i += 1

File_SiPM.close()
File_Pixel.close()
