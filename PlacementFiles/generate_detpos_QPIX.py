xlen = 1400
ylen = 600
zlen = 360

distance = 10

x_num = xlen / distance
y_num = ylen / distance
z_num = zlen / distance
i = j = 0;
num = 0


# # z= 360 plane
# while( i<x_num ):
#     j = 0
#     while( j<y_num ):
#         print( num, distance/2 + i*distance, distance/2 + j*distance, zlen , 1, 3 )
#         j += 1
#         num += 1
#     i += 1

#z= 0 plane
i = j = 0;
while( i<x_num ):
    j = 0
    while( j<y_num ):
        print( num, distance/2 + i*distance, distance/2 + j*distance, 0 , 1, 3 )
        j += 1
        num += 1
    i += 1


# # y = 0 plane
# i = j = 0;
# while( i<x_num ):
#     j = 0
#     while( j<z_num ):
#         print( num, distance/2 + i*distance, 0 , distance/2 + j*distance, 1, 2 )
#         j += 1
#         num += 1
#     i += 1

# # # y = 600  plane
# i = j = 0;
# while( i<x_num ):
#     j = 0
#     while( j<z_num ):
#         print( num, distance/2 + i*distance, ylen , distance/2 + j*distance, 1, 2 )
#         j += 1
#         num += 1
#     i += 1


# # x = 0 plane
# i = j = 0;
# while( i<y_num ):
#     j = 0
#     while( j<z_num ):
#         print( num, 0, distance/2 + i*distance , distance/2 + j*distance, 1, 1 )
#         j += 1
#         num += 1
#     i += 1

# # # x = 230  plane
# i = j = 0;
# while( i<y_num ):
#     j = 0
#     while( j<z_num ):
#         print( num, xlen,  distance/2 + i*distance, distance/2 + j*distance, 1, 1 )
#         j += 1
#         num += 1
#     i += 1
