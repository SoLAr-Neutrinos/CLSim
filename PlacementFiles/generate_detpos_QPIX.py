xlen = 1400
ylen = 600
zlen = 360

pixel_size = 10

x_num = xlen / pixel_size
y_num = ylen / pixel_size
z_num = zlen / pixel_size
i = j = 0;
num = 0

# # z= 360 plane
# while( i<x_num ):
#     j = 0
#     while( j<y_num ):
#         print( num, pixel_size/2 + i*pixel_size, pixel_size/2 + j*pixel_size, zlen , 1, 3 )
#         j += 1
#         num += 1
#     i += 1

#z= 0 plane
i = j = 0;
while( i<x_num ):
    j = 0
    while( j<y_num ):
        print( num, pixel_size/2 + i*pixel_size, pixel_size/2 + j*pixel_size, 0 , 1, 3 )
        j += 1
        num += 1
    i += 1


# # y = 0 plane
# i = j = 0;
# while( i<x_num ):
#     j = 0
#     while( j<z_num ):
#         print( num, pixel_size/2 + i*pixel_size, 0 , pixel_size/2 + j*pixel_size, 1, 2 )
#         j += 1
#         num += 1
#     i += 1

# # # y = 600  plane
# i = j = 0;
# while( i<x_num ):
#     j = 0
#     while( j<z_num ):
#         print( num, pixel_size/2 + i*pixel_size, ylen , pixel_size/2 + j*pixel_size, 1, 2 )
#         j += 1
#         num += 1
#     i += 1


# # x = 0 plane
# i = j = 0;
# while( i<y_num ):
#     j = 0
#     while( j<z_num ):
#         print( num, 0, pixel_size/2 + i*pixel_size , pixel_size/2 + j*pixel_size, 1, 1 )
#         j += 1
#         num += 1
#     i += 1

# # # x = 230  plane
# i = j = 0;
# while( i<y_num ):
#     j = 0
#     while( j<z_num ):
#         print( num, xlen,  pixel_size/2 + i*pixel_size, pixel_size/2 + j*pixel_size, 1, 1 )
#         j += 1
#         num += 1
#     i += 1
