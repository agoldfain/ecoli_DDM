def get_difference_fits(filepath, outfolder = None, delta_t = 1, constscalingsave = True, crop=[None,None,None,None,None,None], return_difference = False, return_header = False, abssq = False):
    
    import numpy as np
    import pyfits as ft
    from os.path import exists
    from os import makedirs
    from AaronFunctions import A_imsave

    # gets differences between all images in a series with time lag of delta_t frames
    # outfolder is where to save output images to. If outfolder == None, output images are not saved.
    # delta_t is time lag
    # constscalingsave = True writes images with the same contrast scaling for the whole data set (False rescales each image).
    # crop = (xmin,xmax,ymin,ymax,zmin,zmax) = crop values in x,y,z.
    # returns_divided returns divided images
    # return_header returns data header
    # abssq = True takes absolute value squared of difference

    if outfolder != None:
        if not exists(outfolder):
            makedirs(outfolder)
    
    #open data
    data = ft.open(filepath, memmap=True)[0]
    data_header = data.header
    
    #do cropping
    for i in [0,2,4]:
        if crop[i]==None:
            crop[i] = 0
    if crop[1] == None:
        crop[1] = data.shape[1]
    if crop[3] == None:
        crop[3] = data.shape[2]
    if crop[5] == None:
        crop[5] = data.shape[0]
    data = data.data[ crop[4]:crop[5], crop[0]:crop[1], crop[2]:crop[3] ]


    #initialize difference array
    difference = np.empty( (data.shape[1],data.shape[2], data.shape[0]-delta_t) )

    #loop through each image
    for i in np.arange(0,data.shape[0]-delta_t):
        result = data[i+delta_t,:,;] - data[i,:,:]
        
        if abssq:
            result = abs(result)**2
            
        difference[:,:,i] = np.array(result)
        
    

    if constscalingsave:
        im_max = difference.max()
        im_min = difference.min()
    else:
        im_max = None
        im_min = None
        
    if outfolder != None:
        for i in np.arange(0,difference.shape[-1]):
            A_imsave(outfolder + 'difference_image' + str(i).zfill(5) +'.tif', difference[:,:,i], im_max, im_min)

    if return_difference:
        if return_header:
            return [difference, data_header]
        else:
            return difference
    else:
     if return_header:
        return data_header
                   


