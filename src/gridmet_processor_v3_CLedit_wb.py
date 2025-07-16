'''
Takes an nc4 file and reprojects, resamples to match MODIS VI data, calling gdal directly
(espg 4326, pixels size ~ 250m)

using python 3.5.5

script originally by Mike Tercek, last edited by Carolyn Lober 3/12/24
'''

import sys,os,datetime,glob

def reproject(fn, out_fn = 'nps_maca_reprojected.tif'):
    global input_path, output_path
    print('Reprojecting.')
    print(datetime.datetime.now())
    #reproject_command = 'gdalwarp -t_srs "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" -of netCDF -co "FORMAT=NC4" -multi {i}{fn} {i}{out_fn}'.format(fn = fn, out_fn = out_fn, i = input_path)
    reproject_command = 'gdalwarp -t_srs "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0" -multi -co "BIGTIFF=YES" {i}{fn} {i}{out_fn}'.format(fn = fn, out_fn = out_fn, i = input_path)
    print(reproject_command)
    os.system(reproject_command)
    #os.remove(input_path + fn)
    return out_fn

def set_bounding_box(in_name, out_name): 
    global input_path, output_path
    print('Setting bounding box . . . ')
    print(datetime.datetime.now())
    #del_list = []
    #the_nums = list(range(1,13)) # This is more than needed. One for each layer
    command1 = 'gdalwarp -s_srs "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0" -te -124.75 24.5 -66.9 49.4 -co "BIGTIFF=YES" {i}{in_name} {i}{out_name}'
    
    command1 = command1.format(i = input_path, in_name = in_name, out_name = out_name)
    os.system(command1)
    
    ''' 
    for num in the_nums:
        pre_command = 'gdal_translate -b {num} {i}{in_name} {i}tempfoo.tif'.format(num = num, in_name = in_name, i = input_path)
        print(pre_command)
        os.system(pre_command)
        this_command = command1.format(in_name = 'tempfoo.tif' ,out_name = out_name, num2 = num, i = input_path)
        print(this_command)
        os.system(this_command)
        try:
            os.remove(input_path + 'tempfoo.tif')
            del_name = out_name + str(num)
            del_list.append(del_name)
        except:
            continue
    '''
    
    os.remove(input_path + in_name)
    return out_name

def resample(in_name, out_name):
    global input_path, output_path
    print('Resampling . . .')
    print(datetime.datetime.now())
    command2 = 'gdal_translate -outsize 30330 13055 -co "BIGTIFF=YES" -co "COMPRESS=DEFLATE" {i}{in_fn} {o}{out_fn}'
    
    command2 = command2.format(i = input_path, in_fn = in_name, o = output_path, out_fn = out_name)
    os.system(command2)
    
    try:
        os.remove(input_path + in_name)
    except:
        print('Failed to remove:', in_name)
    
    '''
    for fn in del_list:
        this_outname = fn + '_resampled.tif'
        this_command = command2.format(in_fn = fn, out_fn = this_outname, i = input_path, o = output_path)
        print(this_command)
        os.system(this_command)
        try:
            #print("not removing file")
            os.remove(input_path + fn + '.tif') 
        except:
            print('Failed to remove:', fn)
            continue
    '''
  
if __name__ == "__main__":
    web = False
    if web == False: # change me for different folders/variables
        input_path = 'data/results/CONUS/gs_dates/' 
        output_path = 'data/results/CONUS/gs_dates/projected/'
    else:
        input_path = './'
        output_path = './'
        #output_path = './results/'
    start_name = 'copy1.nc'
    if "gs_dates" in input_path:
        file_pattern = os.path.join(input_path, '*20*.tif')
        fl = glob.glob(file_pattern)
        fl = [f.split('/')[4] for f in fl]
    else:
        fl = os.listdir(input_path)
    
    done_list = os.listdir(output_path)
    done_numbers = [x.split('_')[3] for x in done_list]
    i = 0
    for start_name in fl:
        if "gs_dates" not in input_path:
            if start_name.split('.')[-1] != 'nc4' : continue
            file_no = start_name.split('_')[3]
        
        
            if file_no in done_numbers:
                print('Skipping: ', start_name)
                continue
        
        
        print(start_name)
        #ft_name = make_time_unlimited(start_name)
        #converted_name = convert_to_nc3(ft_name)
        reprojected_name = reproject(start_name, start_name[:-4] + "_reprojected.tif")
        bbox_name = set_bounding_box(reprojected_name, start_name[:-4] + '_reprojected_with_extent.tif')
        resample(bbox_name, start_name[:-4] + '_resampled.tif')
        i +=1
        #if i > 0 : break

