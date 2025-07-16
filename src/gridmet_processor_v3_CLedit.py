'''
Takes an nc4 file and reprojects, resamples to match MODIS VI data, calling gdal directly
(espg 4326, pixels size ~ 250m)

using python 3.5.5

script originally by Mike Tercek, last edited by Carolyn Lober 3/12/24
'''

import sys,os,datetime,glob


def set_bounding_box(in_name, out_name): 
    global input_path, output_path
    print('Setting bounding box . . . ')
    print(datetime.datetime.now())
    #del_list = []
    #the_nums = list(range(1,13)) # This is more than needed. One for each layer
    command1 = 'gdalwarp -te -124.75 24.5 -66.9 49.4 -r bilinear -co "BIGTIFF=YES" {i}{in_name} {i}{out_name}'
    
    command1 = command1.format(i = input_path, in_name = in_name, out_name = out_name)
    os.system(command1)
    
    #os.remove(input_path + in_name)
    return out_name

def resample(in_name, out_name):
    global input_path, output_path
    print('Resampling . . .')
    print(datetime.datetime.now())
    command2 = 'gdal_translate -outsize 30330 13055 -co "BIGTIFF=YES" -co "COMPRESS=DEFLATE" {i}{in_fn} {o}{out_fn}'
    
    command2 = command2.format(i = input_path, in_fn = in_name, o = output_path, out_fn = out_name)
    print(command2)
    os.system(command2)
    
    try:
        os.remove(input_path + in_name)
    except:
        print('Failed to remove:', in_name)

if __name__ == "__main__":
    input_path = 'data/raw_data/GRIDmet/wateryear/precip/' 
    output_path = 'data/results/CONUS/gm_proj/bilinear/precip/'
    
    start_name = 'copy1.nc'
    fl = os.listdir(input_path)
    fl.sort()
    
    done_list = os.listdir(output_path)
    done_numbers = [x.split('.')[0] for x in done_list]
    i = 0
    for start_name in fl:
        file_no = start_name.split('.')[0]
        
        
        if file_no in done_numbers:
            print('Skipping: ', start_name)
            continue
        
        
        print(start_name)
        bbox_name = set_bounding_box(start_name, start_name[:-3] + '_bbox.tif')
        resample(bbox_name, start_name[:-3] + '_resampled.tif')
        i +=1
        #if i > 0 : break

