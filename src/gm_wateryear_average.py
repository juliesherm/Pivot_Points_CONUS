import xarray as xr
import pandas as pd
import os
import numpy as np

def process_two_years(file_paths, my_year):
    # Load data from two consecutive NetCDF files
    datasets = [xr.open_dataset(path) for path in file_paths]
    
    # Concatenate datasets along the time dimension
    full_dataset = xr.concat(datasets, dim='day')
    
    # Convert time index to pandas datetime index
    full_dataset['day'] = pd.to_datetime(full_dataset.day.values)
    
    # Group data by water year
    full_dataset = full_dataset.assign_coords(water_year = full_dataset.day.dt.year.where(full_dataset.day.dt.month < 10, full_dataset.day.dt.year + 1))
    
    # select only my wateryear
    full_dataset = full_dataset.sel(day = full_dataset.water_year == np.int64(my_year))
    
    # Average data for each water year
    water_year_avg = full_dataset.groupby('water_year').sum()
    
    return water_year_avg

def main(input_dir, output_dir):
    # List all NetCDF files in the input directory
    files = [f for f in os.listdir(input_dir) if f.endswith('.nc')]
    files.sort()  # Ensure files are processed in chronological order
    
    # Process files in pairs
    for i in range(len(files) - 1):  # Avoid going out of index range
        print(f"Processing {files[i]} and {files[i+1]}...") 
        file_paths = [os.path.join(input_dir, files[i]), os.path.join(input_dir, files[i+1])]
        #get wateryear 
        my_year = files[i+1].split('.')[0].split('_')[1]
        
        water_year_avg = process_two_years(file_paths, my_year)
        
        # Save averaged data to new NetCDF files for each water year
        output_path = os.path.join(output_dir, files[i+1])
        water_year_avg.to_netcdf(output_path)
        print(f"Saved averaged data for water year {my_year} to {output_path}")

# call functions
input_directory = 'data/raw_data/GRIDmet/precip'
output_directory = 'data/raw_data/GRIDmet/wateryear/precip'
main(input_directory, output_directory)
