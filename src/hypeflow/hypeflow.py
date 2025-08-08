# Load needed packages
import xarray as xr
import pint_xarray
import pint
import glob
import netCDF4 as nc4
import os
import cdo
import pandas as pd
from   easymore import Easymore
import numpy       as      np
import geopandas   as      gpd
import sys
from   itertools   import  product
import datetime
from alive_progress import alive_bar #progress bar


# sort geodata from upstream to downstream
def sort_geodata(geodata):
    # find the subbasin order
    geodata['n_ds_subbasins'] = 0
    for index, row in geodata.iterrows():
        n_ds_subbasins = 0
        nextID = row['maindown']
        nextID_index = np.where(geodata['subid']==nextID)[0]

        while (len(nextID_index>0) or nextID > 0) :
            n_ds_subbasins += 1
            nextID = geodata['maindown'][nextID_index[0]]
            nextID_index = np.where(geodata['subid']==nextID)[0]
            
        geodata.loc[index, 'n_ds_subbasins']= n_ds_subbasins

    # Sort the DataFrame by 'n_ds_subbasins' in descending order
    geodata = geodata.sort_values(by='n_ds_subbasins', ascending=False, ignore_index=True)
    geodata = geodata.drop(columns=['n_ds_subbasins'])
    return geodata


# write HYPE forcing from easymore nc files
def write_hype_forcing(easymore_output, timeshift, forcing_units, geofabric_mapping, path_to_save):
    """Initialize forcing files."""
    # check if path_to_save exists and if not create it
    if not os.path.exists(path_to_save):
        os.makedirs(path_to_save)
    # initalize pint registry
    _ureg = pint.UnitRegistry(force_ndarray_like=True)

    # read a list of netcdf files
    easymore_output_list = sorted(glob.glob(os.path.join(easymore_output, '*.nc*')))

    # read the forcing files using xarray and create a self.forcing
    datasets = [xr.open_dataset(p) for p in easymore_output_list]
    ds = xr.concat(datasets, 'time')

    # adjust the model time zone
    ds = ds.assign_coords({
            'time': ds.time.to_index().tz_localize('UTC').tz_convert('America/Edmonton').tz_localize(None)
        })

    # rename the self.forcing variables to match the forcing_vars
    # and assign pint units to the self.forcing variables
    rename_vars_dict = {}
    for key, value_dict in forcing_units.items():
        rename_vars_dict[value_dict['in_varname']] = key
    ds = ds.rename(rename_vars_dict)

    # drop the variable not in self.forcing_vars
    ds = ds[['temperature', 'precipitation']]

    # assign pint units to the self.forcing variables
    renamed_forcing_units = {}
    for key, value_dict in forcing_units.items():
        renamed_forcing_units[key] = value_dict['in_units']
    ds = ds.pint.quantify(units=renamed_forcing_units, unit_registry=_ureg)

    # convert the self.forcing units to the default forcing units
    renamed_to_forcing_units = {}
    for key, value_dict in forcing_units.items():
        renamed_to_forcing_units[key] = value_dict['out_units']
    ds = ds.pint.to(units=renamed_to_forcing_units)

    # after unit conversion, dequantify the self.forcing
    ds = ds.pint.dequantify()

    # print Tobs.txt, TMINobs.txt, TMAXobs.txt, and Pobs.txt files
    # first Pobs.txt
    pobs_file = os.path.join(path_to_save, 'Pobs.txt')
    pobs = ds['precipitation'].resample(time='1D').mean().to_dataframe()
    pobs_unstacked = pobs['precipitation'].unstack(level='COMID')
    pobs_unstacked.columns = [str(int(col)) for col in pobs_unstacked.columns]
    pobs_unstacked.to_csv(pobs_file, float_format="%.3f", sep='\t')

    # Tobs.txt
    tobs_file = os.path.join(path_to_save, 'Tobs.txt')
    tobs = ds['temperature'].resample(time='1D').mean().to_dataframe()
    tobs_unstacked = tobs['temperature'].unstack(level='COMID')
    tobs_unstacked.columns = [str(int(col)) for col in tobs_unstacked.columns]
    tobs_unstacked.to_csv(tobs_file, float_format="%.3f", sep='\t')

    # TMAXobs.txt
    tmaxobs_file = os.path.join(path_to_save, 'TMAXobs.txt')
    tmaxobs = ds['temperature'].resample(time='1D').max().to_dataframe()
    tmaxobs_unstacked = tmaxobs['temperature'].unstack(level='COMID')
    tmaxobs_unstacked.columns = [str(int(col)) for col in tmaxobs_unstacked.columns]
    tmaxobs_unstacked.to_csv(tmaxobs_file, float_format="%.3f", sep='\t')

    # TMINobs.txt
    tminobs_file = os.path.join(path_to_save, 'TMINobs.txt')
    tminobs = ds['temperature'].resample(time='1D').min().to_dataframe()
    tminobs_unstacked = tminobs['temperature'].unstack(level='COMID')
    tminobs_unstacked.columns = [str(int(col)) for col in tminobs_unstacked.columns]
    tminobs_unstacked.to_csv(tminobs_file, float_format="%.3f", sep='\t')

    return


# write GeoData and GeoClass files
def write_hype_geo_files(gistool_outputs, subbasins_shapefile, rivers_shapefile, frac_threshold, geofabric_mapping, path_to_save):
    
    if not os.path.isdir(path_to_save):
        os.makedirs(path_to_save)
    
    # extract geofabric mapping values
    basinID = geofabric_mapping['basinID']['in_varname']
    NextDownID = geofabric_mapping['nextDownID']['in_varname']

    # load the information from the gistool for soil and land cover and find the number of geoclass
    soil_type = pd.read_csv(gistool_outputs['soil'])
    landcover_type = pd.read_csv(gistool_outputs['landcover'])
    elevation_mean = pd.read_csv(gistool_outputs['elevation'])

    soil_type = soil_type.sort_values(by=basinID).reset_index(drop=True)
    landcover_type = landcover_type.sort_values(by=basinID).reset_index(drop=True)
    elevation_mean = elevation_mean.sort_values(by=basinID).reset_index(drop=True)

    # find the combination of the majority soil and land cover
    combinations_set_all = set()
    for index, row in landcover_type.iterrows():
        # get the fraction for land cover for each row
        fractions = [col for col in landcover_type.columns if col.startswith('frac') and row[col] > frac_threshold]
        # remove frac_ from the list
        fractions = [col.split('_')[1] for col in fractions]
        fractions = [int(name) for name in fractions]

        # get the majority soil type for each row
        majority_soil = [soil_type['majority'].iloc[index].item()]

        # Combine as combination of soil and land cover and keep as a set
        combinations = list(product(fractions, majority_soil))
        combinations_set = set(combinations)
        combinations_set_all.update(combinations_set)

    data_list = [{'landcover': item[0], 'soil': item[1]} for item in combinations_set_all]

    # Create a pandas DataFrame from the list of dictionaries
    combination = pd.DataFrame(data_list)

    combination ['SLC'] = 0
    combination ['SLC'] = np.arange(len(combination))+1
    
    #######################
    landcover_type_prepared = landcover_type.copy()

    for i in range(1, len(combination)+1):
        column_name = f'SLC_{i}'
        landcover_type_prepared[column_name] = 0.00

    landcover_type_prepared['soil'] = soil_type['majority']

    def get_non_zero_columns(row):
        return [col for col in row.index if col.startswith('frac_') and row[col] > frac_threshold]

    # Apply the function to each row
    landcover_type_prepared['non_zero_columns'] = landcover_type_prepared.apply(get_non_zero_columns, axis=1)


    for index, row in landcover_type_prepared.iterrows():
        # get the soil type
        soil_type_value = soil_type['majority'].iloc[index]

        for i in row['non_zero_columns']:

            # remove frac from column name 
            land_cover_value = i.replace("frac_", "")

            # get the SLC value
            result = combination[(combination['landcover'] == int(land_cover_value)) & (combination['soil'] == int(soil_type_value))]['SLC']
            column_name = 'SLC_'+str(result.values[0])
            landcover_type_prepared.loc[index, column_name] = landcover_type_prepared[i].iloc[index]
#######################
    riv = gpd.read_file(rivers_shapefile)
    riv.sort_values(by=basinID).reset_index(drop=True)
    riv['lengthm'] = 0.00
    
    rivlen_name = geofabric_mapping['rivlen']['in_varname']
    length_in_units = geofabric_mapping['rivlen']['in_units']
    length_out_units = geofabric_mapping['rivlen']['out_units']

    # Initialize a unit registry
    ureg = pint.UnitRegistry()
    # Quantify the DataFrame column with the original units
    lengthm = riv[rivlen_name].values * ureg(length_in_units)
    # Convert to the desired units
    riv['lengthm'] = lengthm.to(length_out_units).magnitude

    cat = gpd.read_file(subbasins_shapefile)
    cat.sort_values(by=basinID).reset_index(drop=True)
    cat['area'] = 0.00
    # cat['area'] = cat['unitarea'] * 1000000 # km2 to m2
    area_name = geofabric_mapping['area']['in_varname']
    area_in_units = geofabric_mapping['area']['in_units']
    area_out_units = geofabric_mapping['area']['out_units']

    # Initialize a unit registry
    ureg = pint.UnitRegistry()
    # Quantify the DataFrame column with the original units
    area_m2 = cat[area_name].values * ureg(area_in_units)
    # Convert to the desired units
    cat['area'] = area_m2.to(area_out_units).magnitude

    cat['latitude'] = cat.centroid.y
    cat['longitude'] = cat.centroid.x
    
    # add information to the geodata dataframe
    
    landcover_type_prepared[NextDownID] = riv[NextDownID]
    landcover_type_prepared['area'] = cat['area']
    landcover_type_prepared['latitude'] = cat['latitude']
    landcover_type_prepared['longitude'] = cat['longitude']
    landcover_type_prepared['elev_mean'] = elevation_mean['mean']
    landcover_type_prepared['slope_mean'] = riv['slope']
    landcover_type_prepared['rivlen'] = riv['lengthm']
    # landcover_type_prepared['uparea'] = riv['uparea']
    
    column_name_mapping = {
    basinID: 'subid',
    NextDownID: 'maindown',
    'area': 'area',
    'latitude': 'latitude',
    'longitude': 'longitude',
    'elev_mean': 'elev_mean',
    'slope_mean': 'slope_mean',
    'rivlen': 'rivlen'
    }

    # Rename the columns based on the dictionary
    landcover_type_prepared = landcover_type_prepared.rename(columns=column_name_mapping)

    slc_columns = [col for col in landcover_type_prepared.columns if col.startswith('SLC_')]

    # Sort the columns as per your requirements
    column_order = ['subid', 'maindown', 'area', 'latitude', 'longitude', 'elev_mean', 'slope_mean', 'rivlen'] + slc_columns

    landcover_type_prepared = landcover_type_prepared[column_order]
    #######################
    # sort geodata file from upstream to downstream
    landcover_type_prepared = sort_geodata(landcover_type_prepared)
    
    # normalize fracs
    # Identify columns starting with 'SLC_'
    slc_columns = [col for col in landcover_type_prepared.columns if col.startswith('SLC_')]

    # Normalize SLC values so that they sum to 1 for each row
    landcover_type_prepared[slc_columns] = landcover_type_prepared[slc_columns].div(landcover_type_prepared[slc_columns].sum(axis=1), axis=0)


    landcover_type_prepared.to_csv(path_to_save+'GeoData.txt', sep='\t', index=False)
    #######################
    
    # write geoclass file

    combination = combination.rename(columns={'landcover': 'LULC'})
    combination = combination.rename(columns={'soil': 'SOIL TYPE'})
    combination = combination[['SLC','LULC','SOIL TYPE']]
    combination['Main crop cropid'] = 0
    combination['Second crop cropid'] = 0
    combination['Crop rotation group'] = 0
    combination['Vegetation type'] = 1
    combination['Special class code'] = 0
    combination['Tile depth'] = 0
    combination['Stream depth'] = 2.296
    combination['Number of soil layers'] = 3
    combination['Soil layer depth 1'] = 0.091
    combination['Soil layer depth 2'] = 0.493
    combination['Soil layer depth 3'] = 2.296

    # combination
    #######################
    # Add commented lines
    commented_lines = [
    """! MODIS landcover													
! Add legend (raster value) and discription													
!	original legend (raster_value)	description											
!   1: 'Temperate or sub-polar needleleaf forest',
!   2: 'Sub-polar taiga needleleaf forest',
!   3: 'Tropical or sub-tropical broadleaf evergreen forest',
!   4: 'Tropical or sub-tropical broadleaf deciduous forest',
!   5: 'Temperate or sub-polar broadleaf deciduous forest',
!   6: 'Mixed forest',
!   7: 'Tropical or sub-tropical shrubland',
!   8: 'Temperate or sub-polar shrubland',
!   9: 'Tropical or sub-tropical grassland',
!   10: 'Temperate or sub-polar grassland',
!   11: 'Sub-polar or polar shrubland-lichen-moss',
!   12: 'Sub-polar or polar grassland-lichen-moss',
!   13: 'Sub-polar or polar barren-lichen-moss',
!   14: 'Wetland',
!   15: 'Cropland',
!   16: 'Barren lands',
!   17: 'Urban',
!   18: 'Water',
!   19: 'Snow and Ice',											
!													
!													
!													
! ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------													
!	SoilGrid V1												
!		original legend (raster_value)	description										
!	 C 	    1	 clay										
!	 SIC 	2	 silty clay										
!	 SC 	3	 sandy clay										
!	 CL 	4	 clay loam										
!	 SICL 	5	 silty clay loam										
!	 SCL 	6	 sandy clay loam										
!	 L   	7	 loam										
!	 SIL 	8	 silty loam										
!	 SL 	9	 sandy loam										
!	 SI 	10	 silt										
!	 LS 	11	 loamy sand										
!	 S  	12	 sand										
!													
!													
! ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	"""
    ]

    # Open the file in write mode
    with open(os.path.join(path_to_save, 'GeoClass.txt'), 'w') as file:
        # Write the commented lines
        for line in commented_lines:
            file.write(line + '\n')

    # re-number landcover and soil   
    for i in ['LULC', 'SOIL TYPE']:
        # Identify unique values and create a mapping from old values to new sequential values
        unique_values = pd.unique(combination[i])
        value_mapping = {value: idx + 1 for idx, value in enumerate(unique_values)}

        # Apply the new numbering to the LULC column
        combination[i] = combination[i].map(value_mapping)

        # Create a DataFrame for the mapping and write to geoclass
        mapping_df = pd.DataFrame(list(value_mapping.items()), columns=['Old Value', 'New Value'])

        with open(path_to_save+'GeoClass.txt', 'a') as f:
            f.write('! changes (reclassification) to '+i+'\n')
            for _, row in mapping_df.iterrows():
                # Write each row to the file with "!" at the beginning
                f.write(f"! {row['Old Value']} -> {row['New Value']}\n")


    # writing the `GeoClass.txt` file
    with open(os.path.join(path_to_save, 'GeoClass.txt'), 'a') as file:
            file.write("""! ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	
!          SLC	LULC	SOIL TYPE	Main crop cropid	Second crop cropid	Crop rotation group	Vegetation type	Special class code	Tile depth	Stream depth	Number of soil layers	Soil layer depth 1	Soil layer depth 2	Soil layer depth 3 \n""")
            combination.to_csv(file, sep='\t', index=False, header=False)
            

################################################################

# write par.txt file
def write_hype_par_file(path_to_save):

    output_file = os.path.join(path_to_save, 'par.txt')

    if os.path.isfile(output_file):
        os.remove(output_file)
    par_file = """!!	=======================================================================================================									
!! Parameter file for:										
!! HYPE -- Generated by the Model Agnostic Framework (hypeflow)									
!!	=======================================================================================================									
!!										
!!	------------------------									
!!										
!!	=======================================================================================================									
!!	"SNOW - MELT, ACCUMULATION, AND DISTRIBUTION; sublimation is sorted under Evapotranspiration"									
!!	-----									
!!	"General snow accumulation and melt related parameters (baseline values from SHYPE, unless noted otherwise)"									
ttpi	1.7083	!! width of the temperature interval with mixed precipitation								
sdnsnew	0.13	!! density of fresh snow (kg/dm3)								
snowdensdt	0.0016	!! snow densification parameter								
fsceff	1	!! efficiency of fractional snow cover to reduce melt and evap								
cmrefr	0.2	"!! snow refreeze capacity (fraction of degreeday melt factor) - baseline value from HBV (pers comm Barbro Johansson, but also in publications)"								
!!	-----									
!!	Landuse dependent snow melt parameters									
!!LUSE:	LU1	LU2	LU3	LU4	LU5					
ttmp	 -0.9253	 -1.5960	 -0.9620	 -2.7121	  2.6945    -0.9253	 -1.5960	 -0.9620	 -2.7121	  2.6945    -0.9253	 -1.5960	 -0.9620	 -2.7121	  2.6945    -0.9253	 -1.5960	 -0.9620	 -2.7121	  2.6945    !! Snowmelt threshold temperature (deg), baseline zero for all landuses"				
cmlt	   9.6497	   9.2928	   9.8897	   5.5393	   2.5333   9.6497	   9.2928	   9.8897	   5.5393	   2.5333   9.6497	   9.2928	   9.8897	   5.5393	   2.5333   9.6497	   9.2928	   9.8897	   5.5393	   2.5333	!! Snowmelt degree day coef (mm/deg/timestep)							
!!	-----									
!!	=======================================================================================================									
!!	EVAPOTRANSPIRATION PARAMETERS									
!!	-----									
!!	General evapotranspiration parameters									
lp	    0.6613	!! Threshold for water content reduction of transpiration (fraction of field capacity) - baseline value from SHYPE because its more realistic with a value slightly below field capacity								
epotdist	   4.7088	!! Coefficient in exponential function for potential evapotranspiration's depth dependency - baseline from EHYPE and/or SHYPE (very similar)																					
!!	-----									
!!										
!!LUSE:	LU1	LU2	LU3	LU4	LU5					
cevp	  0.4689	  0.7925	  0.6317	  0.1699	  0.4506    0.4689	  0.7925	  0.6317	  0.1699	  0.4506    0.4689	  0.7925	  0.6317	  0.1699	  0.4506    0.4689	  0.7925	  0.6317	  0.1699	  0.4506
ttrig	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	!! Soil temperature threshold to allow transpiration - disabled if treda is set to zero				
treda	0.84	0.84	0.84	0.84	0.95	0.84	0.84	0.84	0.84	0.95	0.84	0.84	0.84	0.84	0.95	0.84	0.84	0.84	0.84	0.95	"!! Coefficient in soil temperature response function for root water uptake, default value from �gren et al, set to zero to disable the function"				
tredb	0.4	0.4	0.4	0.4	0.4	0.4	0.4	0.4	0.4	0.4	0.4	0.4	0.4	0.4	0.4	0.4	0.4	0.4	0.4	0.4	"!! Coefficient in soil temperature response fuction for root water uptake, default value from �gren et al"				
fepotsnow	0.8	0.8	0.8	0.8	0.8	0.8	0.8	0.8	0.8	0.8	0.8	0.8	0.8	0.8	0.8	0.8	0.8	0.8	0.8	0.8	!! Fraction of potential evapotranspiration used for snow sublimation				
!!										
!! Frozen soil infiltration parameters										
!! SOIL:	S1	S2								
bfroznsoil  3.7518  3.2838  3.7518  3.2838  3.7518  3.2838  3.7518  3.2838  3.7518  3.2838  3.7518  3.2838  3.7518  3.2838  3.7518  3.2838  3.7518  3.2838  3.7518  3.2838								
logsatmp	1.15	1.15	1.15	1.15	1.15	1.15	1.15	1.15	1.15	1.15	1.15	1.15	1.15	1.15	1.15	1.15	1.15	1.15	1.15	1.15								
bcosby	    11.2208	    19.6669	    11.2208	    19.6669	    11.2208	    19.6669	    11.2208	    19.6669	    11.2208	    19.6669	    11.2208	    19.6669	    11.2208	    19.6669	    11.2208	    19.6669	    11.2208	    19.6669	    11.2208	    19.6669								
!!	=======================================================================================================									
!!	"SOIL/LAND HYDRAULIC RESPONSE PARAMETERS - recession coef., water retention, infiltration, macropore, surface runoff; etc."									
!!	-----									
!!	Soil-class parameters									
!!	S1	S2								
rrcs1   0.4345  0.5985   0.4345  0.5985   0.4345  0.5985   0.4345  0.5985   0.4345  0.5985   0.4345  0.5985   0.4345  0.5985   0.4345  0.5985   0.4345  0.5985   0.4345  0.5985	!! recession coefficients uppermost layer (fraction of water content above field capacity/timestep)							
rrcs2   0.1201  0.1853   0.1201  0.1853   0.1201  0.1853   0.1201  0.1853   0.1201  0.1853   0.1201  0.1853   0.1201  0.1853   0.1201  0.1853   0.1201  0.1853   0.1201  0.1853	!! recession coefficients bottom layer (fraction of water content above field capacity/timestep)							
rrcs3	    0.0939	!! Recession coefficient (upper layer) slope dependance (fraction/deg)								
sfrost  1   1  1   1  1   1  1   1  1   1  1   1  1   1  1   1  1   1  1   1	!! frost depth parameter (cm/degree Celsius) soil-type dependent							
wcwp    0.1171  0.0280    0.1171  0.0280    0.1171  0.0280    0.1171  0.0280    0.1171  0.0280    0.1171  0.0280    0.1171  0.0280    0.1171  0.0280    0.1171  0.0280    0.1171  0.0280	!! Soil water content at wilting point (volume fraction)											
wcfc    0.3771  0.2009    0.3771  0.2009    0.3771  0.2009    0.3771  0.2009    0.3771  0.2009    0.3771  0.2009    0.3771  0.2009    0.3771  0.2009    0.3771  0.2009    0.3771  0.2009	!! Field capacity, layerOne (additional to wilting point) (volume fraction)"										
wcep    0.4047  0.4165    0.4047  0.4165    0.4047  0.4165    0.4047  0.4165    0.4047  0.4165    0.4047  0.4165    0.4047  0.4165    0.4047  0.4165    0.4047  0.4165    0.4047  0.4165	!! Effective porosity, layerOne (additional to wp and fc) (volume fraction)"							
!!	-----									
!!	Landuse-class parameters	parameters								
!!LUSE:	LU1	LU2	LU3	LU4	LU5					
srrcs   0.0673  0.1012  0.1984  0.0202  0.0202   0.0673  0.1012  0.1984  0.0202  0.0202   0.0673  0.1012  0.1984  0.0202  0.0202   0.0673  0.1012  0.1984  0.0202  0.0202	!! Runoff coefficient for surface runoff from saturated overland flow of uppermost soil layer (fraction/timestep)				
!!	-----									
!!	Regional groundwater outflow									
rcgrw	0	!! recession coefficient for regional groundwater outflow from soil layers								
!!	=======================================================================================================									
!!	SOIL TEMPERATURE AND SOIL FROST DEPT									
!!	-----									
!!	General									
deepmem	1000	!! temperature memory of deep soil (days)								!! temperature memory of deep soil (days)							
!!-----										
!!LUSE:	LU1	LU2	LU3	LU4	LU5					
surfmem 17.8	17.8	17.8	17.8	5.15 17.8	17.8	17.8	17.8	5.15 17.8	17.8	17.8	17.8	5.15 17.8	17.8	17.8	17.8	5.15	!! upper soil layer soil temperature memory (days)				
depthrel	1.1152	1.1152	1.1152	1.1152	2.47    1.1152	1.1152	1.1152	1.1152	2.47	1.1152	1.1152	1.1152	1.1152	2.47	1.1152	1.1152	1.1152	1.1152	2.47	!! depth relation for soil temperature memory (/m)				
frost	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	!! frost depth parameter (cm/degree Celsius) soil-type dependent				
!!	-----									
!!	=======================================================================================================									
!!	LAKE DISCHARGE									
!!	-----									
!!	-----									
!!	"ILAKE and OLAKE REGIONAL PARAMETERS (1 ilakeregions , defined in geodata)"									
!!	ILAKE parameters																	
!! ilRegion	PPR 1									
ilratk  149.9593						
ilratp  4.9537						
illdepth    0.33					
ilicatch    1.0								
!!										
!!	=======================================================================================================									
!!	RIVER ROUTING									
!!	-----									
damp	   0.2719	!! fraction of delay in the watercourse which also causes damping								
rivvel	     9.7605	!! celerity of flood in watercourse (rivvel>0)								
qmean 	200	!! initial value for calculation of mean flow (mm/yr) - can also be given in LakeData								"""

    with open(output_file, 'w') as file:
            file.write(par_file)
    
################################################################
# write info and filedir files
def write_hype_info_filedir_files(path_to_save, spinup_days):
    # write filedir file
    output_file = os.path.join(path_to_save, 'filedir.txt')

    if os.path.isfile(output_file):
        os.remove(output_file)

    with open(output_file, 'w') as file:
            file.write('./')
    # create results directory

    if not os.path.isdir(os.path.join(path_to_save, 'results')):
        os.makedirs(os.path.join(path_to_save, 'results'))

    ###########

    # Output par to a .txt file

    output_file = os.path.join(path_to_save, 'info.txt')
    if os.path.isfile(output_file):
        os.remove(output_file)


    # define start time, end time, based on input forcing
    # spinup period is a user defined inputs

    Pobs = pd.read_csv(os.path.join(path_to_save, 'Pobs.txt'), sep='\t', parse_dates=['time'])
    Pobs['time'] = Pobs['time'].dt.date
    start_date = Pobs['time'].iloc[0]
    end_date = Pobs['time'].iloc[-1]
    spinup_date = start_date + pd.Timedelta(days=spinup_days)


    # write out first text section
    s1= [
    """!! ----------------------------------------------------------------------------							
!!							
!! HYPE - Model Agnostic Framework
!!							
!! -----------------------------------------------------------------------------							
!! Check Indata during first runs (deactivate after first runs) 
indatacheckonoff 	2						
indatachecklevel	2		
!! -----------------------------------------------------------------------------							"""
    ]

    # write s1 in output file
    with open(output_file, 'w') as file:
        # Write the commented lines
        for line in s1:
            file.write(line + '\n')

    # write out first text section
    s2= [
    """!!
!! -----------------------------------------------------------------------------							
!!						
!! Simulation settings:							
!!							
!! -----------------	 """
    ]

    # write s2 in output file
    with open(output_file, 'a') as file:
        # Write the commented lines
        for line in s2:
            file.write(line + '\n')

    # create df2
    df2_row=['bdate', 'cdate', 'edate', 'resultdir', 'instate', 'warning']
    df2_val=[start_date, spinup_date, end_date, './results/', 'n', 'y']
    df2=pd.DataFrame(df2_val, index=df2_row, columns=None)

    # append df2
    with open(output_file, 'a') as file:
        # Write the DataFrame to the file
        df2.to_csv(file, sep='\t', index=True, header=False)#, line_terminator='\n')

    # write out s3
    s3= [
    """readdaily 	y						
submodel 	n						
calibration	n						
readobsid   n							
soilstretch	n						
!! Soilstretch enable the use of soilcorr parameters (strech soildepths in layer 2 and 3)
steplength	1d							
!! -----------------------------------------------------------------------------							
!!							
!! Enable/disable optional input files
!!							
!! -----------------							
readsfobs	n	!! For observed snowfall fractions in SFobs.txt							
readswobs	n	!! For observed shortwave radiation in SWobs.txt
readuobs	n	!! For observed wind speeds in Uobs.txt
readrhobs	n	!! For observed relative humidity in RHobs.txt					
readtminobs	y	!! For observed min air temperature in TMINobs.txt				
readtmaxobs	y	!! For observed max air temperature in TMAXobs.txt
soiliniwet	n	!! initiates soil water to porosity instead of field capacity which is default (N). Set Y to use porosity.
usestop84	n	!! flag to use the old return code 84 for a successful run					
!! -----------------------------------------------------------------------------							
!!							
!! Define model options (optional)
!!							
!! -----------------							
!!snowfallmodel:								
!!                  0 threshold temperature model							
!!                  1 inputdata (SFobs.txt)							
!!snowmeltmodel:							
!!                  0,1 temperature index             (with/without snowcover scaling)							
!!                  2   temperature + radiation index (with/without snowcover scaling)							
!!							
!!  snowevapmodel   0 off							
!!                  1 on							
!!                   							
!!  petmodel:  (potential evapotranspiration) (is shown in geodata for WWH)							
!!                  0 original HYPE temperature model (with Xobs epot replacement)							
!!                  1 original HYPE temperature model (without Xobs epot replacement)							
!!                  2 Modified Jensen-Haise 							
!!                  3 Modified Hargreaves-Samani							
!!                  4 Priestly-Taylor							
!!                  5 FAo Penman-Monteith							
!!							
!! lakeriverice:							
!!                  0 off							
!!                  1 on, old (simple) air-water heat exchange              (requires T2 water temperature model)							
!!                  2 on, new heatbalance model for air-water heat exchange (requires T2 water temperature model)							
!!							
!! substance T2     switching on the new water temperature trace model							
!!							
!! deepground       0   off    Deep groundwater (Aquifer) model options							
!!                  1,2 on
!! Glacierini	0 off 1 on	(1 used for statefile preparation)	
!! Floodplain		0, 1, 2, 3 (3 used for WWH)					
!! -----------------							
modeloption snowfallmodel	0						
modeloption snowdensity	0
modeloption snowfalldist	2
modeloption snowheat	0
modeloption snowmeltmodel	0	
modeloption	snowevapmodel	1				
modeloption snowevaporation	1					
modeloption lakeriverice	0									
modeloption deepground	0 	
modeloption glacierini	1
modeloption floodmodel 0
modeloption frozensoil 2
modeloption infiltration 3
modeloption surfacerunoff 0
modeloption petmodel	1
modeloption wetlandmodel 2		
modeloption connectivity 0					
!! ------------------------------------------------------------------------------------							
!!							
!! Define outputs
!!							
!! -----------------							
!! meanperiod 1=daymean, 2=weekmean, 3=monthmean, 4=yearmean, 5=total period mean							
!! output variables: see http://www.smhi.net/hype/wiki/doku.php?id=start:hype_file_reference:info.txt:variables 
!! -----------------							
!! BASIN outputs 
!! The present basins are some large rivers distributed over different continents
!! -----------------
!! basinoutput variable	rout	cout	cilv	evap	fnca	fcon
!! basinoutput meanperiod	1						
!! basinoutput decimals	3
!! basinoutput subbasin	basinID1    basinID2    basinID3		
!! -----------------							
!! TIME outputs 
!! -----------------	
timeoutput variable cout	evap	snow
timeoutput meanperiod	1
timeoutput decimals	3					
!! -----------------							
!! MAP outputs
!! -----------------							
!! mapoutput variable	cout cprc ctmp
!! mapoutput decimals	3						
!! mapoutput meanperiod	5						
!! ------------------------------------------------------------------------------------							
!!							
!! Select criteria for model evaluation and automatic calibration
!!							
!! -----------------							
!! General settings
!! -----------------			
!! crit meanperiod	1
!! crit datalimit	30
!! crit subbasin	outletBasinID
!! -----------------			
!! Criterion-specific settings
!! -----------------				
!! crit 1 criterion	MKG
!! crit 1 cvariable	cout
!! crit 1 rvariable	rout
!! crit 1 weight	1"""
    ]

    # write s3 in output file
    with open(output_file, 'a') as file:
        # Write the commented lines
        for line in s3:
            file.write(line + '\n')