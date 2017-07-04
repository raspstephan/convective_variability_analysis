import pygrib

### returns grib data, given a prescribed absolute path and a 'field' 
def grbdat(path, field):
    grib_load   = path
    grib_file   = pygrib.open(grib_load)
    grib        = grib_file.select(name=field)[0]
    data        = grib.values
    grib_file.close()
    return data  

### returns lat, lon
def latlon(path):
    grib_load   = path
    grib_file   = pygrib.open(grib_load)
    grib        = grib_file.select()[0]
    lat, lon    = grib.latlons()
    return lat, lon
