# GetMap
A simple Flask app to clip and ship a variety of grazing land management products in response to a GET request.
Relies on seasonal groundcover composites.

Call using:

  * __Site__ is a label for the site. Will go into the final file name
  * __Strata__ is the field in the polygon attributes that is unique for each paddock if doing paddock based analysis
  * __startDate__ is the starting date of the analysis
  * __endDate__ is the finishing date of the analysis
  * __compareDate__ is the date used for comparison in the mean and rank anomaly products
  * __jsonPoly__ is the polygon of the property or paddocks in geojson format
  * __atype__ is the analysis type. One of _mean, rank, diff, paddock_ or _ property_. See below for details
  * __seasonal__ is a _True_ or _False_ flag indicating whether seasonal or whole of timeseries processing is required.

## Analysis Types
### Mean

### Rank

### Diff

### Property

### Paddock


## Deployment


##
