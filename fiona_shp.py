from shapely.geometry import Point, Polygon, LineString, MultiPoint, MultiPolygon, mapping
import fiona

def openShp(in_f):

    with fiona.open(path + in_f, 'r') as a:
        sh = []
        if list(a)[0]['geometry']['type'] == 'Polygon':
            for i in list(a):
                if len(i['geometry']['coordinates']) == 1:
                    sh.append(Polygon(i['geometry']['coordinates'][0]))
                else:
                    
                    sh.append(Polygon(i['geometry']['coordinates'][0],i['geometry']['coordinates'][1:]))
        elif list(a)[0]['geometry']['type'] == 'LineString':
            for i in list(a):
                sh.append(LineString(i['geometry']['coorindates']))       
        elif list(a)[0]['geometry']['type'] == 'Point':
            for i in list(a):
                sh.append(Point(i['geometry']['coordinates']))
        
        return sh

def writeShp(in_shp, out_name):
    """
    in_shp must be list 
    """
    if in_shp[0].type == 'Polygon':
        schema = schema_poly
    elif in_shp[0].type == 'LineString':
        schema = schema_line
    elif in_shp[0].type == 'Point':
        schema = schema_point
    else:
        raise TypeError
    with fiona.open(path + out_name, 'w', 'ESRI Shapefile', schema) as c:
        for i in in_shp:
            c.write({
                'geometry': mapping(i),
                'properties': {'id': in_shp.index(i)},
            })            