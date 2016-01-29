
import pysal, shapefile, networkx, time, cPickle, random, math, copy, Convexpath_module
from shapely.geometry import Point, Polygon, LineString, MultiPoint, MultiPolygon
from collections import defaultdict
from shapely.ops import cascaded_union

path = "F:\\Dropbox\\research\\Convexpath Approach\\Convexpath buffer\\data\\"

facil = "pt1.shp"
obstacles_f = "ASU_part1"
coverage = 500

def generateGeometry(in_shp):
    resultingGeometry = []
    if in_shp.header['Shape Type'] == 1:
        for i in range(len(in_shp)):
            resultingGeometry.append(Point(in_shp.get_shape(i)['X'], in_shp.get_shape(i)['Y']))
    elif in_shp.header['Shape Type'] == 3:
        for i in range(len(in_shp)):
            resultingGeometry.append(LineString(in_shp.get_shape(i)['Vertices']))
    elif in_shp.header['Shape Type'] == 5:
        for i in range(len(in_shp)):
            resultingGeometry.append(Polygon(in_shp.get_shape(i)['Vertices']))
    return resultingGeometry  

def createConvexhull(poly, endPoints = []):
    convexVerticex = []
    if poly.type == 'MultiPolygon':
        for i in poly.geoms:
            convexVerticex.extend(list(i.exterior.coords))
    else:
        convexVerticex.extend(list(poly.exterior.coords))
    convexVerticex.extend(endPoints)
    convex_hull = MultiPoint(convexVerticex).convex_hull
    
    return convex_hull

def split_circle(center, poly, circle):
    #center = [(x,y)]
    ch = createConvexhull(poly, center)
    minX, minY, maxX, maxY = circle.bounds
    vertices = list(ch.exterior.coords)
    vertices = list(set(vertices))

    sector_lines = []
    for i in vertices:    
        a = LineString([center[0], i])
        if a.crosses(poly) == False:
            if a.length != 0:
                qurt = 0
                xx = center[0][0] - i[0]
                yy = center[0][1] - i[1]
                if xx * yy > 0:
                    if xx < 0:
                        qurt = 1
                    else:
                        qurt = 3
                else:
                    if xx > 0:
                        qurt = 2
                    else:
                        qurt = 4
                sector_lines.append((a, qurt, i))    
    split_lines = []
    for line in sector_lines: 
        seg_coords = list(line[0].coords)
        slope = (seg_coords[0][1]-seg_coords[1][1])/(seg_coords[0][0]-seg_coords[1][0])
        b = seg_coords[0][1] - slope * seg_coords[0][0]
                
        if line[1] == 1:
            if slope < 1:
                c = (maxX, slope * maxX + b)
            elif slope > 1:
                c = ((maxY-b)/slope, maxY)
        elif line[1] == 2:
            if slope < -1:
                c = ((maxY-b)/slope, maxY)
            elif slope > -1:
                c = (minX, slope * minX + b)
        elif line[1] == 3:
            if slope < 1:
                c = (minX, slope * minX + b)
            elif slope > 1:
                c = ((minY-b)/slope, minY)
        elif line[1] == 4:
            if slope < -1:
                c = ((minY-b)/slope, minY)
            elif slope > -1:
                c = (maxX, slope * maxX + b)            
                
        split_lines.append(LineString([center[0], c]))
    qurts = [x[1] for x in sector_lines]
    
    centroid = poly.centroid 
    
    
    return split_lines
    
        

facil_pysal = pysal.IOHandlers.pyShpIO.shp_file(path+facil)
obstacles_pysal = pysal.IOHandlers.pyShpIO.shp_file(path + obstacles_f)
obstacles_shp = generateGeometry(obstacles_pysal)
facil_shp = generateGeometry(facil_pysal)
facil_coords = list(facil_shp[0].coords)


circle_cover = facil_shp[0].buffer(coverage)
c_obs = [x for x in obstacles_shp if circle_cover.intersects(x)]
clip_obs = [circle_cover.intersection(x) for x in c_obs]
obs_covered = []
cent_c_obs = [(x, x.centroid) for x in c_obs]
dis_cent = [(x[1].distance(facil_shp[0]), x[0]) for x in cent_c_obs]
dis_cent.sort()

ch = createConvexhull(dis_cent[0][1], facil_coords )
chs = []
for i in dis_cent:
    chs.append(createConvexhull(i[1], facil_coords))

vertices = list(ch.exterior.coords)

splits = split_circle(facil_coords, dis_cent[0][1], circle_cover)

#needs to deal with dulplication in sector_lines         


w = shapefile.Writer(shapefile.POLYLINE)
w.field('nem')
for line in splits:
    w.line(parts=[[ list(x) for x in list(line.coords) ]])
    w.record('ff')
w.save(path + "lines")

w = shapefile.Writer(shapefile.POLYGON)
w.field('net')
for obs in chs:
    w.poly(parts=[[list(x) for x in list(obs.exterior.coords)]])
    w.record('ff')
w.save(path + "obs") 

    
w = shapefile.Writer(shapefile.POLYGON)
w.field('net')
for obs in clip_obs:
    if obs.type == 'MultiPolygon':
        for i in obs.geoms:
            w.poly(parts=[[list(x) for x in list(i.exterior.coords)]])
            w.record('ff')
    else:
        w.poly(parts=[[list(x) for x in list(obs.exterior.coords)]])
        w.record('ff')
w.save(path + "clip_obs") 

