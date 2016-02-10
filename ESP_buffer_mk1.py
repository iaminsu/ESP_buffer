
import pysal, shapefile, networkx, time, cPickle, random, math, copy, Convexpath_module
from shapely.geometry import Point, Polygon, LineString, MultiPoint, MultiPolygon
from shapely.ops import cascaded_union
from arcpy import env, Dissolve_management, Union_analysis, Intersect_analysis, Erase_analysis, MultipartToSinglepart_management, Buffer_analysis
import networkx


env.overwriteOutput = True
path = "E:\\data\\esp_cover\\data2\\"
env.workspace = path

facil = "pt1.shp"
obstacles_f = "convex_obstacles.shp"
coverage_distance = 500
cover = "pt_buffer.shp"

clip_obs = "clip_obs.shp"
clip_obs_single = "clip_obs_single.shp"
def generateGeometry(in_shp):
    """
    Returns objects in list
    """
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

def saveShp (inShp, inType, inPath, shpName):
    if type(inShp) != type([]):
        inShp = list(inShp)
    if inType == 'line':
        w = shapefile.Writer(shapefile.POLYLINE)
        w.field('nem')
        for line in inShp:
            w.line(parts=[[ list(x) for x in list(line.coords) ]])
            w.record('ff')
        w.save(inPath + shpName)    
    elif inType == "polygon":
        w = shapefile.Writer(shapefile.POLYGON)
        w.field('net')
        for obs in inShp:
            w.poly(parts=[[list(x) for x in list(obs.exterior.coords)]])
            w.record('ff')
        w.save(inPath + shpName)         
    elif inType == 'point':
        w = shapefile.Writer(shapefile.POINT)
        w.field('net')
        for pt in inShp:
            w.point(pt.x, pt.y)
            w.record('ff')
        w.save(inPath + shpName)

class SplitLine:
    """
    Only for single segment LineString object 
    fromNode is tuple of coordinates
    Intersecting line is LineString Object 
    """
    def __init__(self, line, fromNode):
         
        self.fromNode= fromNode
        if type(line) is LineString:
            self.lineObj = line
        else:
            self.lineObj = LineString(line)
        self.qurt = 0
        self.a = 0 
        self.b = 0 
        self.octa = 'I'
        self.MBR = "N"
        self.cuttingLine = ''  #line for cutbox
        self.radiusLine = ''   #radius line for sector 
        coords = list(self.lineObj.coords)
        coords = [x for x in coords if Point(x) != Point(self.fromNode)]
        self.toNode = coords[0]
        xx = self.toNode[0] - self.fromNode[0]
        
        yy = self.toNode[1] - self.fromNode[1]
        
        self.a = yy/xx
        self.b = self.fromNode[1] - self.a * self.fromNode[0]
        #I: maxX; II: maxY; III: maxY; IV: minX; 
        #V: minX; VI: minY; VII: minY; VIII: maxX
                
        if xx * yy > 0:
            if xx > 0:
                self.qurt = 1
                if self.a < 1:
                    self.octa = "I"
                    self.MBR = "E"
                elif self.a > 1:
                    self.octa = "II"
                    self.MBR = "N"
            else:
                self.qurt = 3
                if self.a < 1: 
                    self.octa = "V"
                    self.MBR = "W"
                elif self.a > 1:
                    self.octa = "VI"
                    self.MBR = "S"
        else:
            if xx < 0: 
                self.qurt = 2
                if self.a < -1:
                    self.octa = "III"
                    self.MBR = "N"
                elif self.a > -1:
                    self.octa = "IV"
                    self.MBR = "W"
            else:
                self.qurt = 4
                if self.a < -1: 
                    self.octa = "VII"
                    self.MBR = "S"
                elif self.a > -1:
                    self.octa = "VIII"
                    self.MBR = "E"
                    
        
    def intersectPoint(self, intersectingLine):
        coords = list(intersectingLine.coords)
        xx = coords[1][0] - coords[0][0]
        yy = coords[1][1] - coords[0][1]
        if xx == 0:
            interX = coords[0][0]
            interY = self.a * interX + self.b
        elif yy == 0:
            inA = 0
            inB = coords[0][1]
            interX = (inB - self.b)/(self.a - inA)
            interY = inB
        else:
            inA = xx/yy
            inB = coords[0][1] - inA * coords[0][0] 
            interX = (inB - self.b)/(self.a - inA)
            interY = inA * interX + inB 
        self.cuttingLine = LineString([self.fromNode, (interX, interY)])
        return (interX, interY)
    
    def radius(self, inCircle):
        if self.cuttingLine == '':
            raise 
        else:
            circleBoundary = inCircle.boundary
            inter = self.cuttingLine.intersection(circleBoundary)
            self.radiusLine = LineString([self.fromNode, (inter.x, inter.y)])
        
    
    
        



def cutting_box (sLine1, sLine2, inCircle, circleCenter):
    minX, minY, maxX, maxY = inCircle.bounds 
    maxXLine = LineString([[maxX, maxY], [maxX, minY]])
    maxYLine = LineString([[minX, maxY], [maxX, maxY]])
    minXLine = LineString([[minX, maxY], [minX, minY]])
    minYLine = LineString([[minX, minY], [maxX, minY]])
    boundingBox = Polygon([(minX, minY), (minX, maxY), (maxX, maxY), (maxX, minY)])    
    if sLine1.octa == "I":
        sLin1_interPoint = sLine1.intersectPoint(maxXLine)
    elif sLine1.octa == "II":
        sLin1_interPoint = sLine1.intersectPoint(maxYLine)
    elif sLine1.octa == "III":
        sLin1_interPoint = sLine1.intersectPoint(maxYLine)
    elif sLine1.octa == "IV":
        sLin1_interPoint = sLine1.intersectPoint(minXLine)
    elif sLine1.octa == "V":
        sLin1_interPoint = sLine1.intersectPoint(minXLine)
    elif sLine1.octa == "VI":
        sLin1_interPoint = sLine1.intersectPoint(minYLine)
    elif sLine1.octa == "VII":
        sLin1_interPoint = sLine1.intersectPoint(minYLine)
    elif sLine1.octa == "VIII":
        sLin1_interPoint = sLine1.intersectPoint(maxXLine)
    
    if sLine2.octa == "I":
        sLine2_interPoint = sLine2.intersectPoint(maxXLine)
    elif sLine2.octa == "II":
        sLine2_interPoint = sLine2.intersectPoint(maxYLine)
    elif sLine2.octa == "III":
        sLine2_interPoint = sLine2.intersectPoint(maxYLine)
    elif sLine2.octa == "IV":
        sLine2_interPoint = sLine2.intersectPoint(minXLine)
    elif sLine2.octa == "V":
        sLine2_interPoint = sLine2.intersectPoint(minXLine)
    elif sLine2.octa == "VI":
        sLine2_interPoint = sLine2.intersectPoint(minYLine)
    elif sLine2.octa == "VII":
        sLine2_interPoint = sLine2.intersectPoint(minYLine)
    elif sLine2.octa == "VIII":
        sLine2_interPoint = sLine2.intersectPoint(maxXLine)
    
    #cutting box 
    #coordinates should be in order 
    #determine whether two split lines interset same MBR boundary segment 
    cutting_box_1 = []
    cutting_box_2 = []
    
    if sLine1.MBR == "N":        
        if sLine2.MBR == "E":
            cutting_box_1 = Polygon([circleCenter, sLin1_interPoint, (maxX, maxY), sLine2_interPoint])
            cutting_box_2 = boundingBox.difference(cutting_box_1)
        elif sLine2.MBR == "S":
            cutting_box_1 = Polygon([circleCenter, sLin1_interPoint, (maxX, maxY), (maxX, minY), sLine2_interPoint])
            cutting_box_2 = boundingBox.difference(cutting_box_1)
        elif sLine2.MBR == "W":
            cutting_box_1 = Polygon([circleCenter, sLin1_interPoint, (minX, maxY), sLine2_interPoint])
            cutting_box_2 = boundingBox.difference(cutting_box_1)
        elif sLine2.MBR == "N":
            cutting_box_1 = Polygon([circleCenter, sLin1_interPoint, sLine2_interPoint])
            cutting_box_2 = boundingBox.difference(cutting_box_1)
    if sLine1.MBR == "E":
        if sLine2.MBR == "S":
            cutting_box_1 = Polygon([circleCenter, sLin1_interPoint, (maxX, minY), sLine2_interPoint])
            cutting_box_2 = boundingBox.difference(cutting_box_1)
        elif sLine2.MBR == "W":
            cutting_box_1 = Polygon([circleCenter, sLin1_interPoint, (maxX, minY), (minX, minY), sLine2_interPoint])
            cutting_box_2 = boundingBox.difference(cutting_box_1)
        elif sLine2.MBR == "N":
            cutting_box_1 = Polygon([circleCenter, sLin1_interPoint, (maxX, maxY), sLine2_interPoint])
            cutting_box_2 = boundingBox.difference(cutting_box_1)
        elif sLine2.MBR == "E":
            cutting_box_1 = Polygon([circleCenter, sLin1_interPoint, sLine2_interPoint])
            cutting_box_2 = boundingBox.difference(cutting_box_1)
    if sLine1.MBR == "S":
        if sLine2.MBR == "W":
            cutting_box_1 = Polygon([circleCenter, sLin1_interPoint, (minX, minY), sLine2_interPoint])
            cutting_box_2 = boundingBox.difference(cutting_box_1)
        elif sLine2.MBR == "N":
            cutting_box_1 = Polygon([circleCenter, sLin1_interPoint, (minX, minY), (minX, maxY), sLine2_interPoint])
            cutting_box_2 = boundingBox.difference(cutting_box_1)
        elif sLine2.MBR == "E":
            cutting_box_1 = Polygon([circleCenter, sLin1_interPoint, (maxX, minY), sLine2_interPoint])
            cutting_box_2 = boundingBox.difference(cutting_box_1)
        elif sLine2.MBR == "S":
            cutting_box_1 = Polygon([circleCenter, sLin1_interPoint, sLine2_interPoint])
            cutting_box_2 = boundingBox.difference(cutting_box_1)
    if sLine1.MBR == "W":
        if sLine2.MBR == "N":
            cutting_box_1 = Polygon([circleCenter, sLin1_interPoint, (minX, maxY), sLine2_interPoint])
            cutting_box_2 = boundingBox.difference(cutting_box_1)
        elif sLine2.MBR == "E":
            cutting_box_1 = Polygon([circleCenter, sLin1_interPoint, (minX, maxY), (maxX, maxY), sLine2_interPoint])
            cutting_box_2 = boundingBox.difference(cutting_box_1)
        elif sLine2.MBR == "S":
            cutting_box_1 = Polygon([circleCenter, sLin1_interPoint, (minX, minY), sLine2_interPoint])
            cutting_box_2 = boundingBox.difference(cutting_box_1)
        elif sLine2.MBR == "W":
            cutting_box_1 = Polygon([circleCenter, sLin1_interPoint, sLine2_interPoint])
            cutting_box_2 = boundingBox.difference(cutting_box_1)        
    return cutting_box_1, cutting_box_2

#initializing objects 
def file2Shp(in_file):
    f_pysal = pysal.IOHandlers.pyShpIO.shp_file(path + in_file)
    f_shp = generateGeometry(f_pysal)
    return f_shp

def rrd(in_coords):
    if type(in_coords) == type(tuple()):
        out = (round(in_coords[0], 5), round(in_coords[1], 5))
        return out
    elif type(in_coords) == type(list()):
        out = [(round(x[0], 5), round(x[1], 5)) for x in in_coords]
        return out
    else:
        raise TypeError
    
def createGraph(inShp):
    results = networkx.Graph()
    lineList = [list(x.coords) for x in inShp]
    for line in lineList:
        results.add_edge(line[0], line[1], weight = LineString(list(line)).length)
    return results

obstacles_shp = file2Shp(obstacles_f)
facil_shp = file2Shp(facil)
facil_coords = list(facil_shp[0].coords)

#generate initial, circular coverage 
Buffer_analysis(facil, 'cover_shp.shp', coverage_distance)
#cover_shp = facil_shp[0].buffer(coverage_distance)
#saveShp([cover_shp], 'polygon', path, 'cover_shp.shp')
cover_shp = file2Shp('cover_shp.shp')[0]
IminX, IminY, ImaxX, ImaxY = cover_shp.bounds 
IminX = rrd(IminX)
IminY = rrd(IminY)
ImaxX = rrd(ImaxX)
ImaxY = rrd(ImaxY)
ImaxXLine = LineString([[ImaxX, ImaxY], [ImaxX, IminY]])
ImaxYLine = LineString([[IminX, ImaxY], [ImaxX, ImaxY]])
IminXLine = LineString([[IminX, ImaxY], [IminX, IminY]])
IminYLine = LineString([[IminX, IminY], [ImaxX, IminY]])
boundingBox = Polygon([(IminX, IminY), (IminX, ImaxY), (ImaxX, ImaxY), (ImaxX, IminY)])

#generate set of obstacles for circular cover 
Intersect_analysis(['cover_shp.shp', obstacles_f], 'initial_obs.shp')
clip_obs_shp = file2Shp('initial_obs.shp')
#clip_obs_shp = [x.intersection(cover_shp) for x in obstacles_shp]
#clip_obs_shp = [x for x in clip_obs_shp if not x.is_empty]
init_clip_obs_shp = file2Shp('initial_obs.shp')
#init_clip_obs_shp = [x for x in clip_obs_shp if not x.is_empty]

terminate = False
#begin main loop 
dealt_obs = []

while terminate == False:
    #first loop: generate 
    #For each iteration, list of obstacles is updated using updated coverage. 
    #Remove obstacles that cannot be covered with the updated coverage 
    clip_obs_shp = [x for x in clip_obs_shp if not x in dealt_obs]   
    saveShp(clip_obs_shp, 'polygon', path, 'clip_obs.shp')
    Intersect_analysis(['cover_shp.shp', 'clip_obs.shp'], 'clip_obs_inter.shp')
    clip_obs_multi = file2Shp('clip_obs_inter.shp')
    MultipartToSinglepart_management('clip_obs_inter.shp', 'clip_obs.shp')
    clip_obs_shp = file2Shp('clip_obs.shp')
    #if len(clip_obs_shp) == 0:
        #break
    
    #clip_obs_shp = [x.intersection(cover_shp) for x in clip_obs_shp]
    #clip_obs_shp = [x for x in clip_obs_shp if not x.is_empty]    
    
    #convex hulls for clipped obstacles
    #sort them as distance between facility and centorids 
    cent_clip_obs = [(x.centroid, x) for x in clip_obs_shp]
    cent_dist_clip_obs = [(x[0].distance(facil_shp[0]), x[1]) for x in cent_clip_obs]
    cent_dist_clip_obs.sort()
    dealt_obs.append(cent_dist_clip_obs[0][1])
    current_obs = cent_dist_clip_obs[0][1]
    saveShp([current_obs], 'polygon', path, 'current_obstacle')
    ch = createConvexhull(current_obs, facil_coords)
    
    #pick closest one to facility 
    ch_vertices = list(ch.exterior.coords)
    ch_vertices = list(set(ch_vertices))    
    ch_vertices = [x for x in ch_vertices if not Point(x) == facil_shp[0]]
    #ch_vertices = [(round(x[0],5), round(x[1], 5)) for x in ch_vertices]
    ch_v_esp = [(x, Convexpath_module.Convexpath_shapely(path, facil_shp[0], Point(x), init_clip_obs_shp).esp) for x in ch_vertices]    
    ch_v_esp = [(x[1].length, x) for x in ch_v_esp]
    ch_v_esp.sort()
    #Initial split vertex (ISV) 1
    #ISLs: two last segment of ESPs that from center to ISVs 
    v1 = ch_v_esp[0][1][0]
    #ESPs that via ISV1 
    v1_esps = [x for x in ch_v_esp if x[1][1].intersects(Point(v1).buffer(1))]
    v1_esps.sort()
    #Initial Split Line (ISL) I 
    lastSeg = list(v1_esps[0][1][1].coords)[-2:]   #This should be fixed 
    isl1 = SplitLine(lastSeg, lastSeg[0])
    v1_remain_dis = coverage_distance - v1_esps[0][0]
    
    #ISV 2 
    v2_esps = [x for x in ch_v_esp if not x in v1_esps]
    v2_esps.sort()
    v2 = v2_esps[0][1][0]
    lastSeg = list(v2_esps[0][1][1].coords)[-2:]  #This should be fixed
    isl2 = SplitLine(lastSeg, lastSeg[0])    
    v2_remain_dis = coverage_distance - v2_esps[0][0]
    
    v1_coords = list(v1_esps[0][1][1].coords)
    v2_coords = list(v2_esps[0][1][1].coords)
    
    missing_piece = v1_coords
    obs_boundary = dealt_obs[-1].boundary 
    bound_coords = list(obs_boundary.coords)
    v1_p = [x for x in v1_coords if x in bound_coords]
    v2_p = [x for x in v2_coords if x in bound_coords]  
    gg = bound_coords[bound_coords.index(v1_p[0]):bound_coords.index(v2_p[0]) + 1]
    for i in gg:
        if not i in missing_piece:
            missing_piece.append(i)
    for i in v2_coords:
        if not i in missing_piece:
            missing_piece.append(i)
        
    missing_poly = Polygon(missing_piece)
    
    cutting_box_1, cutting_box_2 = cutting_box(isl1, isl2, cover_shp, facil_coords[0])
    cutbox = [x for x in (cutting_box_1, cutting_box_2) if cent_dist_clip_obs[0][1].centroid.intersects(x)]
    other_cutbox = [x for x in (cutting_box_1, cutting_box_2) if x not in cutbox]
    saveShp([cutbox[0]], 'polygon', path, 'cutbox_0')
    saveShp(other_cutbox, 'polygon', path, 'othercutbox')
    Intersect_analysis(['cover_shp.shp', 'cutbox_0.shp'], 'initial_sector.shp')
    Intersect_analysis(['cover_shp.shp', 'othercutbox.shp'], 'other_sector.shp')
    #intial_sector = cover_shp.intersection(cutbox[0])  
    #saveShp([intial_sector], 'polygon', path, 'intial_sector')
    #other_sector = cover_shp.intersection(other_cutbox[0])
    #saveShp([other_sector], 'polygon', path, 'other_sector')
    saveShp([missing_poly], 'polygon', path, 'missing_poly')
    Union_analysis(['other_sector.shp', 'missing_poly.shp'], 'union.shp')
    Dissolve_management('union.shp', 'union_dissol.shp')
    v1_sLines = [isl1]
    v2_sLines = [isl2]
    v1_overlap = [list(v1_esps[-1][1][1].coords)[-1] ,list(v2_esps[-1][1][1].coords)[-1]]
    v2_overlap = [list(v2_esps[-1][1][1].coords)[-1] ,list(v1_esps[-1][1][1].coords)[-1]]
    v2_esps.remove(v2_esps[0])
    v1_esps.remove(v1_esps[0])
    v1_sectors = []
    v2_sectors = []
    #what about area between esp1 esp2 and obstacle? it must be included 
    
    while len(v1_esps) != 0:
        if v1_remain_dis <= 0:
            break
        
        lastSeg = list(v1_esps[0][1][1].coords)[-2:]
        #new ciclar partial cover, centered at fromNode of lastSeg 
        temp_bf_shp = Point(lastSeg[-2]).buffer(v1_remain_dis)
        #generate split line 
        sl1 = SplitLine(lastSeg, lastSeg[0])
        #split line from prior split 
        v1_sLines[-1].radius(temp_bf_shp) 
        prior_radi = v1_sLines[-1].radiusLine
        sl2 = SplitLine((lastSeg[0], list(prior_radi.coords)[1]), lastSeg[0])
        cutting_box_1, cutting_box_2 = cutting_box(sl1, sl2, temp_bf_shp, lastSeg[0])
        cutbox = [x for x in (cutting_box_1, cutting_box_2) if not cent_dist_clip_obs[0][1].centroid.intersects(x)]
        sector = temp_bf_shp.intersection(cutbox[0]) 
        v1_remain_dis -= sl1.lineObj.length
        v1_esps.remove(v1_esps[0])
        v1_sectors.append(sector)
        v1_sLines.append(sl1)
    if v1_remain_dis > 0:
        temp_bf_shp = Point(v1_overlap[0]).buffer(v1_remain_dis)
        sl1 = SplitLine(v1_overlap, v1_overlap[0])
        v1_sLines[-1].radius(temp_bf_shp) 
        prior_radi = v1_sLines[-1].radiusLine
        sl2 = SplitLine((v1_overlap[0], list(prior_radi.coords)[1]), v1_overlap[0])
        cutting_box_1, cutting_box_2 = cutting_box(sl1, sl2, temp_bf_shp, v1_overlap[0])
        cutbox = [x for x in (cutting_box_1, cutting_box_2) if not cent_dist_clip_obs[0][1].centroid.intersects(x)]
        sector = temp_bf_shp.intersection(cutbox[0]) 
        v1_sectors.append(sector)
        
    
        
    while len(v2_esps) != 0:
        if v2_remain_dis <= 0:
            break
        
        lastSeg = list(v2_esps[0][1][1].coords)[-2:]
        #new ciclar partial cover, centered at fromNode of lastSeg 
        temp_bf_shp = Point(lastSeg[-2]).buffer(v2_remain_dis)
        #generate split line 
        sl1 = SplitLine(lastSeg, lastSeg[0])
        #split line from prior split 
        v2_sLines[-1].radius(temp_bf_shp) 
        prior_radi = v2_sLines[-1].radiusLine
        sl2 = SplitLine((lastSeg[0], list(prior_radi.coords)[1]), lastSeg[0])
        cutting_box_1, cutting_box_2 = cutting_box(sl1, sl2, temp_bf_shp, lastSeg[0])
        cutbox = [x for x in (cutting_box_1, cutting_box_2) if not cent_dist_clip_obs[0][1].centroid.intersects(x)]
        sector = temp_bf_shp.intersection(cutbox[0]) 
        v2_remain_dis -= sl1.lineObj.length
        v2_esps.remove(v2_esps[0])
        v2_sectors.append(sector)
        v2_sLines.append(sl1)
    if v2_remain_dis > 0:
        temp_bf_shp = Point(v2_overlap[0]).buffer(v2_remain_dis)
        sl1 = SplitLine(v2_overlap, v2_overlap[0])
        v2_sLines[-1].radius(temp_bf_shp) 
        prior_radi = v2_sLines[-1].radiusLine
        sl2 = SplitLine((v2_overlap[0], list(prior_radi.coords)[1]), v2_overlap[0])
        cutting_box_1, cutting_box_2 = cutting_box(sl1, sl2, temp_bf_shp, v2_overlap[0])
        cutbox = [x for x in (cutting_box_1, cutting_box_2) if not cent_dist_clip_obs[0][1].centroid.intersects(x)]
        sector = temp_bf_shp.intersection(cutbox[0]) 
        v2_sectors.append(sector)
        
        
        
    intial_sector = file2Shp("initial_sector.shp")[0]
    saveShp(v1_sectors, 'polygon', path, 'v1_sectors')
    saveShp(v2_sectors, 'polygon', path, 'v2_sectors')
    Union_analysis(['v1_sectors.shp', 'v2_sectors.shp'], 'sectors.shp')
    #sectors = intial_sector.intersection(cascaded_union(v1_sectors + v2_sectors))
    #sectors_bd = sectors.boundary
    #sbd_vertices = list(sectors_bd.coords)
    #fixed_sbd_v = []
    #for i in sbd_vertices:
        #if sbd_vertices.count(i) == 1:
            #fixed_sbd_v.append(i)
        #else:
            #if not i in fixed_sbd_v:
                #fixed_sbd_v.append(i)
    #sectors = Polygon(fixed_sbd_v)
    #saveShp([sectors], 'polygon', path, 'sectors')
    saveShp(dealt_obs, 'polygon', path, 'dealtobs.shp')
    Union_analysis(['sectors.shp', 'union_dissol.shp'], 'union_2.shp')
    Dissolve_management('union_2.shp', 'union_2_dissol.shp')
    Erase_analysis('union_2_dissol.shp', 'dealtobs.shp', 'cover_shp.shp')
    cover_shp = file2Shp('cover_shp.shp')[0]
    
    

"""
1. Union: other sector + missing poly
2. Dissolve: 1
3. Union: 2 + sectors
4. Dissolve: 3

"""


