#Mk1: cannot handle obstacles beyond other obstacles. Need to change the approach fundamentally, not based on cutting box
#Mk2: Mk1 bugs are fixed. But cannot handle nonconvexity that arises during the process. Need to move toward nonconvexity. 
#     We cannot do this without consideration of non-convex obstacles 
#Mk3: non-convex obstacles 
import networkx, time, cPickle, random, math, copy, Convexpath_module, fiona
from shapely.geometry import Point, Polygon, LineString, MultiPoint, MultiPolygon, mapping
from shapely.ops import cascaded_union
from arcpy import env, Dissolve_management, Union_analysis, Intersect_analysis, Erase_analysis, MultipartToSinglepart_management, Buffer_analysis, CopyFeatures_management



env.overwriteOutput = True
path = "E:\\data\\esp_cover\\data3\\"
#"F:\\Dropbox\\research\\Distance restricted covering model\\Locating recharging station\\data7\\"
#"E:\\data\\esp_cover\\data3\\"
env.workspace = path
schema_poly = {
    'geometry': 'Polygon',
    'properties': {'id': 'int'}, }
schema_line = {
    'geometry': 'LineString',
    'properties': {'id': 'int'}, }
schema_point = {
    'geometry': 'Point',
    'properties': {'id': 'int'}, }

#user parameters 
facil = "pt1.shp"
obstacles_f = "nconvex.shp"
coverage_distance = 500  # 17424
cover = "pt_buffer.shp"

clip_obs = "clip_obs.shp"
clip_obs_single = "clip_obs_single.shp"

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
        
    
    
class SplitLine_mk2:
    """
    Only for single segment LineString object 
    fromNode is tuple of coordinates
    Intersecting line is LineString Object 
    """
    def __init__(self, line, fromNode, min_distance, inCircle):
         
        self.fromNode= fromNode
        if type(line) is LineString:
            self.lineObj = line
        else:
            self.lineObj = LineString(line)
        self.qurt = 0
        self.a = 0 
        self.b = 0 
        self.interPts = ''
        self.interLine = ''
        coords = list(self.lineObj.coords)
        coords = [x for x in coords if Point(x) != Point(self.fromNode)]
        self.toNode = coords[0]
        
        xx = self.toNode[0] - self.fromNode[0]
        yy = self.toNode[1] - self.fromNode[1]
        self.a = yy/xx
        self.b = self.fromNode[1] - self.a * self.fromNode[0]
        if xx > 0:
            self.new_toNode = (self.toNode[0] + min_distance, self.a * (self.toNode[0] + min_distance) + self.b)
        elif xx < 0:
            self.new_toNode = (self.toNode[0] - min_distance, self.a * (self.toNode[0] - min_distance) + self.b)
        ext_line = LineString([self.fromNode, self.new_toNode])
        minX, minY, maxX, maxY = inCircle.bounds 
        maxXLine = LineString([[maxX, maxY], [maxX, minY]])
        maxYLine = LineString([[minX, maxY], [maxX, maxY]])
        minXLine = LineString([[minX, maxY], [minX, minY]])
        minYLine = LineString([[minX, minY], [maxX, minY]])        
        if ext_line.intersects(maxXLine):
            self.interPts = self.intersectPoint(maxXLine)
            self.interLine = 'maxXLine'
        elif ext_line.intersects(maxYLine):
            self.interPts = self.intersectPoint(maxYLine)
            self.interLine = 'maxYLine'
        elif ext_line.intersects(minXLine):
            self.interPts = self.intersectPoint(minXLine)        
            self.interLine = 'minXLine'
        elif ext_line.intersects(minYLine):
            self.interPts = self.intersectPoint(minYLine)        
            self.interLine = 'minYLine'
            
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
    
                



def cutting_box_old (sLine1, sLine2, inCircle, circleCenter):
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
    try:
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
    except:
        writeShp([cutting_box_1], 'cutting_box_1.shp')
        raise Exception
    return cutting_box_1, cutting_box_2

#initializing objects 

def cutting_box_mk2(sLine1, sLine2, inCircle):
    minX, minY, maxX, maxY = inCircle.bounds
    boundingBox = Polygon([(minX, minY), (minX, maxY), (maxX, maxY), (maxX, minY)])
    cutting_box = []
    if sLine1.interLine == 'maxXLine':
        if sLine2.interLine == 'maxXLine':
            cutting_box = Polygon([sLine1.fromNode, sLine1.interPts, sLine2.interPts, sLine2.fromNode])
        elif sLine2.interLine == 'minYLine':
            cutting_box = Polygon([sLine1.fromNode, sLine1.interPts, (maxX, minY), sLine2.interPts, sLine2.fromNode])
        elif sLine2.interLine == 'minXLine':
            cutting_box = Polygon([sLine1.fromNode, sLine1.interPts, (maxX, minY), (minX, minY), sLine2.interPts, sLine2.fromNode])
        elif sLine2.interLine == 'maxYLine':
            cutting_box = Polygon([sLine1.fromNode, sLine1.interPts, (maxX, maxY), sLine2.interPts, sLine2.fromNode])
    elif sLine1.interLine == 'minYLine':
        if sLine2.interLine == 'minYLine':
            cutting_box = Polygon([sLine1.fromNode, sLine1.interPts, sLine2.interPts, sLine2.fromNode])
        elif sLine2.interLine == 'minXLine':
            cutting_box = Polygon([sLine1.fromNode, sLine1.interPts, (minX, minY), sLine2.interPts, sLine2.fromNode])
        elif sLine2.interLine == 'maxYLine':
            cutting_box = Polygon([sLine1.fromNode, sLine1.interPts, (minX, minY), (minX, maxY), sLine2.interPts, sLine2.fromNode])
        elif sLine2.interLine == 'maxXLine':
            cutting_box = Polygon([sLine1.fromNode, sLine1.interPts, (maxX, minY), sLine2.interPts, sLine2.fromNode])
    elif sLine1.interLine == 'minXLine':
        if sLine2.interLine == 'minXLine':
            cutting_box = Polygon([sLine1.fromNode, sLine1.interPts, sLine2.interPts, sLine2.fromNode])
        elif sLine2.interLine == 'maxYLine':
            cutting_box = Polygon([sLine1.fromNode, sLine1.interPts, (minX, maxY), sLine2.interPts, sLine2.fromNode])
        elif sLine2.interLine == 'maxXLine':
            cutting_box = Polygon([sLine1.fromNode, sLine1.interPts, (minX, maxY), (maxX, maxY), sLine2.interPts, sLine2.fromNode])
        elif sLine2.interLine == 'minYLine':
            cutting_box = Polygon([sLine1.fromNode, sLine1.interPts, (minX, minY), sLine2.interPts, sLine2.fromNode])
    elif sLine1.interLine == 'maxYLine':
        if sLine2.interLine == 'maxYLine':
            cutting_box = Polygon([sLine1.fromNode, sLine1.interPts, sLine2.interPts, sLine2.fromNode])
        elif sLine2.interLine == 'maxXLine':
            cutting_box = Polygon([sLine1.fromNode, sLine1.interPts, (maxX, maxY), sLine2.interPts, sLine2.fromNode])
        elif sLine2.interLine == 'minYLine':
            cutting_box = Polygon([sLine1.fromNode, sLine1.interPts, (maxX, maxY), (maxX, minY), sLine2.interPts, sLine2.fromNode])
        elif sLine2.interLine == 'minXLine':
            cutting_box = Polygon([sLine1.fromNode, sLine1.interPts, (minX, maxY), sLine2.interPts, sLine2.fromNode])
    
    cutting_box_2 = boundingBox.difference(cutting_box)
    return cutting_box, cutting_box_2    

def init_cuuting_box (sLine1, sLine2, inCircle):
    minX, minY, maxX, maxY = inCircle.bounds 
    boundingBox = Polygon([(minX, minY), (minX, maxY), (maxX, maxY), (maxX, minY)])
    cutting_box = []
    if sLine1.interLine == 'maxXLine':
        if sLine2.interLine == 'maxXLine':
            cutting_box = Polygon([sLine1.toNode, sLine1.interPts, sLine2.interPts, sLine2.toNode])
        elif sLine2.interLine == 'minYLine':
            cutting_box = Polygon([sLine1.toNode, sLine1.interPts, (maxX, minY), sLine2.interPts, sLine2.toNode])
        elif sLine2.interLine == 'minXLine':
            cutting_box = Polygon([sLine1.toNode, sLine1.interPts, (maxX, minY), (minX, minY), sLine2.interPts, sLine2.toNode])
        elif sLine2.interLine == 'maxYLine':
            cutting_box = Polygon([sLine1.toNode, sLine1.interPts, (maxX, maxY), sLine2.interPts, sLine2.toNode])
    elif sLine1.interLine == 'minYLine':
        if sLine2.interLine == 'minYLine':
            cutting_box = Polygon([sLine1.toNode, sLine1.interPts, sLine2.interPts, sLine2.toNode])
        elif sLine2.interLine == 'minXLine':
            cutting_box = Polygon([sLine1.toNode, sLine1.interPts, (minX, minY), sLine2.interPts, sLine2.toNode])
        elif sLine2.interLine == 'maxYLine':
            cutting_box = Polygon([sLine1.toNode, sLine1.interPts, (minX, minY), (minX, maxY), sLine2.interPts, sLine2.toNode])
        elif sLine2.interLine == 'maxXLine':
            cutting_box = Polygon([sLine1.toNode, sLine1.interPts, (maxX, minY), sLine2.interPts, sLine2.toNode])
    elif sLine1.interLine == 'minXLine':
        if sLine2.interLine == 'minXLine':
            cutting_box = Polygon([sLine1.toNode, sLine1.interPts, sLine2.interPts, sLine2.toNode])
        elif sLine2.interLine == 'maxYLine':
            cutting_box = Polygon([sLine1.toNode, sLine1.interPts, (minX, maxY), sLine2.interPts, sLine2.toNode])
        elif sLine2.interLine == 'maxXLine':
            cutting_box = Polygon([sLine1.toNode, sLine1.interPts, (minX, maxY), (maxX, maxY), sLine2.interPts, sLine2.toNode])
        elif sLine2.interLine == 'minYLine':
            cutting_box = Polygon([sLine1.toNode, sLine1.interPts, (minX, minY), sLine2.interPts, sLine2.toNode])
    elif sLine1.interLine == 'maxYLine':
        if sLine2.interLine == 'maxYLine':
            cutting_box = Polygon([sLine1.toNode, sLine1.interPts, sLine2.interPts, sLine2.toNode])
        elif sLine2.interLine == 'maxXLine':
            cutting_box = Polygon([sLine1.toNode, sLine1.interPts, (maxX, maxY), sLine2.interPts, sLine2.toNode])
        elif sLine2.interLine == 'minYLine':
            cutting_box = Polygon([sLine1.toNode, sLine1.interPts, (maxX, maxY), (maxX, minY), sLine2.interPts, sLine2.toNode])
        elif sLine2.interLine == 'minXLine':
            cutting_box = Polygon([sLine1.toNode, sLine1.interPts, (minX, maxY), sLine2.interPts, sLine2.toNode])
    writeShp([cutting_box] , 'cutting_box.shp')    
    cutting_box_2 = boundingBox.difference(cutting_box)
    return cutting_box, cutting_box_2
        
    

    
 

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
    all elements must be same type 
    """
    for i in in_shp:
        if i.is_empty:
            print str(i) + ' is empty'
            raise Exception
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

obstacles_shp = openShp(obstacles_f)
facil_shp = openShp(facil)
facil_coords = list(facil_shp[0].coords)

#generate initial, circular coverage 

cover_shp = facil_shp[0].buffer(coverage_distance)
initial_cover_shp = facil_shp[0].buffer(coverage_distance)
writeShp([cover_shp], 'cover_shp.shp')
writeShp([initial_cover_shp], 'initial_cover.shp')
IminX, IminY, ImaxX, ImaxY = cover_shp.bounds 
ImaxXLine = LineString([[ImaxX, ImaxY], [ImaxX, IminY]])
ImaxYLine = LineString([[IminX, ImaxY], [ImaxX, ImaxY]])
IminXLine = LineString([[IminX, ImaxY], [IminX, IminY]])
IminYLine = LineString([[IminX, IminY], [ImaxX, IminY]])
boundingBox = Polygon([(IminX, IminY), (IminX, ImaxY), (ImaxX, ImaxY), (ImaxX, IminY)])

#generate set of obstacles for circular cover 


clip_obs_shp = [x.intersection(cover_shp) for x in obstacles_shp]
clip_obs_shp = [x for x in clip_obs_shp if not x.is_empty]
init_clip_obs_shp = [x for x in clip_obs_shp if not x.is_empty]

terminate = False
#begin main loop 
dealt_obs = []
count = 0

while terminate == False:
    #first loop: generate 
    #For each iteration, list of obstacles is updated using updated coverage. 
    #Remove obstacles that cannot be covered with the updated coverage 

    cover_shp = openShp('cover_shp.shp')[0]
    clip_obs_shp = [x for x in clip_obs_shp if not x in dealt_obs]   
    clip_obs_shp = [x.intersection(cover_shp) for x in clip_obs_shp]
    clip_obs_shp = [x for x in clip_obs_shp if not x.is_empty]

    if len(clip_obs_shp) == 0:
        break
    
    #convex hulls for clipped obstacles
    #sort them as distance between facility and centorids 
    cent_clip_obs = [(x.centroid, x) for x in clip_obs_shp]
    cent_dist_clip_obs = [(x[0].distance(facil_shp[0]), x[1]) for x in cent_clip_obs]
    cent_dist_clip_obs.sort()
    dealt_obs.append(cent_dist_clip_obs[0][1])
    current_obs = cent_dist_clip_obs[0][1]
    writeShp([current_obs], 'current_obs.shp')

    ch = createConvexhull(current_obs, facil_coords)
    writeShp([ch], 'ch.shp')
    
    #ch_vertices for processing 
    ch_vertices = list(ch.exterior.coords)
    ch_vertices = list(set(ch_vertices))    
    ch_vertices = [x for x in ch_vertices if not Point(x) == facil_shp[0]]
    ch_v_esp = [(x, Convexpath_module.Convexpath_shapely(path, facil_shp[0], Point(x), clip_obs_shp + dealt_obs).esp) 
                for x in ch_vertices] 
    ch_v_esp = [(x[1].length, x) for x in ch_v_esp]
    ch_v_esp.sort()    
    obs_vertices = list(current_obs.boundary.coords)
    #pick two vertices (ISVs)
    ch_boundary = list(ch.boundary.coords)
    ch_lines = []
    for i in range(len(ch_boundary) -1):
        ch_lines.append(LineString((ch_boundary[i], ch_boundary[i+1])))
    ch_lines = [x for x in ch_lines if x.intersects(facil_shp[0].buffer(1))]
    isl_vertices = []    
    for line in ch_lines:
        l = list(line.coords)
        isl_vertices.extend(l)
    isl_vertices = [x for x in isl_vertices if not Point(x) == facil_shp[0]]
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #extract vertices of current_obs that are interior of convex hull 
    #c_obs_vertices = list(current_obs.exterior.coords)
    
    #c_obs_v_esps = [(x, Convexpath_module.Convexpath_shapely(path, facil_shp[0], clip_obs+dealt_obs).esp) 
                    #for x in c_obs_vertices]    




















    
    v1 = isl_vertices[0]
    v1_esp = Convexpath_module.Convexpath_shapely(path, facil_shp[0], Point(v1), clip_obs_shp + dealt_obs).esp
    v1_esp_vertices = list(v1_esp.coords)
    v1_esp_vertices = [x for x in v1_esp_vertices if x in obs_vertices]
    if len(v1_esp_vertices) >= 2:
        v1_esp_vertices = [x for x in v1_esp_vertices if not Point(v1).equals(Point(x))]
        esps = [(x, Convexpath_module.Convexpath_shapely(path, facil_shp[0], Point(x), clip_obs_shp + dealt_obs).esp) for x in v1_esp_vertices]
        esps = [(x[1].length, x[0], x[1]) for x in esps]
        esps.sort()
        v1 = esps[0][1]
        v1_esp = esps[0][2]
        for i in esps:
            ch_v_esp.append((i[0], (i[1], i[2])))
    v1_esps = [x for x in ch_v_esp if x[1][1].intersects(Point(v1).buffer(1))]
    v1_esps.sort()
   
    #Initial split vertex (ISV) 1
    #ESPs that via ISV1 
    #Initial Split Line (ISL) I 
    lastSeg = list(v1_esps[0][1][1].coords)[-2:]  
    isl1 = SplitLine_mk2(lastSeg, lastSeg[0], coverage_distance, initial_cover_shp)
    v1_remain_dis = coverage_distance - v1_esps[0][0]
    writeShp([isl1.lineObj], 'isl1.shp')
    
    
    #ISV 2 
    v2 = isl_vertices[1]
    v2_esp = Convexpath_module.Convexpath_shapely(path, facil_shp[0], Point(v2), clip_obs_shp + dealt_obs).esp
    v2_esp_vertices = list(v2_esp.coords)
    v2_esp_vertices = [x for x in v2_esp_vertices if x in obs_vertices]
    if len(v2_esp_vertices) >= 2:
        v2_esp_vertices = [x for x in v2_esp_vertices if not Point(v2).equals(Point(x))]
        esps = [(x, Convexpath_module.Convexpath_shapely(path, facil_shp[0], Point(x), clip_obs_shp + dealt_obs).esp) for x in v2_esp_vertices]
        esps = [(x[1].length, x[0], x[1]) for x in esps]
        esps.sort()
        v2 = esps[0][1]
        v2_esp = esps[0][2]
        for i in esps:
            ch_v_esp.append((i[0], (i[1], i[2])))
            
            
    v2_esps = [x for x in ch_v_esp if x[1][1].intersects(Point(v2).buffer(1))]
    v2_esps.sort()
    lastSeg = list(v2_esps[0][1][1].coords)[-2:]  
    
    isl2 = SplitLine_mk2(lastSeg, lastSeg[0], coverage_distance, initial_cover_shp)
    v2_remain_dis = coverage_distance - v2_esps[0][0]
    writeShp([isl2.lineObj], 'isl2.shp')
    v1_coords = list(v1_esps[0][1][1].coords)
    v2_coords = list(v2_esps[0][1][1].coords)
        
    cutting_box_1, cutting_box_2 = init_cuuting_box(isl1, isl2, initial_cover_shp)
    #need to be fixed: centroid is not suitable for this. 
    writeShp([cutting_box_1, cutting_box_2], 'cutting_boxes.shp')
    cutbox = [x for x in (cutting_box_1, cutting_box_2) if not facil_shp[0].intersects(x)]
    other_cutbox = [x for x in (cutting_box_1, cutting_box_2) if x not in cutbox]    
    initial_sector = cover_shp.intersection(cutbox[0])
    bounding_sector = initial_sector.union(ch)
    writeShp([bounding_sector],'bounding_sector.shp')
    other_sector = cover_shp.intersection(other_cutbox[0])
    writeShp([initial_sector], 'initial_sector.shp')
    if other_sector.type == 'MultiPolygon':
        writeShp(list(other_sector), 'other_sector.shp')
    else:
        writeShp([other_sector], 'other_sector.shp')
    
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
        writeShp([temp_bf_shp], 'buffer_v1.shp')
        #generate split line 
        sl1 = SplitLine_mk2(lastSeg, lastSeg[0], coverage_distance, temp_bf_shp)
        
        #split line from prior split 
        
        prior_radi = v1_sLines[-1].cuttingLine  
        sl2 = SplitLine_mk2((lastSeg[0], list(prior_radi.coords)[1]), lastSeg[0], coverage_distance, temp_bf_shp)
        
        cutting_box_1, cutting_box_2 = cutting_box_mk2(sl1, sl2, temp_bf_shp)
        writeShp([cutting_box_1, cutting_box_2], 'small_cuttingBox.shp')
        
        s_cutbox = [cutting_box_1, cutting_box_2]
        sector = temp_bf_shp.intersection(s_cutbox[0]) 
        v1_remain_dis -= sl1.lineObj.length
        v1_esps.remove(v1_esps[0])
        v1_sectors.append(sector)
        v1_sLines.append(sl1)
        #print "end of sub loop"
        #raw_input()
    if v1_remain_dis > 0:
        temp_bf_shp = Point(v1_overlap[0]).buffer(v1_remain_dis)
        sl1 = SplitLine_mk2(v1_overlap, v1_overlap[0], coverage_distance, temp_bf_shp)
        prior_radi = v1_sLines[-1].cuttingLine  
        sl2 = SplitLine_mk2((v1_overlap[0], list(prior_radi.coords)[1]), v1_overlap[0], coverage_distance, temp_bf_shp)
        cutting_box_1, cutting_box_2 = cutting_box_mk2(sl1, sl2, temp_bf_shp)
        
        s_cutbox = [cutting_box_1, cutting_box_2]
        sector = temp_bf_shp.intersection(s_cutbox[0]) 
        v1_sectors.append(sector)
        
    while len(v2_esps) != 0:
        if v2_remain_dis <= 0:
            break
        
        lastSeg = list(v2_esps[0][1][1].coords)[-2:]
        #new ciclar partial cover, centered at fromNode of lastSeg 
        temp_bf_shp = Point(lastSeg[-2]).buffer(v2_remain_dis)
        writeShp([temp_bf_shp], 'buffer_v2.shp')
        #generate split line 
        sl1 = SplitLine_mk2(lastSeg, lastSeg[0], coverage_distance, temp_bf_shp)
        prior_radi = v2_sLines[-1].cuttingLine  
        sl2 = SplitLine_mk2((lastSeg[0], list(prior_radi.coords)[1]), lastSeg[0], coverage_distance, temp_bf_shp)
        cutting_box_1, cutting_box_2 = cutting_box_mk2(sl1, sl2, temp_bf_shp)
        writeShp([cutting_box_1, cutting_box_2], 'small_cuttingBox.shp')
        s_cutbox = [cutting_box_1, cutting_box_2]        
        sector = temp_bf_shp.intersection(s_cutbox[0]) 
        v2_remain_dis -= sl1.lineObj.length
        v2_esps.remove(v2_esps[0])
        v2_sectors.append(sector)
        v2_sLines.append(sl1)
        #print "end of sub loop"
        #raw_input()        
    if v2_remain_dis > 0:
        temp_bf_shp = Point(v2_overlap[0]).buffer(v2_remain_dis)
        sl1 = SplitLine_mk2(v2_overlap, v2_overlap[0], coverage_distance, temp_bf_shp)
        prior_radi = v2_sLines[-1].cuttingLine  
        sl2 = SplitLine_mk2((v2_overlap[0], list(prior_radi.coords)[1]), v2_overlap[0], coverage_distance, temp_bf_shp)
        cutting_box_1, cutting_box_2 = cutting_box_mk2(sl1, sl2, temp_bf_shp)
        s_cutbox = [cutting_box_1, cutting_box_2]        
        sector = temp_bf_shp.intersection(s_cutbox[0]) 
        v2_sectors.append(sector)
            
 
    
    v1_sectors = [initial_sector.intersection(x) for x in v1_sectors]    
    v2_sectors = [initial_sector.intersection(x) for x in v2_sectors]    
    v1_sectors = [x for x in v1_sectors if not x.is_empty]
    v2_sectors = [x for x in v2_sectors if not x.is_empty]    
    v1_sectors = [x for x in v1_sectors if x.type == 'Polygon']
    v2_sectors = [x for x in v2_sectors if x.type == 'Polygon']
    union_sectors = []
    if len(v1_sectors) != 0:
        writeShp(v1_sectors, 'v1_sectors.shp')
        union_sectors.append('v1_sectors.shp')
    if len(v2_sectors) != 0:
        writeShp(v2_sectors, 'v2_sectors.shp')
        union_sectors.append('v2_sectors.shp')
    if len(v1_sectors) == 0:
        if len(v2_sectors) == 0:
            writeShp(dealt_obs, 'dealtobs.shp')
            CopyFeatures_management('cover_shp.shp', 'cover_1_shp.shp')
            Erase_analysis('cover_1_shp.shp', 'dealtobs.shp', 'cover_shp.shp')
            print 'end of loop, no sectors'
            #raw_input()
        else:
            Union_analysis(union_sectors, 'sectors.shp')
            
            writeShp(dealt_obs, 'dealtobs.shp')
            
            Union_analysis(['sectors.shp', 'other_sector.shp'], 'union_2.shp')
            Dissolve_management('union_2.shp', 'union_2_dissol.shp')
            Erase_analysis('union_2_dissol.shp', 'dealtobs.shp', 'cover_shp.shp')
            print "end of loop"
        
            #raw_input()            
    else:
        Union_analysis(union_sectors, 'sectors.shp')
        
        writeShp(dealt_obs, 'dealtobs.shp')
        
        Union_analysis(['sectors.shp', 'other_sector.shp'], 'union_2.shp')
        Dissolve_management('union_2.shp', 'union_2_dissol.shp')
        Erase_analysis('union_2_dissol.shp', 'dealtobs.shp', 'cover_shp.shp')
        print "end of loop"
    
        #raw_input()         
    
   
    

"""
1. Union: other sector + missing poly
2. Dissolve: 1
3. Union: 2 + sectors
4. Dissolve: 3

"""


