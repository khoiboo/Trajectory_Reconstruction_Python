'''
Method of coplanarity to determine a line in 3d space
given direction vectors from multiple stations.

Created on 19.01.2021

@author: Anastasios Margonis & Georgios Malissiovas
'''

import numpy as np
import csv
#import pandas
import math
#import geopandas as gpd
from shapely.geometry import Point, Polygon, MultiPoint, LineString
from pyproj import Geod
#import time
from io import StringIO   
import sys 
import argparse

import Line_intersection
import Functions
np.set_printoptions(precision=22)     

SEARCH_RADIUS_AZIMUTH = 5000
GRID_SIZE_2D = 400 #400 points on each side
GRID_SPACING_2D = 5000 # unit in meters
a = 6378137
e = 0.081819190842622
f = 0.0034
b = 6.3568e6
e2 = math.sqrt((math.pow(a, 2) - math.pow(b, 2)) / math.pow(b, 2))
EPSILON = 0.0001
EARTH_RADIUS_METER = 6371000
SEARCH_RADIUS_AZIMUTH = 5000
OBSERVATION_RANGE = 700000
MAX_TRAJECTORY_LENGTH = 250 # kilometers 

# For the determination of the X_MAX, X_MIN, Y_MAX, Y_MIN, Z_MAX and Z_MIN, 
# the following 4 reference points (lat, long) are used: 
# (62.041, -5.972), (36.127, -5.352), (37.983, 23.732) and (61.679, 27.283)
X_MIN = 2696 #2090
X_MAX = 5135 #5139
Y_MIN = -481 #-1564
Y_MAX = 2025 #2067
Z_MIN = 3739 #3426
Z_MAX = 5610 #5735
GRID_SPACING_3D = 6
MAXIMUM_FIREBALL_ALTITUDE = 100
MINIMUM_FIREBALL_ALTITUDE = 30
EARTH_RADIUS = 6371
SEARCH_RADIUS_ELEVATION = 15
PARAMETER = 3000

class Observation(object):
    def __init__(self, observation_ID, observer_name, lat, lon,
                 start_azimuth, start_elevation, end_azimuth, end_elevation):
        assert type(lat) is float or int, "Latitude value is not a float"
        assert type(lon) is float or int, "Longitude value is not a float"
        assert type(start_azimuth) is float or int, "Start azimuth value is not a float"
        assert type(end_azimuth) is float or int, "End azimuth value is not a float"
        assert type(start_elevation) is float or int, "Start elevation value is not a float"
        assert type(end_elevation) is float or int, "End elevation value is not a float"        
        
        self.observation_ID = observation_ID
        self.observer_name = observer_name
        self.lat = lat
        self.lon = lon
        self.start_azimuth   = start_azimuth
        self.start_elevation = start_elevation
        self.end_azimuth   = end_azimuth
        self.end_elevation = end_elevation
        
    def get_observation_ID(self):
        return self.observation_ID
    
    def get_latitude(self):
        return self.lat
    
    def get_longitude(self):
        return self.lon
    
    def get_start_azimuth(self):
        return self.start_azimuth
    
    def get_start_elevation(self):
        return self.start_elevation
    
    def get_end_azimuth(self):
        return self.end_azimuth
    
    def get_end_elevation(self):
        return self.end_elevation
    

class Point3D(object):
    def __init__(self,x,y,z):
        assert type(x) is float or int, "X coordinate is not a float"
        assert type(y) is float or int, "Y coordiante is not a float"
        assert type(z) is float or int, "Z coordinate is not a float"
        
        self.X = x
        self.Y = y
        self.Z = z
    
    def get_X(self):
        return self.X
    
    def get_Y(self):
        return self.Y
    
    def get_Z(self):
        return self.Z   

    def __repr__(self):
        return (str(self.X) + ", " + str(self.Y) + ", " + str(self.Z))


class Vector3D(object):
    def __init__(self,x,y,z):
        assert type(x) is float or int, "X coordinate is not a float"
        assert type(y) is float or int, "Y coordiante is not a float"
        assert type(z) is float or int, "Z coordinate is not a float"
        
        self.X = x
        self.Y = y
        self.Z = z
    
    def get_X(self):
        return self.X
    
    def get_Y(self):
        return self.Y
    
    def get_Z(self):
        return self.Z
    

class Line3D(object):
    def __init__(self, point, vector):
        self.point = point
        self.directionVector = vector
        
    def get_point(self):
        return self.point
    
    def get_direction_vector(self):
        return self.direction_vector
    

class Plane3D(object):
    def __init__(self, a, b, c, d):
        self.A = a
        self.B = b
        self.C = c
        self.D = d
        
    def get_A(self):
        return self.A
    
    def get_B(self):
        return self.B
    
    def get_C(self):
        return self.C
    
    def get_D(self):
        return self.D
    
def to_radian(theta):    
    return (theta * math.pi)/180

def to_degree(theta):    
    return (theta*180)/math.pi

# This function identifies a new point based on distance 
# and direction from the origin point    
def move_geodetically(origin, azimuth, distance_in_meter):
    geoid = Geod(ellps='WGS84')
    long_new, lat_new, return_az = geoid.fwd(origin.x, origin.y, azimuth, distance_in_meter)
    return Point(long_new, lat_new)

# This function creates a regular polygon in the shape of a circle 
# surrounding a point. The parameter "number_of_points" determines 
# the number of vertices
def create_buffer(center, radius, number_of_points):
    assert number_of_points >=4, "not enough points to form a proper circle"    
    geoid = Geod(ellps='WGS84')
    step = 360/number_of_points
    coords = []
    
    for i in range(0, number_of_points):
        point = move_geodetically(center, i*step,radius)
        coords.append((point.x,point.y))
        
    buffer = Polygon(coords)    
    return buffer  
    
# This function computes the rank of each point on the grid by 
# counting the number of observation triangles that overlay this point
def update_rank_of_gridpoint(point, list_observation_ball, polygon_array, map_ball_rank):    
    geoid = Geod(ellps='WGS84')    
    buffer = create_buffer(point,SEARCH_RADIUS_AZIMUTH,10)    
    
    for i in range(0, len(polygon_array)):
        if (polygon_array[i].intersects(buffer)):            
            if (map_ball_rank.get(point) is None):
                map_ball_rank[point] = 1
            else:
                current_value = map_ball_rank.get(point)
                map_ball_rank[point] = current_value + 1
            list_observation_ball[i].append(point)
            #buffer_point_table[point] = point
    
    return 1                
    
def find_second_max(input_list):
    max_value = max(input_list)
    second_max_value = -1
     
    for i in range(0, len(input_list)):
        val = input_list[i]
        if (val > second_max_value and val < max_value):
            second_max_value = val
     
    return second_max_value     

# Reference point (lon, lat) for the function to rank the observation by azimuth
TOPLEFT_POINT = Point(-9.754214, 55.492227)
# This function returns a list of observations with the highest 
# probabilities of containing the trajectory of the observed 
# fireball by ranking them based on the azimuth values
def rank_observation_by_azimuth_optimized(observation_list, distance, 
                                          points_with_highest_ranking_by_azimuth):
    reliable_observations_azimuth = []    
    list_observation_ball = []
    #buffer_point_table = {}
    polygon_array = []    
    geoid = Geod(ellps='WGS84')
    
    for obs in observation_list:
        list_ball_of_observation = []
        list_observation_ball.append(list_ball_of_observation)
                
        observer = Point(obs.get_longitude(), obs.get_latitude()) # Longitude first, then latitude
        print("location of observer is",observer.y, observer.x)
        endpoint1 = move_geodetically(observer, obs.get_start_azimuth(), distance)
        endpoint2 = move_geodetically(observer, obs.get_end_azimuth(), distance)
        
        triangle = Polygon([(observer.x, observer.y), (endpoint1.x, endpoint1.y), 
                            (endpoint2.x, endpoint2.y)])        
        
        poly_area, poly_perimeter = geoid.geometry_area_perimeter(triangle)
        print("Observation triangle is", triangle, "area", poly_area )
        polygon_array.append(triangle)        
    
    map_ball_rank = {}
    list_point = []
    list_point.append(TOPLEFT_POINT)
    
    print("Initially, Size of mapBallRank is ", len(map_ball_rank))
    #intersectionCount = 0
    
    for j in range(0, GRID_SIZE_2D):
        print("j is now", j,"/", GRID_SIZE_2D)        
        for i in range(0, GRID_SIZE_2D - 1):
            new_point = move_geodetically(list_point[len(list_point)-1], 90, GRID_SPACING_2D)    
            #print("new_point is", new_point.y, new_point.x)
            list_point.append(new_point)
            update_rank_of_gridpoint(new_point,list_observation_ball,polygon_array,map_ball_rank)
        
        current_row_start = list_point[len(list_point) - GRID_SIZE_2D]
        
        if (j%2 == 0):            
            new_point1 = move_geodetically(current_row_start, 135, GRID_SPACING_2D)
            list_point.append(new_point1)
            update_rank_of_gridpoint(new_point1,list_observation_ball,polygon_array,map_ball_rank)
        else:            
            new_point2 = move_geodetically(current_row_start, 225, GRID_SPACING_2D)
            list_point.append(new_point2)
            update_rank_of_gridpoint(new_point2,list_observation_ball,polygon_array,map_ball_rank)
    
    print("Number of points generated is", len(list_point))
    print("Size of mapBallRank is ", len(map_ball_rank))
    #print("Number of intersection is", intersectionCount)
    
    highest_rank_overall = max(map_ball_rank.values())    
    second_highest_rank_overall = find_second_max(list(map_ball_rank.values()))
    print("Highest ranking overall is", highest_rank_overall, ", second highest is", second_highest_rank_overall)
        
    for i in range(0, len(observation_list)):
        max_rank_of_observation = 0
        if (len(list_observation_ball[i]) > 0):
            for j in range(0, len (list_observation_ball[i])):
                if (map_ball_rank[list_observation_ball[i][j]] > max_rank_of_observation):
                    max_rank_of_observation = map_ball_rank[list_observation_ball[i][j]]
        
        if ((max_rank_of_observation == highest_rank_overall) 
            or (max_rank_of_observation == second_highest_rank_overall)):
            reliable_observations_azimuth.append(observation_list[i])
            
    for point in map_ball_rank.keys():
        if (map_ball_rank[point] == highest_rank_overall):
            points_with_highest_ranking_by_azimuth.append(point)
    
    print("Number of points with highest azimuth ranking is", len(points_with_highest_ranking_by_azimuth))
    
    return reliable_observations_azimuth
    
# This function converts a local direction vector defined by azimuth 
# and elevation angle from a point on Earth into a vector 
# in the ECEF reference system 
def convert_angular_direction_vector_to_cartesian_vector(lat, lon, azimuth, elevation):
    azimuth_south = ((180 - azimuth) * math.pi) / 180
    elevation_positive_z = ((90 - elevation) * math.pi) / 180
    
    x = math.cos(azimuth_south) * math.sin(elevation_positive_z)
    y = math.sin(azimuth_south) * math.sin(elevation_positive_z)
    z = math.cos(elevation_positive_z)
    
    alpha = to_radian(lon)
    beta = to_radian(lat)
    
    component1 = x * math.sin(beta) * math.cos(alpha) - y * math.sin(alpha) + z * math.cos(alpha) * math.cos(beta)
    component2 = x * math.sin(beta) * math.sin(alpha) + y * math.cos(alpha) + z * math.cos(beta) * math.sin(alpha)
    component3 = (-1) * x * math.cos(beta)            + 0                   + z * math.sin(beta)
    
    sum_of_square = math.pow(component1, 2) + math.pow(component2, 2) + math.pow(component3, 2)
    component1 = component1/sum_of_square
    component2 = component2/sum_of_square
    component3 = component3/sum_of_square
    
    return Vector3D(component1, component2, component3)

# This function converts a point determined by lat, long and altitude (in meters) 
# into a point in the ECEF reference system
def convert_lat_long_alt_to_cartesian_coordinates(lat_degree, long_degree, height_in_meters):
    lat = to_radian(lat_degree)
    lon = to_radian(long_degree)
    alt = height_in_meters
    
    N = a / math.sqrt(1 - math.pow(math.sin(lat) * e, 2))
    
    xx = (N + alt) * math.cos(lat)* math.cos(lon)
    yy = (N + alt) * math.cos(lat)* math.sin(lon)
    zz = ((1 - e*e)*N + alt) * math.sin(lat)

    return Point3D(xx/1000, yy/1000, zz/1000)

# This function computes the cross product of 2 vectors    
def cross_product(a, b):
    x = a.get_Y()*b.get_Z() - a.get_Z()*b.get_Y()
    y = a.get_Z()*b.get_X() - a.get_X()*b.get_Z()
    z = a.get_X()*b.get_Y() - a.get_Y()*b.get_X()
    
    return Vector3D(x, y, z)

# This function computes the distance in kilometers between 2 points in 3D space
def compute_distance_between_points(point1, point2):
    distance = math.sqrt( math.pow(point1.get_X() - point2.get_X(), 2) 
                          + math.pow(point1.get_Y() - point2.get_Y(), 2) 
                          + math.pow(point1.get_Z() - point2.get_Z(), 2))    
    return distance 

# This function computes the distance in kilometers between a point and a plane in 3D space
def compute_distance_point_to_plane(point, plane):
    numerator = abs(plane.get_A()*point.get_X() + plane.get_B()*point.get_Y() 
                    + plane.get_C()*point.get_Z() + plane.get_D())
    denominator = math.sqrt( math.pow(plane.get_A(), 2) 
                             + math.pow(plane.get_B(), 2) 
                             + math.pow(plane.get_C(), 2))
    return numerator / denominator

# Define a plane in 3D space in the form of A*x + B*y + C*z + D = 0
def create_plane_from_vector_and_point(vector_a, vector_b, point):
    orthogonal_vector = cross_product(vector_a, vector_b)
    component_D = (-1) * (orthogonal_vector.get_X() * point.get_X() 
                          + orthogonal_vector.get_Y() * point.get_Y() 
                          + orthogonal_vector.get_Z() * point.get_Z())
                          
    plane = Plane3D(orthogonal_vector.get_X(), 
                    orthogonal_vector.get_Y(), 
                    orthogonal_vector.get_Z(), 
                    component_D)
    return plane

def project_point_to_plane(point, plane):
    numerator = (-1)*(plane.get_A()*point.get_X() + plane.get_B()*point.get_Y() 
                      + plane.get_C()*point.get_Z() + plane.get_D())
    denominator = (math.pow(plane.get_A(), 2) + math.pow(plane.get_B(), 2) + math.pow(plane.get_C(), 2))
    t = numerator / denominator
    
    result_point = Point3D(plane.get_A() * t + point.get_X(), 
                           plane.get_B() * t + point.get_Y(), 
                           plane.get_C() * t + point.get_Z())
    
    return result_point

def dot_product(a, b):
    return a.get_X() * b.get_X() + a.get_Y() * b.get_Y() + a.get_Z() * b.get_Z()

def compute_vector_magnitude(vector):
    return math.sqrt(math.pow(vector.get_X(), 2) 
                     + math.pow(vector.get_Y(), 2) 
                     + math.pow(vector.get_Z(), 2))

# This function computes the angle in degrees between 2 vectors in 
# 3D space. The resulted angle is at most 180 degrees
def compute_angle_between_vectors_using_cosine(a, b):
    dot_product_value = dot_product(a,b)
    
    magnitude_a_vector = compute_vector_magnitude(a)
    magnitude_b_vector = compute_vector_magnitude(b)
    
    cos_theta = dot_product_value / (magnitude_a_vector * magnitude_b_vector)
    if (cos_theta > 1 or cos_theta < -1):
        print("cos_theta is", cos_theta)
        print(a.get_X(), a.get_Y(), a.get_Z(),b.get_X(), b.get_Y(), b.get_Z())
        print("dot product is", dot_product_value)
        print("magnitude of a is", magnitude_a_vector)
        print("magnitude of b is", magnitude_b_vector)
        cos_theta = 1
    
    return to_degree(math.acos(cos_theta)) 

#This function determines whether the projection of a point into a plane 
# lies inside a triangle defined by  observer, the start and end vector 
def determine_projected_point_in_observation_triangle(high_point, start_vector, end_vector, observer_3D_location):
    whole_observation_plane = create_plane_from_vector_and_point(start_vector, end_vector, observer_3D_location)
    projected_point = project_point_to_plane(high_point, whole_observation_plane)
    
    angle_between_start_and_end_vector = compute_angle_between_vectors_using_cosine(start_vector, end_vector)
    
    vector_observer_to_projected_point = Vector3D(projected_point.get_X() - observer_3D_location.get_X(), 
                                                  projected_point.get_Y() - observer_3D_location.get_Y(), 
                                                  projected_point.get_Z() - observer_3D_location.get_Z())
    
    first_angle = compute_angle_between_vectors_using_cosine(start_vector, vector_observer_to_projected_point)
    second_angle = compute_angle_between_vectors_using_cosine(vector_observer_to_projected_point, end_vector)
    
    return (abs(angle_between_start_and_end_vector - (first_angle + second_angle)) < EPSILON)

EARTH_CENTER =  Point3D(0.0,0.0,0.0)
# This function returns a list of observations with the highest probabilities of containing 
# the trajectory of the observed fireball by ranking them based on the elevation values
def rank_observation_by_elevation_optimized(observation_list, list_highest_rank_points_by_elevation):
    list_ranked_observation = []
    list_point_of_observation_plane = []    
    point_map = {}    
    observation_plane_array = []
    start_vector_array = []
    end_vector_array = []
    observation_position_3d_array = []
    temp_observation_list = []
            
    for obs in observation_list:
        #obs = observation_list[index]
        point_of_specific_observation_plane = []
        
        
        start_vector = convert_angular_direction_vector_to_cartesian_vector(obs.get_latitude(), obs.get_longitude(), 
                                                                            obs.get_start_azimuth(), obs.get_start_elevation())        
        
        end_vector = convert_angular_direction_vector_to_cartesian_vector(obs.get_latitude(), obs.get_longitude(), 
                                                                         obs.get_end_azimuth(), obs.get_end_elevation())        
        
        observerPosition =  convert_lat_long_alt_to_cartesian_coordinates(obs.get_latitude(), obs.get_longitude(),0)
        #print(observerPosition)
        
                
        orthogonalVector = cross_product(start_vector,end_vector)
        
        A_component = orthogonalVector.get_X()
        B_component = orthogonalVector.get_Y()
        C_component = orthogonalVector.get_Z()
        D_component = (-1) * (orthogonalVector.get_X() * observerPosition.get_X() 
                              + orthogonalVector.get_Y() * observerPosition.get_Y() 
                              + orthogonalVector.get_Z() * observerPosition.get_Z() )
        observationPlane = Plane3D(A_component, B_component, C_component, D_component)
        
        if ( (math.pow(A_component,2) + math.pow(B_component,2) + math.pow(C_component,2)) != 0 ):            
            observation_plane_array.append(observationPlane)
            temp_observation_list.append(obs)
            start_vector_array.append(start_vector)
            end_vector_array.append(end_vector)
            observation_position_3d_array.append(observerPosition)
            list_point_of_observation_plane.append(point_of_specific_observation_plane)           
            print("observation plane is ",observationPlane.get_A(), observationPlane.get_B(), 
                                      observationPlane.get_C(), observationPlane.get_D() )
        else:
            print("This obs", obs.get_observation_ID(), "cannot form an observation plane")
    
    print("Number of observation plane is", len(observation_plane_array), "length of observation_list", len(temp_observation_list))    
    
    countPoint3D = 0
    pointNearTriangle = 0
    for i in range(X_MIN, X_MAX, GRID_SPACING_3D):
        print("i is now", i, " i max is ", X_MAX)
        for j in range(Y_MIN, Y_MAX, GRID_SPACING_3D):
            for k in range(Z_MIN, Z_MAX + MAXIMUM_FIREBALL_ALTITUDE, GRID_SPACING_3D):
                
                point = Point3D(i,j,k)
                if (compute_distance_between_points(EARTH_CENTER, point) >= EARTH_RADIUS + MINIMUM_FIREBALL_ALTITUDE):
                    countPoint3D += 1    
                    for index in range (0, len(observation_plane_array)):                         
                        if (compute_distance_point_to_plane(point, observation_plane_array[index]) < SEARCH_RADIUS_ELEVATION):
                            if (determine_projected_point_in_observation_triangle(point, start_vector_array[index], 
                                                                                  end_vector_array[index], 
                                                                                  observation_position_3d_array[index]) == True):                                
                                if (point_map.get(point) is None ):
                                    point_map[point] = 1
                                else:
                                    value = point_map.get(point)
                                    point_map[point] = value + 1
                                list_point_of_observation_plane[index].append(point)
                                pointNearTriangle += 1
    
    print("Number of point 3D above earth surface is ", countPoint3D)
    print("Number of point near triangle is", pointNearTriangle)                        
    
    map_observation_rank = {}
    for index in range(0, len(observation_plane_array)):
        max_rank_of_observation = 0
        if (len(list_point_of_observation_plane[index]) > 0):
            for i in range(0, len(list_point_of_observation_plane[index])):
                if (point_map.get(list_point_of_observation_plane[index][i]) > max_rank_of_observation):
                    max_rank_of_observation = point_map.get(list_point_of_observation_plane[index][i])
        map_observation_rank[temp_observation_list[index]] = max_rank_of_observation
    
    print("point_map 3D has size of ", len(point_map))
    
    max_rank_overall = max(point_map.values())
    second_max_rank_overall = find_second_max(list(point_map.values()))
    
    print("ELEVATION: Highest ranking is ", max_rank_overall, "and second highest is ", second_max_rank_overall)
    
    for obs in temp_observation_list:
        
        if ( (map_observation_rank.get(obs) == max_rank_overall) 
              or (map_observation_rank.get(obs) == second_max_rank_overall)):
            list_ranked_observation.append(obs)
    
    for temp_point in point_map.keys(): 
        if (point_map.get(temp_point) == max_rank_overall):
            list_highest_rank_points_by_elevation.append(temp_point)
    print("Number of 3D points with highest elevation ranking is", len(list_highest_rank_points_by_elevation))
            
    return list_ranked_observation

def convert_cartesian_coordintes_to_lat_long_alt(x,y,z):
    p = math.sqrt(math.pow(x, 2) + math.pow(y, 2))
    theta = math.atan((z * a) / (p * b));

    lon = math.atan(y / x);

    lat = math.atan(((z + math.pow(e2, 2) * b * math.pow(math.sin(theta), 3)) 
                     / ((p - math.pow(e, 2) * a * math.pow(math.cos(theta), 3)))));
    
    N = a / (math.sqrt(1 - (math.pow(e, 2) * math.pow(math.sin(lat), 2))));

    m = (p / math.cos(lat));
    height = m - N;
    lon = lon * 180/math.pi
    lat = lat * 180/math.pi
    
    result = [lat, lon, height]
    return result                                
                            
# This function determines a point on the trajectory line, which serves 
# as the basis for further determination of the trajectory.                           
def determine_point_on_line(list_points_ranked_by_azimuth, list_points_ranked_by_elevation):
    point_array = []
    selected_point = None
    for point in list_points_ranked_by_azimuth:
        point_array.append(point)
    
    multipoint = MultiPoint(point_array)
    convex_hull = multipoint.convex_hull
    centroid  = multipoint.centroid
    
    print("Centroid of the convex hull is", centroid)
    
    min_distance = sys.float_info.max
    
    geod = Geod(ellps="WGS84")
    for point1 in list_points_ranked_by_azimuth:
        for point2 in list_points_ranked_by_elevation:
            array = convert_cartesian_coordintes_to_lat_long_alt(point2.get_X()*1000, 
                                                                 point2.get_Y()*1000, 
                                                                 point2.get_Z()*1000 )
            point2_projected_earth = Point(array[1], array[0])
            
            distance = geod.inv(point1.x, point1.y, point2_projected_earth.x, point2_projected_earth.y)[2]
            if (distance < min_distance):
                min_distance = distance
                selected_point = point2
    
    altitude_of_selected_point = convert_cartesian_coordintes_to_lat_long_alt(selected_point.get_X()*1000, 
                                                                              selected_point.get_Y()*1000, 
                                                                              selected_point.get_Z()*1000)[2]
    print("Altitude of selected point is",altitude_of_selected_point)
    
    return convert_lat_long_alt_to_cartesian_coordinates(centroid.y, centroid.x, altitude_of_selected_point)
            
# This function determines the trajectory's start and end point on the trajectory line            
def determine_start_point_end_point_improved(point_on_traj_line, direction_vector, observation_ranked_azimuth):
    result = []
    geod = Geod(ellps='WGS84')
    
    point_collection = []
    
    for i in range(-2000,2001):
        point_XYZ = Point3D(point_on_traj_line.get_X() + direction_vector.get_X()*i, 
                            point_on_traj_line.get_Y() + direction_vector.get_Y()*i, 
                            point_on_traj_line.get_Z() + direction_vector.get_Z()*i)
                            
        if (compute_distance_between_points(point_XYZ, point_on_traj_line) < MAX_TRAJECTORY_LENGTH):
            conversion_result = convert_cartesian_coordintes_to_lat_long_alt(point_XYZ.get_X()*1000, 
                                                                             point_XYZ.get_Y()*1000, 
                                                                             point_XYZ.get_Z()*1000)
            projected_point = Point(conversion_result[1], conversion_result[0])
            point_collection.append(projected_point)
    
    sum_distance = 0
    for i in range(1, len(point_collection)):
        distance = geod.inv(point_collection[i-1].x, point_collection[i-1].y, 
                            point_collection[i].x, point_collection[i].y)[2] 
        sum_distance += distance
    
    print("First point is", point_collection[0].x, point_collection[0].y)
    print("Last point is", point_collection[len(point_collection)-1].x, point_collection[len(point_collection)-1].y)
    print("Distance of the projected whole traj line is", sum_distance)    
    
    long_trajectory_line = LineString([point_collection[0], point_collection[len(point_collection)-1]])
    
    point_collection_start_line = []
    point_collection_end_line = []
    
    for observation in observation_ranked_azimuth:
        observer_on_map = Point(observation.get_longitude(), observation.get_latitude())
        
        move_along_start_azimuth = move_geodetically(observer_on_map, observation.get_start_azimuth(), OBSERVATION_RANGE)
        start_viewing_line = LineString([observer_on_map, move_along_start_azimuth])
        
        move_along_end_azimuth = move_geodetically(observer_on_map, observation.get_end_azimuth(), OBSERVATION_RANGE)
        end_viewing_line = LineString([observer_on_map, move_along_end_azimuth])
        
        intersection_start = long_trajectory_line.intersection(start_viewing_line)
        
        if (intersection_start.is_empty is False):
            lon, lat = intersection_start.x, intersection_start.y
            point_collection_start_line.append(Point(lon, lat))
        
        intersection_end = long_trajectory_line.intersection(end_viewing_line)
        
        if (intersection_end.is_empty is False):
            lon, lat = intersection_end.x, intersection_end.y
            point_collection_end_line.append(Point(lon, lat))
    
    multipoint_start = MultiPoint(point_collection_start_line)    
    trajectory_start_point = multipoint_start.centroid  
    print("\n Final result -----------***-----------")
    print("START is", trajectory_start_point.y, trajectory_start_point.x)
    altitude_start = determine_altitude(trajectory_start_point, point_on_traj_line, direction_vector, 1000)
    print("Estimated start altitude is", altitude_start)
    
    multipoint_end = MultiPoint(point_collection_end_line)    
    trajectory_end_point = multipoint_end.centroid    
    altitude_end = determine_altitude(trajectory_end_point, point_on_traj_line, direction_vector, 1000)
    print("END is", trajectory_end_point.y, trajectory_end_point.x)
    print("Estimated end altitude is", altitude_end)
    
    return [trajectory_start_point.y, trajectory_start_point.x, altitude_start, 
            trajectory_end_point.y, trajectory_end_point.x, altitude_end]
        
# This function determines the altitude of a point on the trajectory line whose Earth projection is known. 
# The trajectory line is defined by a point and a direction vector   
def determine_altitude(point_on_map, point_on_traj_line, direction_vector, parameter):    
    geod = Geod(ellps='WGS84')
    
    first_point = Point3D(point_on_traj_line.get_X() + direction_vector.get_X()*(-1)*parameter, 
                          point_on_traj_line.get_Y() + direction_vector.get_Y()*(-1)*parameter, 
                          point_on_traj_line.get_Z() + direction_vector.get_Z()*(-1)*parameter )
    conversion_result = convert_cartesian_coordintes_to_lat_long_alt(first_point.get_X()*1000, 
                                                                     first_point.get_Y()*1000, 
                                                                     first_point.get_Z()*1000)
    first_point_on_map = Point(conversion_result[1], conversion_result[0])
    
    last_distance = geod.inv(point_on_map.x, point_on_map.y, first_point_on_map.x, first_point_on_map.y)[2]
    
    nearest_point_3D = None
    
    for i in range((-1)*parameter+1, parameter+1):        
        moving_point_on_traj_line = Point3D(point_on_traj_line.get_X() + direction_vector.get_X()*i, 
                                            point_on_traj_line.get_Y() + direction_vector.get_Y()*i, 
                                            point_on_traj_line.get_Z() + direction_vector.get_Z()*i )
                                            
        conversion_result = convert_cartesian_coordintes_to_lat_long_alt(moving_point_on_traj_line.get_X()*1000, 
                                                                         moving_point_on_traj_line.get_Y()*1000, 
                                                                         moving_point_on_traj_line.get_Z()*1000)
        moving_point_on_map = Point(conversion_result[1], conversion_result[0])
        
        temp_distance = geod.inv(point_on_map.x, point_on_map.y, 
                                 moving_point_on_map.x, moving_point_on_map.y)[2]        
        
        if (temp_distance < last_distance):
            last_distance = temp_distance
            nearest_point_3D = moving_point_on_traj_line
        else:
            break
    
    altitude = convert_cartesian_coordintes_to_lat_long_alt(nearest_point_3D.get_X()*1000, 
                                                            nearest_point_3D.get_Y()*1000, 
                                                            nearest_point_3D.get_Z()*1000)[2]    
    return altitude        
                
#path = "C:\\Users\\ngo_ma\\Documents\\Nachtlicht-BüHNE2\\Khoi\\IMO France - FRIPON reconstruction\\"
#eventID = input("Type event ID: ")

parser = argparse.ArgumentParser(description='Reconstruct the trajectory of a fireball based on its observations')
parser.add_argument('event_ID', type=str, help='The event ID is required')

args = parser.parse_args()

list_all_observation = []
list_highest_rank_points_by_azimuth = []
list_highest_rank_points_by_elevation = []

with open(args.event_ID + ".csv", mode = "r") as file:
    csv_file = csv.reader(file)
    count_line = 0;
    for lines in csv_file:
        if (count_line != 0):  # to skip the header line in the csv file            
            print(lines)
            if ((float(lines[5]) != float(lines[7])) or (float(lines[4]) != float(lines[6]))): # to avoid errorneous observations (observation that are 0° or 180°)
                obs = Observation(lines[0], lines[1], float(lines[3]), float(lines[2]), 
                                  float(lines[4]), float(lines[5]), float(lines[6]), float(lines[7]))
                list_all_observation.append(obs)            
        count_line += 1 
        
    print("Number of all observations from file is ", count_line-1, "number of observations as input for the reconstruction is ", len(list_all_observation))

reliable_obs_azimuth = rank_observation_by_azimuth_optimized(list_all_observation, OBSERVATION_RANGE, 
                                                             list_highest_rank_points_by_azimuth)    
print(len(reliable_obs_azimuth), "reliable observations by azimuth:")
for observation in reliable_obs_azimuth:
    print(observation.get_observation_ID())

reliable_obs_elevation = rank_observation_by_elevation_optimized(list_all_observation,
                                                                 list_highest_rank_points_by_elevation)
print(len(reliable_obs_elevation), " reliable observation by elevation:")
for observation in reliable_obs_elevation:
    print(observation.get_observation_ID())
    
#Combine the list of observation ranked by azimuth with the list of observation ranked by elevation
reliable_observation_list = set(reliable_obs_azimuth).intersection(set(reliable_obs_elevation))

print(len(reliable_observation_list)," reliable observation (azi + ele):")
for observation in reliable_observation_list:
    print(observation.get_observation_ID())

ai = []
bi = []
ci = []
Xc = []
Yc = []
Zc = []
lat = []
lon = []

string = ""

for observation in reliable_observation_list:
    observer_3D_location = convert_lat_long_alt_to_cartesian_coordinates(observation.get_latitude(), observation.get_longitude(),0)
    
    start_vector = convert_angular_direction_vector_to_cartesian_vector(observation.get_latitude(), observation.get_longitude(), 
                                                                        observation.get_start_azimuth(), observation.get_start_elevation())
                                                                        
    end_vector = convert_angular_direction_vector_to_cartesian_vector(observation.get_latitude(), observation.get_longitude(), 
                                                                      observation.get_end_azimuth(),observation.get_end_elevation())
    
    string += str(start_vector.get_X()) + " " + str(start_vector.get_Y()) + " " + str(start_vector.get_Z()) + " " + str(observer_3D_location.get_X()) + " " + str(observer_3D_location.get_Y()) + " " + str(observer_3D_location.get_Z()) + " " + str(observation.get_latitude()) + " " + str(observation.get_longitude()) + "\n"    
    
    string += str(end_vector.get_X()) + " " + str(end_vector.get_Y()) + " " + str(end_vector.get_Z()) + " " + " " + str(observer_3D_location.get_X()) + " " + str(observer_3D_location.get_Y()) + " " + str(observer_3D_location.get_Z()) + " " + str(observation.get_latitude()) + " " + str(observation.get_longitude()) + "\n"

#Produce the input for the coplanarity method   
input_string = StringIO(string)

# read the observations
ai, bi, ci, Xc, Yc, Zc, lat, lon = np.loadtxt(input_string, skiprows=0, dtype=float, unpack=True)
nL = len(ai)
print(ai)

# Observations and stochastic parameters
# ----------------------------------------------------------------------------------------
nDir = 2                # number of directions from each station
assert nL % nDir == 0, "there should number of measurements in the input file should be dividable by nDir"
nStations = int(nL / nDir)          # number of stations
nX = 4                  # number of unknowns (Xm,Ym,Zm,b,c)  

# standard deviations of the observations
# equal weights to all stations/observers
std = np.full(nStations, 1)

# Compute starting values 
# ---------------------------------------------------------------------------------------
# Parameters am and Xm are taken as fixed in the adjustment
# Input Data from first 2 Stations
a1, b1, c1, Xc1, Yc1, Zc1  = ai[0:nDir], bi[0:nDir], ci[0:nDir], Xc[0:nDir], Yc[0:nDir], Zc[0:nDir]
a2, b2, c2, Xc2, Yc2, Zc2  = ai[nDir:nDir*2], bi[nDir:nDir*2], ci[nDir:nDir*2], Xc[nDir:nDir*2], Yc[nDir:nDir*2], Zc[nDir:nDir*2]
am,bm,cm,Xm,Ym,Zm = Line_intersection.meteor_position(a1, b1, c1, Xc1, Yc1, Zc1, a2, b2, c2, Xc2, Yc2, Zc2)

# normalise to am
bm=bm/am
cm=cm/am
am=1

# Starting valued for the errors
va = np.zeros((nL))
vb = np.zeros((nL))
vc = np.zeros((nL))
v_matrix = np.mat(np.hstack([va,vb,vc]))

# this is just for the Qll matrix

std_ai, std_bi, std_ci, Qll = Functions.calcErrorSphere2Rect(nDir,nStations,ai,bi,ci,std,std)
Qll = np.identity(nL*3)     # 3 --> x,y,z


# Functional model and Design matrices
# ---------------------------------------------------------------------------------------

# Determinant of the direction vectors from a station and the "meteor" direction vector should be zero [Eq.23]
def functional_model(XmS,YmS,ZmS, XcS,YcS,ZcS, aiS,biS,ciS, vaS,vbS,vcS, amS,bmS,cmS):
    g = (XmS-XcS)*(bmS*(ciS-vcS)-cmS*(biS-vbS)) - (YmS-YcS)*(amS*(ciS-vcS)-cmS*(aiS-vaS)) + \
        (ZmS-ZcS)*(amS*(biS-vbS)-bmS*(aiS-vaS))
    return g

B = np.zeros((nL,nL*3))
A = np.zeros((nL,4))

def Design_matrix(Xm,Ym,Zm, Xc,Yc,Zc, ai,bi,ci, va,vb,vc, am,bm,cm):
    B_va = np.mat(np.diag( bm*(Zm-Zc)-cm*(Ym-Yc) )) 
    B_vb = np.mat(np.diag( cm*(Xm-Xc)-am*(Zm-Zc) )) 
    B_vc = np.mat(np.diag( am*(Ym-Yc)-bm*(Xm-Xc) )) 
    B = np.bmat([B_va, B_vb, B_vc])
    
    
    A_Ym = np.mat([cm*(ai-va)-am*(ci-vc)]).T
    A_Zm = np.mat([am*(bi-vb)-bm*(ai-va)]).T
    A_bm = np.mat([(Xm-Xc)*(ci-vc)-(Zm-Zc)*(ai-va)]).T
    A_cm = np.mat([(Ym-Yc)*(ai-va)-(Xm-Xc)*(bi-vb)]).T
    A = np.bmat([A_Ym, A_Zm, A_bm, A_cm])
    return B, A

# Iterative procedure for the adjustment in the GH model
# ---------------------------------------------------------------------------------------
stopping_value = 1
iteration = 0
        
while stopping_value > 1e-10:
    
    # compute design matrices
    B, A = Design_matrix(Xm,Ym,Zm, Xc,Yc,Zc, ai,bi,ci, va,vb,vc, am,bm,cm)
        
    # vector of misclosures w
    g = functional_model(Xm,Ym,Zm, Xc,Yc,Zc, ai,bi,ci, va,vb,vc, am,bm,cm)
    w = g - np.transpose(B*v_matrix.T)
    w = w.T
    
    # normal matrix N
    N = np.bmat([[B*Qll*B.T, A], [A.T, np.zeros((4,4))] ])
    
    # vector of absolute values n 
    n = np.vstack((-w,np.mat(np.zeros((nX,1)))))
    
    Dx = np.linalg.solve(N,n) 
    x_hat = Dx[nL:]
    
    # Langrange multiplier k
    k = Dx[:nL]
    
    vv = Qll * B.T * k
    
    # update approximated values for the unknown parameters
    Ym = float(Ym + x_hat[0])
    Zm = float(Zm + x_hat[1])
    bm = float(bm + x_hat[2])
    cm = float(cm + x_hat[3])
    
    # update approximated values for the residuals
    va = np.array(vv[0:nL]).flatten()           # [array-name].flatten() --> 2D to 1D array
    vb = np.array(vv[nL:nL*2]).flatten()
    vc = np.array(vv[nL*2:]).flatten()
    
    v_matrix = np.mat(np.hstack([va,vb,vc]))
   
    iteration = iteration + 1 
    stopping_value = max(np.abs(x_hat))
    #print(stopping_value)

# direction of the meteor
abc = [am,bm,cm]/np.sqrt(am**2+bm**2+cm**2)
print("abc vector is", abc)

direction_vector = Vector3D(abc[0], abc[1], abc[2])

point_on_traj_line_XYZ = determine_point_on_line(list_highest_rank_points_by_azimuth, list_highest_rank_points_by_elevation)

point_on_traj_line_LLA = convert_cartesian_coordintes_to_lat_long_alt(point_on_traj_line_XYZ.get_X()*1000, 
                                                                      point_on_traj_line_XYZ.get_Y()*1000, 
                                                                      point_on_traj_line_XYZ.get_Z()*1000)

print("point on traj line LLA:", point_on_traj_line_LLA[0], point_on_traj_line_LLA[1], point_on_traj_line_LLA[2], "or XYZ", point_on_traj_line_XYZ.get_X(), point_on_traj_line_XYZ.get_Y(), point_on_traj_line_XYZ.get_Z())

determine_start_point_end_point_improved(point_on_traj_line_XYZ, direction_vector, reliable_obs_azimuth)

# adjusted measurements
ai_adj = ai - va
bi_adj = bi - vb
ci_adj = ci - vc

