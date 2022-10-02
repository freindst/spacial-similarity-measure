from shapely.geometry import Point, MultiPoint, LineString, Polygon
import random
import similaritymeasures as sm
import numpy as np


def get_exterior_edge_from_polygon(polygon: Polygon):
    # this function returns all edges of a polygon as shapely LineString in a list
    coords = polygon.exterior.coords
    total = len(coords)
    lines = []
    for i in range(total - 1):
        lines.append(LineString(coords[i:i + 2]))
    return lines


def get_centralized_edge_from_polygon(polygon: Polygon):
    # this function returns all edges of a polygon as shapely LineString in a list shift to its centroid
    centroid = polygon.centroid
    coords = polygon.exterior.coords
    total = len(coords)
    lines = []
    for i in range(total - 1):
        [(x0, y0), (x1, y1)] = coords[i:i + 2]
        lines.append(LineString([(x0 - centroid.x, y0 - centroid.y), (x1 - centroid.x, y1 - centroid.y)]))
    return lines


# ********** polygon to grid **********
def scan_polygon_to_grid(polygon: Polygon, interval: float):
    # convert polygon into a bit map.
    # the resolution is interval
    # a pixel represents a square polygon whose side length is equal to interval
    # a pixel is 1 when it is intersected with the polygon or any vertex of the polygon is inside the pixel, otherwise it is 0
    bounds = polygon.bounds
    x0, y0, x1, y1 = bounds[0], bounds[1], bounds[2], bounds[3]
    grid = []
    x_curr = x0
    y_curr = y0
    x_next = x_curr + interval
    y_next = y_curr + interval
    while y_next < y1:
        row = []
        while x_next < x1:
            pixel_curr = Polygon([
                (x_curr, y_curr),
                (x_curr, y_next),
                (x_next, y_next),
                (x_next, y_curr)
            ])
            if is_intersect_or_contain(pixel_curr, polygon):
                row.append(1)
            else:
                row.append(0)
            x_curr = x_next
            x_next += interval
        grid.append(row)
        y_curr = y_next
        y_next += interval
        x_curr = x0
        x_next = x_curr + interval
    return grid


def is_intersect_or_contain(p_src: Polygon, p_tar: Polygon):
    # check if p_src intersects with p_tar or any point of p_tar is inside p_src
    for point in p_tar.exterior.coords:
        if p_src.contains(Point(point)):
            return True
    for edge in get_exterior_edge_from_polygon(p_tar):
        if edge.intersects(p_src):
            return True;
    return False


def z_order_curve(array, x_off: int, y_off: int, degree: int):
    # following the pattern of Z-order curve, convert grid into a of string of 0s and 1s
    if degree == 1:
        if len(array) > y_off and len((array[y_off])) > x_off:
            return str(array[y_off][x_off])
        else:
            return '0'
    else:
        quadrant = []
        move = degree // 2
        quadrant.append(z_order_curve(array, x_off, y_off, move))  # quadrant (0, 0)
        quadrant.append(z_order_curve(array, x_off + move, y_off, move))  # quadrant (0, 1)
        quadrant.append(z_order_curve(array, x_off, y_off + move, move))  # quadrant (1, 0)
        quadrant.append(z_order_curve(array, x_off + move, y_off + move, move))  # quadrant (1, 1)
        return ''.join(quadrant)


def grid_to_string(polygon: Polygon, interval: float, model: int):
    # convert grid into a string of 0s and 1s
    # model = 1: sinple append line by line
    # model = 2: use Z-order curve to keep more locality information
    grid = scan_polygon_to_grid(polygon, interval)
    if model == 1:
        result = ""
        for row in grid:
            for cell in row:
                result += str(cell)
        return result
    else:
        degree = 1
        y_max = len(grid)
        x_max = 0
        if y_max > 0:
            x_max = len(grid[0])
        while degree < y_max or degree < x_max:
            degree *= 2
        return z_order_curve(grid, 0, 0, degree)


# ********** polygon to multi-resolution hash **********
def get_canvas_from_polygons(polygons):
    # determine the space where the multi-resolution mask will be scattered
    x0, y0, x1, y1 = 0, 0, 0, 0
    for polygon in polygons:
        bounds = polygon.bounds
        centroid = polygon.centroid
        x0 = min(bounds[0] - centroid.x, x0)
        y0 = min(bounds[1] - centroid.y, y0)
        x1 = max(bounds[2] - centroid.x, x1)
        y1 = max(bounds[3] - centroid.y, y1)
    return {
        "x0": x0,
        "y0": y0,
        "x1": x1,
        "y1": y1
    }


def generate_multi_resolution_mask(polygons, number, radius):
    # generate the masks, and keep it consistent when generating hash from polygons
    disks = []
    canvas = get_canvas_from_polygons(polygons)
    x_span = canvas["x1"] - canvas["x0"]
    y_span = canvas["y1"] - canvas["y0"]
    for i in range(number):
        x = random.random() * (x_span + 2 * radius) + canvas["x0"]
        y = random.random() * (y_span + 2 * radius) + canvas["y0"]
        disks.append(Point(x, y))
    return disks


def get_multi_resolution_hash(disks, polygon, radius):
    # record id of each round mask the polygon interacts with, use as the hash value
    hash_str = ''
    edges = get_centralized_edge_from_polygon(polygon)  # shift the polygon coordinates to the centroid point
    for edge in edges:
        intersections = []
        for i in range(len(disks)):
            point = disks[i]
            disk = point.buffer(radius).boundary
            foundIntersect = edge.intersection(disk)
            if type(foundIntersect) == Point:  # a line only has one intersect with a circle
                distance = foundIntersect.distance(Point(edge.coords[0]))
                intersections.append({"id": str(i), "distance": distance})
            elif type(foundIntersect) == MultiPoint:  # a line has two intersect with a circle
                for p in foundIntersect.geoms:
                    distance = p.distance(Point(edge.coords[0]))
                    intersections.append({"id": str(i), "distance": distance})
        intersections.sort(key=lambda x: x["distance"], reverse=False)
        hash_str += ''.join(p['id'] for p in intersections)
    return hash_str


#********** spatial distance **********
def get_hausdorff_distance(polygon1: Polygon, polygon2: Polygon):
    # use hausdorff distance calculation method built-in shapely package to get hausdorff distance between two polygon
    array1 = []
    centroid1 = polygon1.centroid
    for (x, y) in polygon1.exterior.coords:
        array1.append([x - centroid1.x, y - centroid1.y])
    p1 = Polygon(array1)

    array2 = []
    centroid2 = polygon2.centroid
    for (x, y) in polygon2.exterior.coords:
        array2.append([x - centroid2.x, y - centroid2.y])
    p2 = Polygon(array2)
    return p1.hausdorff_distance(p2)


def get_fretcher_distance(polygon1: Polygon, polygon2: Polygon):
    # use frechet distance method in similarity measures package to calculate frechet distance
    array1 = []
    centroid1 = polygon1.centroid
    for (x, y) in polygon1.exterior.coords:
        array1.append([x - centroid1.x, y - centroid1.y])

    array2 = []
    centroid2 = polygon2.centroid
    for (x, y) in polygon2.exterior.coords:
        array2.append([x - centroid2.x, y - centroid2.y])

    return sm.frechet_dist(np.array(array1), np.array(array2))