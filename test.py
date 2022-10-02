from shapely.geometry import Polygon, Point, LineString
import geopandas as gpd
import shape_utility as su

def get_exterior_boundary_from_polygon(polygon: Polygon):
    #this function returns all edges of a polygon as shapely LineString in a list
    coords = polygon.exterior.coords
    total = len(coords)
    lines = []
    for i in range(total - 1):
        lines.append(LineString(coords[i:i+2]))
    return lines

def is_intersect_or_contain(p_src: Polygon, p_tar: Polygon):
    for point in p_tar.exterior.coords:
        if p_src.contains(Point(point)):
            return True
    for edge in get_exterior_boundary_from_polygon(p_tar):
        if edge.intersects(p_src):
            return True
    return False

def scan_polygon_to_grid(polygon: Polygon, interval: float):
    #convert polygon into a bit map.
    #the resolution is interval
    #a pixel represents a square polygon whose side length is equal to interval
    #a pixel is 1 when it is intersected with the polygon or any vertex of the polygon is inside the pixel, otherwise it is 0
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
        x_next = x_curr+interval
    return grid

def z_order_curve(array, x_off: int, y_off: int, degree: int):
    #following the pattern of Z-order curve, convert grid into a of string of 0s and 1s
    if degree == 1:
        if len(array) > y_off and len((array[y_off])) > x_off:
            return str(array[y_off][x_off])
        else:
            return '0'
    else:
        quadrant = []
        move = degree // 2
        quadrant.append(z_order_curve(array, x_off, y_off, move))    #quadrant (0, 0)
        quadrant.append(z_order_curve(array, x_off + move, y_off, move))    #quadrant (0, 1)
        quadrant.append(z_order_curve(array, x_off, y_off + move, move))    #quadrant (1, 0)
        quadrant.append(z_order_curve(array, x_off + move, y_off + move, move)) #quadrant (1, 1)
        return ''.join(quadrant)

def grid_to_list(polygon: Polygon, interval: float, model: int):
    #convert grid into a string of 0s and 1s
    #model = 1: sinple append line by line
    #model = 2: use Z-order curve to keep more locality information
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


geoshape = gpd.read_file("map.shp")
shape = geoshape.iloc[0]['geometry']
grid = (su.scan_polygon_to_grid(shape, 0.05))

geomatrics = []
for i in range(len(geoshape)):
    g = geoshape.iloc[0]['geometry']
    if type(g) == Polygon:
        geomatrics.append(g)

#disks = su.generate_multi_resolution_mask(geomatrics, 24, 2)
#hash_s = su.get_multi_resolution_hash(disks, shape, 2)
#print(hash_s)

shape2 = geoshape.iloc[1]['geometry']
print(shape.hausdorff_distance(shape2))
print(su.get_hausdorff_distance(shape, shape2))

p = Point(1, 1)
print(p.coords[0])