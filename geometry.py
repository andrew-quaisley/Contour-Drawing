import math

#-------- Basic Math Stuff ---------

#given an x-y coordinate for the Cartesian plane, converts it to window coordinates (top left is (0,0), y-coord increases down the screen)
#includes the option to scale vertically and horizontally, and to set where (in window coordinates) the origin goes.
def window_coord(coord, origin, x_scale, y_scale):
	x, y = coord
	x_0, y_0 = origin
	return (x*x_scale + x_0, -y*y_scale + y_0)


#given a point (center) and a list of other points, returns a list of the same points sorted into clockwise order around the given point.  Keeps first entry of the list the same.
#assumes points given in window coordinates, not Cartesian (otherwise, it goes counter-clockwise)
def sort_edges_clockwise(center, adj_list):
	basis_0 = (adj_list[0][0] - center[0], adj_list[0][1] - center[1])	#subtract center and first entry to get first basis vector
	length = math.sqrt(basis_0[0]**2 + basis_0[1]**2)
	basis_0 = (basis_0[0]/length, basis_0[1]/length)			#normalize so we're comparing vectors of the same length 
	
	basis_1 = (-basis_0[1], basis_0[0])				#find orthogonal second basis vector in clockwise direction (bizarre point: because window coordinates are flipped over the x-axis from Cartesian, need to make this basis vector look like it's going in counter-clockwise direction from Cartesian perspective)

	#print('basis_1 = ' + str(basis_1))
	
	sorted_list = [adj_list[0]]						#sorted list of points
	value_list = [(1, basis_0[0]**2 + basis_0[1]**2)]			#sorted list of values
	for p in adj_list[1:]:
		#print('considering ' + str(p))
		vector = (p[0] - center[0], p[1] - center[1])
		length = math.sqrt(vector[0]**2 + vector[1]**2)
		normal = (vector[0]/length, vector[1]/length)		#normalize so we're comparing vectors of the same length
		
		if basis_1[0]*normal[0] + basis_1[1]*normal[1] >= 0:		#if point lies above basis_0 (in direction of basis_1)
			vertical = 1
		else:
			vertical = -1
			
		value = (vertical, basis_0[0]*normal[0] + basis_0[1]*normal[1])	#value is first measured by side of basis_0 it lies on, then by inner product of the vector with basis_0.  Sorting these tuples will allow us to sort the vectors.
		#print('value for ' + str(p) + ': ' + str(value))
		
		for i in range(len(sorted_list)):
			if correct_order(value_list[i], value) and (i+1 == len(sorted_list) or correct_order(value, value_list[i+1])):	#if value fits between i and i+1 entries
				sorted_list.insert(i+1, p)		#put the point in that index of the list
				#print('inserting ' + str(p) + ' at index ' + str(i+1))
				value_list.insert(i+1, value)		#and put its value in that index of the value list
				break

	return sorted_list


#defines an order on the values computed in sort_edges_clockwise
#returns True if value_1 < value_2, False otherwise
def correct_order(value_1, value_2):
	if value_1[0] > value_2[0]:			#if first point above basis_0 and second point below
		return True
	elif value_1[0] < value_2[0]:			#other way around
		return False
	elif value_1[0] == 1:				#if they are both above
		return value_1[1] >= value_2[1]		#in order if first point closer to basis_0 than second point (and therefore if bigger inner product)
	elif value_1[0] == -1:				#if they are both below
		return value_1[1] <= value_2[1]		#in order if first point further from basis_0 than second point
	else:
		raise ValueError()


#x, y, and z are points in the plane
#returns the angle from the ray yx to the ray yz, travelled in the clockwise direction
#assumes points are in window coordinates, not Cartesian (otherwise, it measures counter-clockwise)
def clockwise_angle(x, y, z):
	basis_0 = (x[0] - y[0], x[1] - y[1])				#subtract center and first entry to get first basis vector
	ray = (z[0] - y[0], z[1] - y[1])

	length_to_x = math.sqrt(basis_0[0]**2 + basis_0[1]**2)		#normalize so everything is on unit circle
	basis_0 = (basis_0[0]/length_to_x, basis_0[1]/length_to_x)

	length_to_z = math.sqrt(ray[0]**2 + ray[1]**2)
	ray = (ray[0]/length_to_z, ray[1]/length_to_z)
	
	basis_1 = (-basis_0[1], basis_0[0])				#find orthogonal second basis vector in clockwise direction (bizarre point: because window coordinates are flipped over the x-axis from Cartesian, need to make this basis vector look like it's going in counter-clockwise direction from Cartesian perspective)

	cosine_of_angle = basis_0[0]*ray[0] + basis_0[1]*ray[1]		#cosine of angle is the inner product of the two vectors (once they have length 1)
	if cosine_of_angle > 1:						#have encountered rounding errors making this number slightly outside of range of cosine
		if cosine_of_angle - 1 < 0.0001:
			cosine_of_angle = 1
	elif cosine_of_angle < -1:
		if -1 - cosine_of_angle < 0.0001:
			cosine_of_angle = -1

	angle = math.acos(cosine_of_angle)
	if basis_1[0]*ray[0] + basis_1[1]*ray[1] < 0:			#if ray is below basis_0 (from the perspective of basis_1), and thus in quad 3 or 4
		angle = 2*math.pi - angle				#appropriate transformation to get angle greater than pi

	return angle


#returns coordinates of the point such that the angle in radians between the ray center-p and the ray center-point is angle, measuring clockwise, and its distance from center is the same as that of p.
#everything in window coordinates
def point_at_angle(p, center, angle):
	#print('finding point at clockwise angle ' + str(angle) + ' from the ray from ' + str(center) + ' to ' + str(p))
	b0 = (p[0] - center[0], p[1] - center[1])
	b1 = (-b0[1], b0[0])						#perpendicular to b0 in the clockwise direction (for window coordinates)
	#print('ray from ' + str(center) + ' to ' + str(p) + ', as a vector: ' + str(b0))
	#print('ray in clockwise direction perpendicular to that one, as a vector: ' + str(b1))

	#vector in correct direction: v = b0cos(angle) + b1sin(angle)
	v = (b0[0]*math.cos(angle) + b1[0]*math.sin(angle), b0[1]*math.cos(angle) + b1[1]*math.sin(angle))
	#print('ray in correct dirction, as a vector: ' + str(v))

	return (center[0] + v[0], center[1] + v[1])


#returns True if point p is on the left side of the line from ray1 to ray2; returns False otherwise (including if p is on the line)
#note that everything is in window coordinates, which are reflected from Cartesian
def on_left_side(p, ray1, ray2):
	b1 = (ray2[1] - ray1[1], ray1[0] - ray2[0])
	pvec = (p[0] - ray1[0], p[1] - ray1[1])
	return b1[0]*pvec[0] + b1[1]*pvec[1] > 0


#returns True if point p is on the right side of the line from ray1 to ray2; returns False otherwise (including if p is on the line)
#note that everything is in window coordinates, which are reflected from Cartesian
def on_right_side(p, ray1, ray2):
	b1 = (ray1[1] - ray2[1], ray2[0] - ray1[0])
	pvec = (p[0] - ray1[0], p[1] - ray1[1])
	return b1[0]*pvec[0] + b1[1]*pvec[1] > 0


#distance between x-y coords p and q in the plane
def distance(p, q):
	p1, p2 = p
	q1, q2 = q
	return math.sqrt((p1-q1)**2 + (p2 - q2)**2)


#finds the point of intersection of line connecting w and x and line connecting y and z
#if they are the same line, returns w
#if they do not intersect, returns False
def intersection(w, x, y, z):
	w1, w2 = w
	x1, x2 = x
	y1, y2 = y
	z1, z2 = z
	if w1 == x1:
		if y1 == z1:
			if w1 != y1:
				return False
			else:
				return w
		else:
			m2 = (y2-z2)/(y1-z1)
			b2 = z2 - z1*m2
			return (w1, m2*w1 + b2)
	elif y1 == z1:
		m1 = (w2-x2)/(w1-x1)
		b1 = x2 - x1*m1
		return (y1, m1*y1 + b1)
	else:
		m1 = (w2-x2)/(w1-x1)
		b1 = x2 - x1*m1
		m2 = (y2-z2)/(y1-z1)
		b2 = z2 - z1*m2
		if m1 == m2:
			if b1 == b2:
				return w
			else:
				return False
		else:
			p1 = (b2-b1)/(m1-m2)
			return (p1, m1*p1 + b1)


#returns vector from point a to point b
def vector(a, b):
	a0, a1 = a
	b0, b1 = b
	return (b0-a0, b1-a1)


#returns normal vector in direction of the vector v
def normal(v):
	v0, v1 = v
	length = math.sqrt(v0**2 + v1**2)
	return ( v0/length, v1/length )


#returns the sum of the two input vectors
def add(v,w):
	v0, v1 = v
	w0, w1 = w
	return (v0 + w0, v1 + w1)

#returns a*v
#a is a number, v is a vector
def mult(a, v):
	v0, v1 = v
	return (a*v0, a*v1)


#returns the dot product of vectors v and w
def dot(v, w):
	v0, v1 = v
	w0, w1 = w
	return v0*w0 + v1*w1