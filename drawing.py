import geometry as g

#draws vKd given by the vertex dictionary on the given canvas
#adds, to each vertex's dictionary in vert_dict, a key 'id', giving the canvas id of the circle drawn at that vertex, a key 'dist id', giving the canvas id of the distance label for that vertex, and a key 'edge ids', giving a dictionary with a key for each adjacent vertex with corresponding value the canvas id of the edge drawn between the two vertices
#gives each distance label the tag 'distance', each vertex the tag 'vertex', and each edge the tag 'edge'
#vertex dictionary contains a key for each vertex, and the corresponding value is a dictionary with keys 'coords' and 'adj', corresponding to the window coordinates of the vertex in the square [0,1]x[0,1].
#vertex dictionary may or may not contain 'distance' keys for each vertex
def draw_vKd(vert_dict, canvas, location, width, height):
	vert_list = list(vert_dict)
	x_0, y_0 = location
	canvas_width = canvas.winfo_width()
	canvas_height = canvas.winfo_height()

	for i in range(len(vert_list)):
		v = vert_list[i]
		x_1 = vert_dict[v]['coords'][0]*width + x_0
		y_1 = vert_dict[v]['coords'][1]*height + y_0
		radius = 1.5					#will be updated to be smaller if vertices are too close together
		vert_dict[v]['dist id'] = canvas.create_text(x_1 + 6, y_1 + 6, text=str(vert_dict[v]['distance']), font=('TimesNewRoman 7'), tags='distance') 	#Labels each vertex with its distance.
		if 'edge ids' not in vert_dict[v]:	#initialize dictionary of edge ids
			vert_dict[v]['edge ids'] = {}
		for w in vert_dict[v]['adj']:
			x_2 = vert_dict[w]['coords'][0]*width + x_0
			y_2 = vert_dict[w]['coords'][1]*height + y_0
			radius = min(radius, g.distance((x_1, y_1), (x_2, y_2))/8)
			if vert_list.index(w) > i:
				edge_id = canvas.create_line(x_1, y_1, x_2, y_2, tags='edge')
				vert_dict[v]['edge ids'][w] = edge_id
				if 'edge ids' not in vert_dict[w]:
					vert_dict[w]['edge ids'] = {v:edge_id}
				else:
					vert_dict[w]['edge ids'][v] = edge_id
		
		v_id = canvas.create_oval(x_1 - radius, y_1 - radius, x_1 + radius, y_1 + radius, fill='blue', tags='vertex')
		vert_dict[v]['id'] = v_id
		canvas.tag_raise(v_id)	#put vertices in front



#given contour dict, draws contour lines (needs vert_dict for actual locations, and other parameters to draw in same spot as vKd)
#def draw_contour_lines(contour_dict, vert_dict, canvas, location, width, height):
#	x_0, y_0 = location
#	for cell in list(contour_dict):
#		for edge in contour_dict[cell]:
#			v_1, v_2 = edge
#			x_1, y_1 = vert_dict[v_1]['coords']
#			x_2, y_2 = vert_dict[v_2]['coords']
#			canvas.create_line(x_1*width + x_0, y_1*height + y_0, x_2*width + x_0, y_2*height + y_0, fill='red')



#draws contour lines using the contour segments given in segments (needs other parameters to draw in same spot as vKd).
#adds a key 'segment ids' to each cell's dictionary, with value a list of canvas ids of the segments drawn in that cell
#each segment drawn is given the tag 'contour segment'
#segments is a dictionary where key n corresponds to a list of segments that make up the contour lines of distance n away from the basepoint.
def draw_contour_lines_better(cell_dict, canvas, location, width, height, cells_to_draw = None):
	x_0, y_0 = location

	if cells_to_draw is None:
		cells_to_draw = list(cell_dict)

	for cell_index in cells_to_draw:
		if cell_index == 'boundary':
			continue
		cell_dict[cell_index]['segment ids'] = []
		segments = cell_dict[cell_index]['segments']
		for n in segments:
			for s in segments[n]:
				((p1, p2), (q1, q2)) = s
				if n%2 == 0:
					color = 'red'
				else:
					color = 'blue'
				cell_dict[cell_index]['segment ids'].append(canvas.create_line(x_0 + width*p1, y_0 + height*p2, x_0 + width*q1, y_0 + height*q2, fill=color, tags='segment'))


#draws the cuts in green for debugging purposes
#each segment drawn is given the tag 'cut'
def draw_cuts(cell_dict, canvas, location, width, height):
	x_0, y_0 = location
	for cell_index in cell_dict:
		if cell_index == 'boundary':	#if you actually have the boundary circuit, just move on to next cell_index
			continue

		for (v_i, v_j) in cell_dict[cell_index]['cut segments']:
			canvas.create_line(v_i[0]*width + x_0, v_i[1]*height + y_0, v_j[0]*width + x_0, v_j[1]*height + y_0, fill='green', tags='cut')



#determines if there is a vertex within a reasonable screen distance of the screen coordinates given, and if so, returns its index in vert_dict and its coordinates on the canvas as a tuple.  If not, returns False.
#"reasonable screen distance" means a distance within which a user might expect to be able to interact with the vertex if the cursor is at the coordinates given.
def closest_vertex(mouse_x, mouse_y, vert_dict, canvas, origin_on_canvas, x_scale, y_scale):
	canvas_mouse_x = canvas.canvasx(mouse_x)	#find where the mouse is in canvas coordinates
	canvas_mouse_y = canvas.canvasy(mouse_y)
	min_square_distance = 50	#start with the biggest (squared) distance from which you might select a vertex
	close_vertex_exists = False

	for v in vert_dict:
		diagram_x, diagram_y = vert_dict[v]['coords']
		canvas_x = origin_on_canvas[0] + x_scale * diagram_x
		square_difference_x = (canvas_x - canvas_mouse_x)**2
		if square_difference_x < min_square_distance:
			canvas_y = origin_on_canvas[1] + y_scale * diagram_y
			square_difference_y = (canvas_y - canvas_mouse_y)**2
			square_distance = square_difference_y + square_difference_x
			if square_distance < min_square_distance:
				close_vertex_exists = True
				closest_vertex = v
				min_square_distance = square_distance
				closest_vertex_coords = (canvas_x, canvas_y)

	if close_vertex_exists:
		return (closest_vertex, closest_vertex_coords)
	else:
		return False


#determines if there is an edge out of the vertex with index v within a reasonable screen distance of the screen coordinates given, and if so, returns the index of its other endpoint in vert_dict, the coordinates of v, and the coordinates of the other endpoint, as a tuple.  If not, returns False.
#"reasonable screen distance" means a distance within which a user might expect to be able to interact with the edge if the cursor is at the coordinates given.
def closest_edge_from(v, mouse_x, mouse_y, vert_dict, canvas, origin_on_canvas, x_scale, y_scale):
	canvas_mouse_x = canvas.canvasx(mouse_x)	#find where the mouse is in canvas coordinates
	canvas_mouse_y = canvas.canvasy(mouse_y)
	v_coords = vert_dict[v]['coords']		#find where v is in canvas coordinates
	canvas_v_coords = g.add(origin_on_canvas, (x_scale*v_coords[0], y_scale*v_coords[1]))
	min_squared_distance = 50	#start with the biggest (squared) distance from which you might select an edge
	close_edge_exists = False

	for w in vert_dict[v]['adj']:
		w_coords = vert_dict[w]['coords']
		canvas_w_coords = g.add(origin_on_canvas, (x_scale*w_coords[0], y_scale*w_coords[1]))
		w_to_mouse_vector = g.vector(canvas_w_coords, (canvas_mouse_x, canvas_mouse_y))
		v_to_w_vector = g.vector(canvas_v_coords, canvas_w_coords)
		signed_distance_along_edge = g.dot(w_to_mouse_vector, g.normal(v_to_w_vector))
		w_to_mouse_squared_distance = (canvas_w_coords[0]-canvas_mouse_x)**2 + (canvas_w_coords[1]-canvas_mouse_y)**2
		perpendicular_squared_distance = w_to_mouse_squared_distance - signed_distance_along_edge**2
		if perpendicular_squared_distance < min_squared_distance:
			v_to_w_distance = g.distance(canvas_v_coords, canvas_w_coords)
			if (-v_to_w_distance < signed_distance_along_edge < 0) or w_to_mouse_squared_distance <= 50:
				close_edge_exists = True
				closest_edge_endpoint = w
				closest_edge_endpoint_coords = canvas_w_coords
				min_squared_distance = perpendicular_squared_distance

	if close_edge_exists:
		return (closest_edge_endpoint, canvas_v_coords, closest_edge_endpoint_coords)
	else:
		return False
