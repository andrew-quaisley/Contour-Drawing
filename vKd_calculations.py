import geometry as g
import math

#Global constants:
max_angle = math.pi - 0.01	#largest angle (in radians) at which a cut is (just barely) necessary

#returns a tuple whose first value is a vertex dictionary for the 1-skeleton of the standard van Kampen diagram for the word w_k = t^kat^(-k)at^ka^(-1)t^(-k)a^(-1) in BS(1,n), and whose second value is the index of the basepoint.
#keys for vertex dictionary are indices of a vertex, given in the form (i,h), where vertex (i,h) is the ith vertex in row h, the rows ranging from -k to k.
#the value in the dictionary for vertex (i,h) is a dictionary, with keys 'coords', 'adj', and 'distance' corresponding to the window coordinates (origin at top left of the window, y-values going top-to-bottom) of the vertex, a list of (indices of) vertices that vertex (i,h) is adjacent to, and the distance of the vertex from the basepoint.
#v['adj'] lists adjacent vertices in clockwise order around v.
#vertices live in the square [0,1]x[0,1] (in window coordinates).
def BS_vert_dict(n,k):
	print('Calculating data for w_' + str(k) + ' in BS(1,' + str(n) + ')...')
	vert_dict = {}
	origin = (0, 0.5)
	x_scale = 1/(1 + 1/n**k)
	y_scale = 0.5

	#points on x-axis are special
	for i in range(n**k+2): 					#i is horizontal displacement times n**k
		v = g.window_coord((i/n**k, 0), origin, x_scale, y_scale)	#coordinates for vertex
		vert_dict[(i,0)] = {'coords':v, 'adj':[]}		#creates entry for (i,0), and puts in coordinates for (i,0)
		if i%n == 0: 						#these vertices are adjacent to the vertex above
			vert_dict[(i,0)]['adj'].append((i/n, 1))
		if i < n**k+1:	 					#these vertices are adjacent to the vertex to the right
			vert_dict[(i,0)]['adj'].append((i+1, 0))
		if i%n == 1:						#these vertices are adjacent to the vertex below
			vert_dict[(i,0)]['adj'].append(((i-1)/n, -1))
		if i > 0: 						#these vertices are adjacent to the vertex to the left
			vert_dict[(i,0)]['adj'].append((i-1, 0))


	for h in range(-k,0): 							#h is "height" times k, all these points are below x-axis
		for i in range(n**(k-abs(h))+1): 				#i is horizontal displacement times n**(k-abs(h))
			v = g.window_coord((i/n**(k-abs(h))+1/n**k, h/k), origin, x_scale, y_scale)	#coordinates for vertex with that index
			if h == -1:
				vert_above = (i*n+1, 0)
			else:
				vert_above = (i*n, h+1)
			vert_dict[(i,h)] = {'coords':v, 'adj':[vert_above]}	#creates entry for (i,h), and puts in coordinates for (i,h) and the vertex above in the adjacency set.
			if i < n**(k-abs(h)): 					#these vertices are adjacent to the vertex to the right
				vert_dict[(i,h)]['adj'].append((i+1, h))
			if i%n == 0 and h > -k: 				#these vertices are adjacent to the vertex below
				vert_dict[(i,h)]['adj'].append((i/n, h-1))
			if i > 0: 						#these vertices are adjacent to the vertex to the left
				vert_dict[(i,h)]['adj'].append((i-1, h))
			
	
	for h in range(1,k+1): 							#h is "height" times k, these are points above x-axis
		for i in range(n**(k-abs(h))+1): 				#i is horizontal displacement times n**(k-abs(h))
			v = g.window_coord((i/n**(k-abs(h)), h/k), origin, x_scale, y_scale)		#coordinates for vertex
			vert_dict[(i,h)] = {'coords':v, 'adj':[(i*n, h-1)]}	#creates entry for (i,h), and puts in coordinates for (i,h) and the vertex below in the adjacency set.
			if i > 0: 						#these vertices are adjacent to the vertex to the left
				vert_dict[(i,h)]['adj'].append((i-1, h))
			if i%n == 0 and h < k: 					#these vertices are adjacent to the vertex above
				vert_dict[(i,h)]['adj'].append((i/n, h+1))
			if i < n**(k-abs(h)): 					#these vertices are adjacent to the vertex to the right
				vert_dict[(i,h)]['adj'].append((i+1, h))

	
	#calculating distances
	print('calculating distances...')
	add_distances(vert_dict, (0,0))		#compute distances from basepoint and add them to vert_dict
	print('done')

	#print(vert_dict)
	return (vert_dict, (0,0))



#given a vert_dict and the index of the basepoint of the diagram, it adds (or updates, if it already exists) a 'distance' entry to vert_dict[v] for each vertex v, which is the distance in the vKd from basepoint to v
def add_distances(vert_dict, basepoint):
	vert_list = list(vert_dict)
	num_vertices = len(vert_list)
	vert_dict[basepoint]['distance'] = 0
	distance_list = [[basepoint]]				#distance_list[n] is the list of vertices that is a distance n away from the basepoint
	num_seen_vertices = 1					#keeps track of the number of vertices we've updated.
	
	d = 0
	while num_seen_vertices < num_vertices:			#while there may still be more vertices to reach
		distance_list.append([])
		for v in distance_list[d]:			#for each vertex of distance d
			for w in vert_dict[v]['adj']:		#look at all its adjacent vertices
				if 'distance' not in vert_dict[w] or vert_dict[w]['distance'] > d+1 or w not in distance_list[vert_dict[w]['distance']]:
					#if the adjacent vertex doesn't have a 'distance' value, or the distance value is wrong (either too big or the vertex isn't in our distance_list), then we haven't seen this vertex yet
					vert_dict[w]['distance'] = d + 1		#so its distance is d + 1
					distance_list[d+1].append(w)			#it goes in the distance list
					num_seen_vertices = num_seen_vertices + 1	#and we've seen it now
		d = d + 1





#returns a list of 2-cells
#each 2-cell is a list of vertices, such that the list is a path around the boundary of the 2-cell counter-clockwise (excluding the repeated end vertex, which is the same as the start vertex)
#assumes that vert_dict[v]['adj'] gives vertices in clockwise order around v (starting anywere)
def find_cells(vert_dict):
	edges_seen = []
	cells = []
	
	for u in vert_dict:
		for v in vert_dict[u]['adj']:						#for each directed edge
			if (u,v) not in edges_seen:					#if we haven't seen it
				path = [u]						#initialize the path
				w = u							#w and x will change, with w the current vertex and x the next one
				x = v
				total_angle = 0						#need to keep track of angle inside 2cell so we can distinguish exterior 2cell
			
				while x != u:						#while we haven't reached the starting vertex
					edges_seen.append((w, x))			#we travel (w, x) and label it as seen
					path.append(x)					#add the next vertex to the path
					edge_index = vert_dict[x]['adj'].index(w)	#find the index of the edge (x, w); the next vertex in the adjacency list for x will be the next one in the clockwise direction (thus traveling around the 2-cell in the counter-clockwise direction).
					num_edges = len(vert_dict[x]['adj'])		#needed to wrap around in next line if w is last on list
					y = vert_dict[x]['adj'][(edge_index + 1) % num_edges]		#make next vertex the next one in clockwise direction
					
					total_angle = total_angle + g.clockwise_angle(vert_dict[w]['coords'], vert_dict[x]['coords'], vert_dict[y]['coords'])	#add angle inside 2cell at x

					w = x						#reset current vertex
					x = y						#reset next vertex

				edges_seen.append((w,u))				#we also need to mark the last edge, closing the boundary, as having been seen
				
				total_angle = total_angle + g.clockwise_angle(vert_dict[w]['coords'], vert_dict[u]['coords'], vert_dict[v]['coords'])	#also need to add in the angle at u, which we never did at the beginning

				#print('total angle: ' + str(total_angle) + '; predicted number of vertices based on angle: ' + str(total_angle/math.pi + 2))

				if total_angle < math.pi * len(path):		#larger angle sum means going around outside (boundary of vKd), not inside (2cell)
					cells.append(path)			#the path is complete, so we add it to the set
				#else:
					#print('angle too big, must be boundary path: ' + str(path))
	
	return cells


	


#given a vert_dict and basepoint (index), returns a dictionary of 2-cells, whose keys are the 2-cells (as tuples, since lists can't be keys) and values are a list of tuples of contour partners in the 2-cell
#includes the necessary information to determine how to draw contour segments even if 2-cells are not simply-bounded, since the order of the tuple plus the orientation of the 2-cell determines the <n-path around the boundary of the 2-cell.
#assumes that 'distance' entry is included in vert_dict
def contour_dict(vert_dict):
	cells = find_cells(vert_dict)							#gets list of 2-cells
	contour_dict = {}
	for cell in cells:								#for each cell
		contour_dict[tuple(cell)] = []						#add a dictionary entry
		
		for i in range(len(cell)):						#for each vertex in that cell
			v = cell[i]
			v_dist = vert_dict[v]['distance']				
			length = len(cell)
			j = (i+1) % length						#index of next vertex around boundary

			w = cell[j]
			w_dist = vert_dict[w]['distance']
			if w_dist < v_dist:						#special case for first vertex: distance must go down to have a partner in the direction you're going
				j = (j+1) % length						#go to next vertex around boundary
				
				while j != i:							#while we haven't gotten back to the given vertex
					w = cell[j]
					w_dist = vert_dict[w]['distance']
					if w_dist == v_dist:					#if you find a vertex along boundary with same distance
						contour_dict[tuple(cell)].append((v,w))		#it's a contour partner with v
						break						#we're done with v, since it can only have one partner this direction
					elif w_dist > v_dist:					#if we come to a vertex bigger than v, no partner this direction
						break
					j = (j+1) % length					#if still here, then distance is still smaller than that of v, so go to next vertex around boundary
	
	return contour_dict
			
			
	


#returns a 2-tuple: first entry a dictionary of 2-cells, each of which has its own dictionary, second entry a dictionary of cell indices for the first dictionary
#keys for cell_dict are integers (so that each 2-cell has a unique key, which doesn't change if any are removed), plus a 'boundary' key
#keys for each 2-cell dictionary are 'circuit' (giving the indices of each vertex in the boundary circuit in order), 'partners' (giving, for each vertex, the index in circuit of its contour partner in the counter-clockwise direction, or -1 if it doesn't have one), 'cuts' (cuts[i] is a list of indices j where there is a cut between vertices i and j in circuit, in reverse numerical order mod length (thus going clockwise around the 2-cell) starting from i-1 and going "down" to i+1), 'cut segments' (a list of the cuts in the cell as line segments (2-tuples of points)), and 'segments' (a dictionary, where for a distance n away from basepoint, segments[n] is a list of line segments (as 2-tuples of points) giving the contour segments at distance n).
#the following were removed, since they relied on circuit including midpoints between adjacent partners: 'bad angles' (a list of indices in circuit at which the interior angle is too big, so there will need to be at least one cut).
#'boundary' 2-cell dictionary only has the 'circuit' key.
#keys for cell_index_dict are directed edges, in the form of a tuple of adjacent vertex indices
#value for (v1, v2) in cell_index_dict is index for the 2-cell with edge (v1, v2) on its counterclockwise boundary circuit
#assumes that distance entries already exist in vert_dict
def cell_dict(vert_dict):

	edges_seen = []
	cell_dict = {}
	cell_index_dict = {}
	cell_index = 0

	#First find the cells, and gather up some other information on the way	
	for u in vert_dict:
		for v in vert_dict[u]['adj']:						#for each directed edge
			if (u,v) not in cell_index_dict:				#if we haven't seen it
				path = [u]						#will be the list of vertices in the boundary circuit in counter-clockwise order.

				w = u							#w and x will change, with w the current vertex and x the next one
				x = v
				total_angle = 0						#need to keep track of angle inside 2cell so we can distinguish exterior 2cell
			
				while x != u:						#while we haven't reached the starting vertex
					cell_index_dict[(w, x)] = cell_index		#we travel (w, x) and put it in dictionary
					path.append(x)					#add the next vertex to the path
					edge_index = vert_dict[x]['adj'].index(w)	#find the index of the edge (x, w); the next vertex in the adjacency list for x will be the next one in the clockwise direction (thus traveling around the 2-cell in the counter-clockwise direction).
					num_edges = len(vert_dict[x]['adj'])		#needed to wrap around in next line if w is last on list
					y = vert_dict[x]['adj'][(edge_index + 1) % num_edges]		#make next vertex the next one in clockwise direction
					
					angle_at_x = g.clockwise_angle(vert_dict[w]['coords'], vert_dict[x]['coords'], vert_dict[y]['coords'])	#interior angle at x
					total_angle = total_angle + angle_at_x		#add angle inside 2cell at x

					w = x						#reset current vertex
					x = y						#reset next vertex

				cell_index_dict[(w, u)] = cell_index			#we also need to put the last edge, closing the boundary, in the dictionary
				
				angle_at_u = g.clockwise_angle(vert_dict[w]['coords'], vert_dict[u]['coords'], vert_dict[v]['coords'])	#interior angle at u
				total_angle = total_angle + angle_at_u			#also need to add in the angle at u, which we never did at the beginning

				#print('total angle: ' + str(total_angle) + '; predicted number of vertices based on angle: ' + str(total_angle/math.pi + 2))

				#larger angle sum means going around outside (boundary of vKd), not inside (2cell)
				if total_angle < math.pi * len(path):		#if 2-cell
					#the path is complete, so we add it to the dictionary
					cell_dict[cell_index] = {'circuit':path}
					cell_index = cell_index + 1
				else:
					cell_dict['boundary'] = {'circuit':path}	#this is the boundary, so add it to cell_dict as boundary
					#cell_index_dict is wrong; we have to go back and note, for every edge in path, that it is on the boundary
					length = len(path)
					for i in range(length):
						cell_index_dict[(path[i], path[(i+1)%length])] = 'boundary'
					#print('angle too big, must be boundary path: ' + str(path))

	#now that the circuits have been figured out, add in the rest of the stuff
	cell_index_set = set(list(cell_dict))
	#print('cell_dict calling fix_cells')
	fix_cells(cell_index_set, cell_dict, vert_dict)

	#print('cell index dict: ' + str(cell_index_dict))

	return (cell_dict, cell_index_dict)



#given a *set* of cell indices, a cell_dict, and a vert_dict, assumes that cells in cells_to_fix have correct circuits, but may have incorrect or missing partners, bad_angles, cuts, and segments, and edits cell_dict to fill these entries in correctly
def fix_cells(cells_to_fix, cell_dict, vert_dict):
	#print('entering fix_cells')
	for cell_index in cells_to_fix:
		if cell_index == 'boundary':
			continue
		cell = cell_dict[cell_index]['circuit'].copy()	#don't want to change circuit in cell_dict
		#print('working on cell ' + str(cell))
		length = len(cell)
		coords = [vert_dict[v]['coords'] for v in cell]
		#print('coords for circuit: ' + str(coords))
		#find bad_angles, the list of indices in cell such that interior angle is too big (recall that it needs to be in order around the cell)
		bad_angles = []
		for i in range(length):
			if g.clockwise_angle(coords[i-1], coords[i], coords[(i+1)%length]) >= max_angle:
				bad_angles.append(i)
		#print('found bad angles')

		#find contour partners
		partners = []			#will be the list of indices of contour partners, with -1 for vertices that don't have a partner in the cell
		adj_partners = -1		#special case: if partners share an edge, we will want to add a vertex in middle of that edge (and mark it as having too big an angle) so that we will cut at it.  Value of adj_partners will be updated to index of 1st vertex in counter-clockwise order in cell (so increment for index of partner).
		#print('length: ' + str(length))
		for i in range(length):						#for each vertex in that cell
			v = cell[i]
			v_dist = vert_dict[v]['distance']				
			j = (i+1) % length					#index of next vertex around boundary

			w = cell[j]
			w_dist = vert_dict[w]['distance']
			if w_dist >= v_dist:					#special case for first vertex: distance must go down to have a partner in the direction you're going
				partners.append(-1)
			else:
				j = (j+1) % length				#go to next vertex around boundary
				
				partner_found = False
				while j != i:					#while we haven't come back to where we started
					w = cell[j]
					w_dist = vert_dict[w]['distance']
					if w_dist == v_dist:			#if you find a vertex along boundary with same distance
						partners.append(j)		#it's a contour partner with v, with j as its index

						if (j+1) % length == i:		#if partners are adjacent
							adj_partners = j	#mark it to be dealt with after loop, so we aren't changing cell as we loop over it

						partner_found = True		#we're done with v, since it can only have one partner this direction
						break
					elif w_dist > v_dist:			#if we come to a vertex bigger than v, no partner this direction
						break
					j = (j+1) % length			#if still here, then distance is still smaller than that of v, so go to next vertex around boundary
				if not partner_found:
					partners.append(-1)			#if made it out of loop without a partner, it doesn't have one.

		#partners is now correct for what circuit is, so add it to cell_dict
		cell_dict[cell_index]['partners'] = partners.copy()

		#print('found contour partners')

		#after loop, check if we have adjacent partners.  If so, we will need to insert vertex into the cell, insert value of new vertex in coords and partners (-1), and increment values bigger than adj_partners in partners and bad_angles (to keep them pointing at the same vertex in the cell)
		if adj_partners >= 0:
			v_1 = cell[adj_partners]
			v_2 = cell[(adj_partners+1) % length]
			cell.insert(adj_partners + 1, 'midpoint')		#add vertex named 'midpoint' to cell
			length = length + 1
			#new vertex is at midpoint of the edge between the partners
			new_coords = ((vert_dict[v_1]['coords'][0] + vert_dict[v_2]['coords'][0])/2, (vert_dict[v_1]['coords'][1] + vert_dict[v_2]['coords'][1])/2)
			coords.insert(adj_partners + 1, new_coords)
			partners.insert(adj_partners + 1, -1)	#new vertex has no partners

			#add entry for new vertex in vert_dict so that we can check coordinates with it as we always do
			vert_dict['midpoint'] = {'coords': new_coords, 'adj':[], 'distance':vert_dict[v_1]['distance'] + 0.5}

			#now we need to fix all the indices we just messed up by moving all the vertices in cell past the new vertex over by one
			for i in range(length):
				if partners[i] >= adj_partners+1:
					partners[i] = partners[i] + 1

			added = False							#will be updated in loop if we add new vertex index to bad_angles (will not happen if bad_angles is empty or does not contain an index bigger than the new vertex's index
			for i in range(len(bad_angles)):
				if bad_angles[i] >= adj_partners+1:	#if the bad angle comes after the first adjacent partner (counterclockwise)
					bad_angles.insert(i, adj_partners+1)	#put new angle in bad angles (and keep it sorted in increasing order)
					added = True
					for j in range(i+1, len(bad_angles)):
						bad_angles[j] = bad_angles[j]+1		#relies on the fact that bad_angles is sorted
					break
			if not added:
				bad_angles.append(adj_partners+1)

		#print('adjusted for adjacent partners.  new cell: ' + str(cell))

		#Now to determine where cuts go.  Remember that bad_angles is in increasing order, counterclockwise around the 2-cell, so it might be more efficient to look in the counterclockwise direction for vertices to finish the cut, as that is more likely to hit a bad vertex that we haven't made a different cut at yet.
		bad_angles = [((i-1)%length, i, (i+1)%length) for i in bad_angles]	#translate bad_angles into 3-tuple form for this loop (doesn't edit cell_dict entry)
		cuts = [[] for i in range(length)]		#cuts[i] is the list of indices j such that there is a cut between vertices i and j
		bad_angle_index = 0					#keep track of the index in bad_angles we are at. increment at the end of loop
		while bad_angle_index < len(bad_angles):	#we'll be editing bad_angles as we go, so we need to keep checking the length
			(h, i, j) = bad_angles[bad_angle_index]
			vertex_before = coords[h]
			bad_vertex = coords[i]
			vertex_after = coords[j]

			angle = g.clockwise_angle(vertex_before, bad_vertex, vertex_after)
			ideal_angle = min(angle/2, max_angle)						#best to be around angle bisector, if possible
			vertex_to_test = g.point_at_angle(vertex_before, bad_vertex, ideal_angle)		#start with arbitrary point at ideal angle, to be cut off by some edge
			index_to_test = -1	#placeholder, to create an error if index_to_test is not replaced
			least_distance = 2	#this is not the actual distance to the above point, it just needs to be bigger than anything to be cut off


			#print('considering bad angle ' + str((h, i, j)))
			#print('considering cut to point ' + str(vertex_to_test))
			
			#This loop tests if index_to_test gives a valid cut, and updates it to a better one to test if not; once an index goes through the loop without being replaced, it gives a valid cut
			#print('looking for edges that cut it off')
			replaced = True
			outer_loop_counter = 0
			while replaced:
				replaced = False					#updated when index_to_test is replaced by a better one

				current_index = i		#start here and go around, testing each edge to see if it cuts off ray to vertex_to_test
				prior_index = h
				next_index = -1			#index of next vertex around (accounting for cuts). value is placeholder so variable is defined before loop starts and checks it

				loop_counter = 0
				while next_index != i:	#loop assumes that it starts with correct current_index and prior_index, but has not updated next_index yet
							#goes around until it would start repeating edges
					#update next_index
					if cuts[current_index] == []:					#if no cuts at current vertex, go to next one around 2-cell
						next_index = (current_index + 1)%length
					elif prior_index in cuts[current_index]:			#if there are cuts and we just came from one
						n = cuts[current_index].index(prior_index)
						if len(cuts[current_index]) >= n+2:			#and there's another one after
							next_index = cuts[current_index][n+1]		#traverse that one
						else:							#if there isn't another after
							next_index = (current_index + 1)%length		#go to next vertex around boundary
					else:								#if there are cuts but we didn't come from one
						next_index = cuts[current_index][0]			#we came from the vertex before around the boundary, so take the first cut

					#print('considering edge from ' + str(current_index) + ' to ' + str(next_index))

					#if this edge can't cut it off because it has vertex_to_test or bad_vertex as endpoint, skip this iteration
					if coords[next_index] in [vertex_to_test, bad_vertex, vertex_after] or coords[current_index] == vertex_to_test:
						#print('skipping this one, because it has ' + str(index_to_test) + ' or ' + str(i) + ' as an endpoint!')
						prior_index = current_index	#update prior and current index before the loop starts over
						current_index = next_index
						continue

					#if next_index and current_index are on opposite sides of ray from bad_vertex to vertex_to_test, then at least one of them is in the allowed region for a cut.  Choose one closest to ideal_angle, if both are allowed
					if (not g.on_right_side(coords[next_index], bad_vertex, vertex_to_test) and g.on_right_side(coords[current_index], bad_vertex, vertex_to_test)) or (not g.on_left_side(coords[next_index], bad_vertex, vertex_to_test) and g.on_left_side(coords[current_index], bad_vertex, vertex_to_test)):
						#print('edge intersects line...')
						intersect_point = g.intersection(bad_vertex, vertex_to_test, coords[current_index], coords[next_index])
						#need to make sure this point is actually in the correct direction to know the line segments intersect
						if g.dot(g.vector(bad_vertex, vertex_to_test), g.vector(bad_vertex, intersect_point)) > 0:
							#print('edge intersects ray...')
							new_distance = g.distance(bad_vertex, intersect_point)

							#if edge between current and next index cuts off line segment from bad_vertex to tested vertex, then tested vertex won't work; find a new one that's in the correct region
							if new_distance < least_distance:
								#print('segment to index ' + str(index_to_test) + ' is cut off by a edge from vertex ' + str(current_index) + ' to vertex ' + str(next_index))
								angle_to_next_index = g.clockwise_angle(vertex_before, bad_vertex, coords[next_index])
								if 0 < angle_to_next_index < max_angle:
									angle_to_current_index = g.clockwise_angle(vertex_before, bad_vertex, coords[current_index])
									if 0 < angle_to_current_index < max_angle:
										#if both could work, choose the one closest to ideal_angle
										if abs(angle_to_next_index - ideal_angle) <= abs(angle_to_current_index - ideal_angle):
											index_to_test = next_index
										else:
											index_to_test = current_index
									else:
										index_to_test = next_index
								elif 0 < g.clockwise_angle(vertex_before, bad_vertex, coords[current_index]) < max_angle:
									index_to_test = current_index
								else:
									print('something went wrong: there is an edge cutting off cut to index ' + str(index_to_test) + ', but neither endpoint works for a cut.  returning...')
									return False

								#print('Trying new index ' + str(index_to_test))
								replaced = True
								least_distance = new_distance

					prior_index = current_index	#update prior and current index before the loop starts over
					current_index = next_index

					loop_counter = loop_counter + 1
					if loop_counter > 2*length:
						print('never got back to bad_vertex. returning...')
						return False

					#print('next_index: ' + str(next_index) + '; i: ' + str(i))

				#once we have gone around the 2-cell and found the closest edge that cuts off the current vertex_to_test, update it
				vertex_to_test = coords[index_to_test]		#will give error if index_to_test was never updated from -1
				least_distance = g.distance(bad_vertex, vertex_to_test)	#need to update least_distance, because it may be different from distance to point of intersection of edge (which is what it was set at)

				#print('updated vertex_to_test to vertex ' + str(index_to_test))

				outer_loop_counter = outer_loop_counter + 1
				if outer_loop_counter > 5:
					print('way too many replacements happening, something has gone wrong.  returning...')
					return False

			#print('settled on cut from index ' + str(i) + ' to index ' + str(index_to_test))
			
			if index_to_test in cuts[i]:
				print('error: tried to add a cut that already exists.  something went wrong...')
				return False


			#since we've left the loop, the current index_to_test and vertex_to_test give a cut
			#add the cut to cuts[i] and cuts[index_to_test]
			#need to put it in appropriate spot in list to keep it sorted
			added = False
			for k in range(len(cuts[i])):
				if (i-cuts[i][k])%length > (i-index_to_test)%length:
					cuts[i].insert(k, index_to_test)
					added = True
					break
			if not added:
				cuts[i].append(index_to_test)

			added = False
			for k in range(len(cuts[index_to_test])):
				if (index_to_test-cuts[index_to_test][k])%length > (index_to_test-i)%length:
					cuts[index_to_test].insert(k, i)
					added = True
					break
			if not added:
				cuts[index_to_test].append(i)

			#if new cut splits a (future) bad angle, remove it from bad_angles and test if you need to add one of the two new angles.  
			k = bad_angle_index + 1
			while k < len(bad_angles):		#for rest of bad angles around the diagram, we might have split one of them...
				angle = bad_angles[k]
				if index_to_test == angle[1]:				#if we hit a vertex with a bad angle...
					#and we split the angle
					if g.clockwise_angle(coords[angle[0]], coords[index_to_test], coords[i]) < g.clockwise_angle(coords[angle[0]], coords[index_to_test], coords[angle[2]]):
						#print('new cut also splits the bad angle ' + str(angle) + '!  removing...')
						bad_angles.remove(angle)			#then angle is replaced by two new angles, at most one of which is bad
						if g.clockwise_angle(coords[angle[0]], vertex_to_test, bad_vertex) >= max_angle:	#if first new angle is bad
							bad_angles.insert(k, (angle[0], index_to_test, i))	#add it to bad angles in place of old one
							#print('added new bad angle ' + str((angle[0], index_to_test, i)))
						elif g.clockwise_angle(bad_vertex, vertex_to_test, coords[angle[2]]) >= max_angle:	#if second new angle is bad
							bad_angles.insert(k, (i, index_to_test, angle[2]))	#add it to bad angles in place of old one
							#print('added new bad angle ' + str((i, index_to_test, angle[2])))
						break	#since we've split this angle, we won't split any others
				k = k+1


			#angle vertex_before-bad_vertex-vertex_to_test is good now, but what about vertex_to_test-bad_vertex-vertex_after?  Need to test this.  If that one is still bad, need to do the whole process again with vertex_before replaced with vertex_to_test.  Need to put this whole process inside a loop, as arbitrarily many cuts may be necessary at a given vertex.
			if g.clockwise_angle(vertex_to_test, bad_vertex, vertex_after) >= max_angle:		#if new angle is bad
				bad_angles.insert(bad_angle_index + 1, (index_to_test, i, j))			#insert it to be the next one checked
				#print('angle ' + str((index_to_test, i, j)) + ' is still bad.  Checking it next...')

			bad_angle_index = bad_angle_index + 1		#increment counter


		#print('found cuts')

		#since cell_dict['circuit'] does not contain potential midpoints, store cuts as just a list of line segments
		cut_list = []
		cuts_seen = []		#keep track of which cuts we've seen, to avoid repeats that just have the order of the endpoints swapped
		for i, index_list in enumerate(cuts):
			for j in index_list:
				if (j, i) not in cuts_seen:	#if haven't seen same segment with endpoints swapped
					cut_list.append((coords[i], coords[j]))
				cuts_seen.append((i, j))
		
		cell_dict[cell_index]['cut segments'] = cut_list

		#now find segments based on everything we've done so far
		distances = [vert_dict[v]['distance'] for v in cell]
		cell_dict[cell_index]['segments'] = segments(cell, partners, bad_angles, cuts, coords, distances)

		#if a midpoint was added between two adjacent partners, edit cuts so that the indices match up with those in the circuit in cell_dict
		if adj_partners >= 0:
			cuts.pop(adj_partners+1)			#this is the index of the midpoint, so remove its list from cuts
			#make a new list for each entry of cuts, that removes any adj_partners+1 and converts indices back into what they should be for the circuit in cell_dict
			for i, index_list in enumerate(cuts):
				new_index_list = []
				for index in index_list:
					if index < adj_partners+1:
						new_index_list.append(index)
					elif index > adj_partners+1:
						new_index_list.append(index-1)
				cuts[i] = new_index_list

		#add cuts to cell_dict
		cell_dict[cell_index]['cuts'] = cuts

		#remember to delete the potential extraneous midpoint from vert_dict
		if adj_partners >= 0:
			del vert_dict['midpoint']
		#print('done with cell ' + str(cell_index))
	#print('new cell dict: ' + str(cell_dict))



#returns a dictionary where key n corresponds to a list of segments that make up the contour lines of distance n away from the basepoint.
#a segment is a tuple (p, q), where p and q are points that we are to draw a segment between
def segments(circuit, partners, bad_angles, cuts, coords, distances):
	segments = {}

	length = len(circuit)

	circum = 0			#circumference of 2-cell
	edge_lengths = []		#edge_lengths[i] is the length of the edge from i to (i+1)%length
	for i in range(length):
		edge_lengths.append(g.distance(coords[i], coords[(i+1)%length]))
		circum = circum + edge_lengths[i]

	for i in range(length):
		p = partners[i]		#p is index of i's counterclockwise partner, if it has one, -1 otherwise
		if p >= 0:			#if it has a partner, we'll put some segments in segments[distance]

			height = distances[i]
			if height not in segments:
				segments[height] = []

			#if these are adjacent partners (with midpoint between them), do something a little different to make the segments look nice
			#still consistent so that these segments won't cross others, given that we will only potentially bring them closer to the midpoint
			#makes 2 segments, with point joining them the midpoint pushed out by at most the minimum from any of the cuts at the midpoint or an appropriate fraction of the length of the edge from p to i
			if (p+2)%length == i and distances[(p+1)%length] == height + 0.5:
				
				p_to_i = g.distance(coords[p], coords[i])
				proportion = p_to_i/circum

				v = g.vector(coords[p], coords[i])
				push_vector = g.normal((v[1], -v[0]))	#normal vector that is vector from p to i rotated counterclockwise by pi/2

				push_distance = p_to_i/4	#first guess at what would look reasonable, if cuts allow
				for c in cuts[(p+1)%length]:
					cut_vector = g.vector(coords[(p+1)%length], coords[c])
					push_distance = min(push_distance, g.dot(push_vector, cut_vector)*proportion)
				
				push_point = g.add(g.mult(push_distance, push_vector), coords[(p+1)%length])
				segments[height].append((coords[p], push_point))
				segments[height].append((push_point, coords[i]))

				continue

			current_vertex = coords[i]	#start at vertex i

			for j in [(i+k)%length for k in range(1, (p-i)%length)]:			#for indices j from i to p
				for c in cuts[j]:									#for each index c with a cut between j and c
					if c in [(p+k)%length for k in range(1, (i-p)%length)]:	#if c is in circuit from p to i (and thus the cut separates i and p)

						#the following is a method of calculating proportion based on how many contour segments might traverse the cut.  It gave pretty ragged edges, so I thought I might try to do better.
						#max_partners_under = min(abs(j-i), abs(j-p)) - 1
						#if min(abs(c-i), abs(c-p)) >= 1 and max(abs(c-i), abs(c-p)) >= 2:	#in this case, can be a vert above with two partners on other side
						#	max_partners_over = min(abs(c-i), abs(c-p))
						#else:
						#	max_partners_over = min(abs(c-i), abs(c-p)) - 1
						#proportion = (max_partners_under + 1)/(max_partners_under + max_partners_over + 2)
							
						i_to_p = sum([edge_lengths[(i+k)%length] for k in range((p-i)%length)]) #distance from i to p over boundary
						proportion = i_to_p/circum

						vert_j = coords[j]
						vert_c = coords[c]
							
						cut_vector = (vert_c[0] - vert_j[0], vert_c[1] - vert_j[1])
						next_vertex = (vert_j[0] + proportion*cut_vector[0], vert_j[1] + proportion*cut_vector[1])

						segments[height].append((current_vertex, next_vertex))

						current_vertex = next_vertex

			#having exited the loop, there are no more cuts to go through, and current_vertex is where the contour segment intersects the last cut.  So draw a segment from current_vertex to p
			segments[height].append((current_vertex, coords[p]))

	return segments



#edits vert_dict, cell_dict, and cell_index_dict to correctly represent a new diagram with an added 2-cell with boundary path \cdot \bar{path}, obtained by cutting the diagram along path and gluing in a 2-cell
#path is assumed to be a list of vertex indices such that adjacent entries are adjacent in vert_dict and there are no repeated vertices
#path cannot be a single edge (since it's assumed that there are no double edges)
def split_path(path, vert_dict, cell_dict, cell_index_dict):
	print('splitting path...')
	
	coords = [vert_dict[v_index]['coords'] for v_index in path]

	verts_to_fix = set()	#a set of vertices whose distances will have to be fixed later
	#not using cells_to_fix for now, because it's too hard to figure out all the cells that need to change
	#cells_to_fix = set()	#a set of cell indices whose contour partners, bad_angles, and cuts will have to be re-computed later (made it a set because it may be that the same cell gets added multiple times, if it has multiple edges in path

	#create the new cell in cell_dict with boundary path \cdot \bar{path}
	right_path = path[0:1] + [(v, 'right') for v in path[1:len(path)-1]] + path[len(path)-1:len(path)]
	reversed_left_path = [(v, 'left') for v in path[1:len(path)-1]]
	reversed_left_path.reverse()
	new_circuit = right_path + reversed_left_path

	cell_index_list = list(cell_dict)
	cell_index_list.remove('boundary')			#need to remove 'boundary' key so that remaining keys are all numbers to find a new index
	new_cell_index = cell_index_list.pop() + 1		#last entry has the biggest index, so add 1 to it and we get a new unique index
	cell_dict[new_cell_index] = {'circuit':new_circuit}
	#cells_to_fix.add(new_cell_index)

	basepoint_split = False		#updated if basepoint is in the middle of the path, and would be copied.  will need to choose which copy is new basepoint

	for i in range(1, len(path)-1):	#won't be changing endpoints
		if i == 1:	#in this case, u is start vertex, so does not have left and right copies
			u_left_index = path[i-1]
			u_right_index = path[i-1]
		else:
			u_left_index = (path[i-1], 'left')
			u_right_index = (path[i-1], 'right')

		v_index = path[i]	#current vertex
		w_index = path[i+1]	#vertex after
		u_coords = coords[i-1]
		v_coords = coords[i]
		w_coords = coords[i+1]
		v_distance = vert_dict[v_index]['distance']	#for now, we'll keep distance the same.  will be updated later
		v_adj = vert_dict[v_index]['adj']

		if v_distance == 0:	#if v is the basepoint of the diagram
			basepoint_split = True
			old_basepoint = v_index

		print('working on splitting vertex ' + str(v_index))

		left_angle = g.clockwise_angle( u_coords, v_coords,  w_coords )
		right_angle = 2*math.pi - left_angle

		#print('left angle: ' + str(left_angle))

		#get vectors that give direction we should push copies of v, to the left and right of path (by bisecting the angle at v)
		#print('point in correct direction: ' + str(point_at_angle( u_coords, v_coords, left_angle/2 )))
		#print('as a vector: ' + str(vector( v_coords, point_at_angle( u_coords, v_coords, left_angle/2 ) )))
		left_push_vector = g.normal( g.vector( v_coords, g.point_at_angle( u_coords, v_coords, left_angle/2 ) ) )
		right_push_vector = g.mult(-1, left_push_vector)
		print('left push vector: ' + str(left_push_vector))
		print('right push vector: ' + str(right_push_vector))

		#guess a good distance to push off to left.  if there are edges on the left side of path, use them to gauge distance.  if not, but there are cuts, use them to gauge distance.  if neither, use u and w.  this is just a heuristic; it is not guaranteed to give push-offs that result in non-intersecting edges.
		if left_angle < math.pi:								#if angle is small enough
			print('angle is less than max angle of ' + str(max_angle) + ' radians')
			left_push_dist = min(g.dot(left_push_vector, g.vector(v_coords, u_coords)), g.dot(left_push_vector, g.vector(v_coords, w_coords)))/2		#start with pushing it halfway out of the corner between u and w
			
		else:
			#in this case, there will be a cut or edge on the left side of v, or it's on the boundary, so just make left_push_dist something way too big to be replaced later
			left_push_dist = 10

		print('initial guess at left push distance: ' + str(left_push_dist))
		
		#now see if there are edges on left side to use to gauge distance
		edges_considered = False			#updated to True in loop if there are edges between u and w
		#will also find all vertices adjacent to v on left side of path (including u and w), so initialize here
		left_v_adj = [u_left_index]
		j = (v_adj.index(u_left_index)+1) % len(v_adj)	#start at first vertex adjacent to v clockwise from u
		end_index = v_adj.index(w_index)		#finish at w
		while j != end_index:				#for edges between u and w
			left_v_adj.append(v_adj[j])		#add vertex to left_v_adj
			#update push distance to 1/3 the way up any edge (since 1/3 on each side gives 2/3 between each vertex)
			left_push_dist = min( left_push_dist, g.dot( left_push_vector, g.vector( v_coords, vert_dict[v_adj[j]]['coords']) ) / 3 )
			edges_considered = True
			j = (j+1) % len(v_adj)			#go to next index around

		left_v_adj.append(w_index) #finish off left_v_adj with w
		#don't do this for w, because the vertex (w_index, 'left') doesn't exist yet
		#if i+1 == len(path)-1:	#if w is last vertex in path, it doesn't have a left copy
		#	left_v_adj.append(w_index)
		#else:		#otherwise, v-left will be adjacent to w-left
		#	left_v_adj.append((w_index, 'left'))


		#if there were no edges between u and w, try for cuts, assuming there is actually a 2-cell to the left, rather than the boundary of the diagram
		if not edges_considered and cell_index_dict[(v_index, w_index)] != 'boundary':
			cell = cell_dict[cell_index_dict[(v_index, w_index)]]
			print('cell for edge ' + str((v_index, w_index)) + ': ' + str(cell))
			circuit = cell['circuit']
			for j in range(len(circuit)):
				if circuit[j] == v_index and circuit[(j+1)%len(circuit)] == w_index:	#v can occur multiple times in the circuit, so we need to make sure this is the index we want
					i_v = j
					break

			cuts_at_v = cell['cuts'][i_v]
			if cuts_at_v != []:
				circum = sum( [g.distance( vert_dict[circuit[j]]['coords'], vert_dict[circuit[(j+1)%len(circuit)]]['coords'] ) for j in range(len(circuit))] )
				proportion = (g.distance(u_coords, v_coords) + g.distance(v_coords, w_coords))/circum	#same proportion used for contour segments
				for j in cuts_at_v:
					#potentially update left_push_dist to push vertex up to the point on the contour segment at that cut
					left_push_dist = min( left_push_dist, proportion*g.dot(left_push_vector, g.vector(v_coords, vert_dict[circuit[j]]['coords'])) )


		left_v_coords = g.add(v_coords, g.mult(left_push_dist, left_push_vector))	#push v off to the left by left_push_dist units
		if left_v_coords[0] < 0 or left_v_coords[0] > 1 or left_v_coords[1] < 0 or left_v_coords[1] > 1:
			left_push_dist = min(v_coords[0], v_coords[1], 1 - v_coords[0], 1 - v_coords[1])	#just make sure it doesn't leave [0,1]x[0,1]
			left_v_coords = g.add(v_coords, g.mult(left_push_dist, left_push_vector))

		print('coords for ' + str((v_index, 'left')) + ': ' + str(left_v_coords))

		vert_dict[(v_index, 'left')] = {'coords':left_v_coords, 'adj':left_v_adj, 'distance':v_distance}	#add new vertex to vert_dict
		verts_to_fix.add((v_index, 'left'))	#distance will need to be fixed later

		#for each vertex x adjacent to left-v, correct its adjacency list, and replace the entries (v, x) and (x, v) in cell_index_dict
		for x_index in left_v_adj:
			x_adj = vert_dict[x_index]['adj']	#going to correct x_adj
			if v_index in x_adj:			#need to check, because some might already be updated, and don't want to think about which ones
				i_v = x_adj.index(v_index)
				x_adj.pop(i_v)
				if i == 1 and x_index == path[0]:	#start of path means we need to add both copies, and in the correct, clockwise order
					x_adj.insert(i_v, (v_index, 'left'))
					x_adj.insert(i_v+1, (v_index, 'right'))
				
				elif x_index == w_index:		#same with w
					x_adj.insert(i_v, (v_index, 'right'))
					x_adj.insert(i_v+1, (v_index, 'left'))

				else:					#otherwise, just add the left one
					x_adj.insert(i_v, (v_index, 'left'))

			#now fix (v,x) and (x,v) entries in cell_index_dict
			#if x = u, then (v, x) corresponds to old cell on right, and should be changed to new cell.
			if x_index == u_left_index:
				cell_index_dict[( (v_index, 'left'), x_index )] = new_cell_index
			else:
				cell_index_dict[( (v_index, 'left'), x_index )] = cell_index_dict[(v_index, x_index)]
			#del cell_index_dict[(v_index, x_index)]	#should not do this yet, since right will need it still

			#if x is w, then (x, v-left) should correspond to the new cell, rather than the old cell on the right
			if x_index == w_index:
				cell_index_dict[( x_index, (v_index, 'left') )] = new_cell_index
			else:
				cell_index_dict[( x_index, (v_index, 'left') )] = cell_index_dict[(x_index, v_index)]
			#del cell_index_dict[(x_index, v_index)]	#should not do this yet, since right will need it still
			

		#now to edit all the cells that once had v in them, and should now have left-v
		for x_index in left_v_adj[1:]:					#for all vertices adjacent to left-v, except for u,
			cell_index = cell_index_dict[((v_index, 'left'), x_index)]	#the directed edge from v to x is on a 2-cell that should contain left-v
			circuit = cell_dict[cell_index]['circuit']
			#if cell_index != 'boundary':				#if there is a 2-cell here, not the boundary of the diagram,
				#cells_to_fix.add(cell_index)			#will need to fix this cell, since it contained v
			
			v_index_in_circuit = circuit.index(v_index)
			#since v might appear multiple times in circuit, need to check this is correct index
			while v_index_in_circuit < len(circuit) - 1 and circuit[v_index_in_circuit + 1] != x_index:
				v_index_in_circuit = circuit[v_index_in_circuit + 1:].index(v_index)

			#now that we have its index, replace v with left-v in circuit
			circuit.pop(v_index_in_circuit)
			circuit.insert(v_index_in_circuit, (v_index, 'left'))


		#now we need to do all this again for the right copy of v:

		#guess a good distance to push off to right.  if there are edges on the right side of path, use them to gauge distance.  if not, but there are cuts, use them to gauge distance.  if neither, use u and w.  this is just a heuristic; it is not guaranteed to give push-offs that result in non-intersecting edges.
		if right_angle < math.pi:								#if angle is small enough
			right_push_dist = min(g.dot( right_push_vector, g.vector(v_coords, w_coords) ), g.dot( right_push_vector, g.vector(v_coords, w_coords) ))/2	#start with pushing it halfway out of corner between u and w	
		else:
			#in this case, either v is on the boundary, or there will be a cut or edge on the right side of v, so just make right_push_dist something way too big to be replaced later
			right_push_dist = 10
		
		#now see if there are edges on right side to use to gauge distance
		edges_considered = False			#updated to True in loop if there are edges between u and w
		right_v_adj = [w_index] #will also find all vertices adjacent to v on right side of path (including u and w), so initialize here
		#don't do this for w, since (w_index, 'right') doesn't exist yet
		#if i+1 == len(path)-1:	#if w is last vertex in path, it doesn't have a right copy
		#	right_v_adj = [w_index]
		#else:		#otherwise, v-right will be adjacent to w-right
		#	right_v_adj = [(w_index, 'right')]

		j = (v_adj.index(w_index)+1) % len(v_adj)	#start at first vertex adjacent to v clockwise from w
		end_index = v_adj.index(u_right_index)		#finish at u
		while j != end_index:				#for edges between w and u
			right_v_adj.append(v_adj[j])		#add vertex to right_v_adj
			#update push distance to 1/3 the way up any edge (since 1/3 on each side gives 2/3 between each vertex)
			right_push_dist = min( right_push_dist, g.dot( right_push_vector, g.vector( v_coords, vert_dict[v_adj[j]]['coords']) ) / 3 )
			edges_considered = True
			j = (j+1) % len(v_adj)			#go to next index around

		#finish off right_v_adj with u
		right_v_adj.append(u_right_index)
		print('verts adjacent to ' + str((v_index, 'right')) + ': ' + str(right_v_adj))

		#if there were no edges between w and u, try for cuts, assuming there is a 2-cell to the right, not the boundary of the diagram
		if not edges_considered and cell_index_dict[(v_index, u_right_index)] != 'boundary':			
			cell = cell_dict[cell_index_dict[(v_index, u_right_index)]]
			circuit = cell['circuit']
			for j in range(len(circuit)):
				if circuit[j] == v_index and circuit[(j+1)%len(circuit)] == u_right_index:	#v can occur multiple times in the circuit, so we need to make sure this is the index we want
					i_v = j
					break

			cuts_at_v = cell['cuts'][i_v]
			if cuts_at_v != []:
				circum = sum( [g.distance( vert_dict[circuit[j]]['coords'], vert_dict[circuit[(j+1)%len(circuit)]]['coords'] ) for j in range(len(circuit))] )
				proportion = (g.distance(w_coords, v_coords) + g.distance(v_coords, u_coords))/circum	#same proportion used for contour segments
				for j in cuts_at_v:
					#potentially update right_push_dist to push vertex up to the point on the contour segment at that cut
					right_push_dist = min( right_push_dist, proportion*g.dot(right_push_vector, g.vector(v_coords, vert_dict[circuit[j]]['coords'])) )

		right_v_coords = g.add(v_coords, g.mult(right_push_dist, right_push_vector))	#push v off to the right by right_push_dist units
		if right_v_coords[0] < 0 or right_v_coords[0] > 1 or right_v_coords[1] < 0 or right_v_coords[1] > 1:
			right_push_dist = min(v_coords[0], v_coords[1], 1 - v_coords[0], 1 - v_coords[1])	#just make sure it doesn't leave [0,1]x[0,1]
			right_v_coords = g.add(v_coords, g.mult(right_push_dist, left_push_vector))

		print('coords for ' + str((v_index, 'right')) + ': ' + str(right_v_coords))

		vert_dict[(v_index, 'right')] = {'coords':right_v_coords, 'adj':right_v_adj, 'distance':v_distance}	#add new vertex to vert_dict
		verts_to_fix.add((v_index, 'right'))	#distance will need to be fixed later

		del vert_dict[v_index]	#having added both copies of v to vert_dict, we can now remove the original


		#for each vertex x adjacent to right-v, correct its adjacency list, and replace the entries (v, x) and (x, v) in cell_index_dict
		for x_index in right_v_adj:
			x_adj = vert_dict[x_index]['adj']	#going to correct x_adj
			if v_index in x_adj:			#need to check, because some might already be updated, and don't want to think about which ones
				i_v = x_adj.index(v_index)
				x_adj.pop(i_v)
				#already added v-right if x is endpoint of path above, so don't need to worry about it here
				x_adj.insert(i_v, (v_index, 'right'))

			#now fix (v,x) and (x,v) entries in cell_index_dict
			#if x = w, then (v, x) corresponds to old cell on left, and should be changed to new cell.
			if x_index == w_index:
				cell_index_dict[( (v_index, 'right'), x_index )] = new_cell_index
			else:
				cell_index_dict[( (v_index, 'right'), x_index )] = cell_index_dict[(v_index, x_index)]
			del cell_index_dict[(v_index, x_index)]

			#if x = u, then (x, v-right) should correspond to the new cell, rather than the old cell on the left
			if x_index == u_right_index:
				cell_index_dict[( x_index, (v_index, 'right') )] = new_cell_index
			else:
				cell_index_dict[( x_index, (v_index, 'right') )] = cell_index_dict[(x_index, v_index)]
			del cell_index_dict[(x_index, v_index)]


		#now to edit all the cells that once had v in them, and should now have right-v
		for x_index in right_v_adj[1:]:					#for all vertices adjacent to right-v, except for w,
			print('updating circuit containing edge ' + str(((v_index, 'right'), x_index)))
			cell_index = cell_index_dict[((v_index, 'right'), x_index)]	#the directed edge from v to x is on a 2-cell that should contain right-v
			circuit = cell_dict[cell_index]['circuit']
			print('circuit: ' + str(circuit))

			#if cell_index != 'boundary':				#if there is a 2-cell here, not boundary of the diagram,
			#	cells_to_fix.add(cell_index)			#will need to fix this cell, since it contained v
			
			v_index_in_circuit = circuit.index(v_index)
			#since v might appear multiple times in circuit, need to check this is correct index
			while v_index_in_circuit < len(circuit) - 1 and circuit[v_index_in_circuit + 1] != x_index:
				v_index_in_circuit = circuit[v_index_in_circuit + 1:].index(v_index)

			#now that we have its index, replace v with left-v in circuit
			circuit.pop(v_index_in_circuit)
			circuit.insert(v_index_in_circuit, (v_index, 'right'))

		print('split vertex ' + str(v_index))
		#print('current vert_dict: ' + str(vert_dict))
		#print('current cell_dict: ' + str(cell_dict))

		
	#now just need to fix distances of some vertices and partners, bad angles, and cuts of some cells

	#first, need to choose basepoint if it was split
	if basepoint_split:
		if (old_basepoint, 'left') in cell_dict['boundary']['circuit']:		#the one on the boundary is the new basepoint
			add_distances(vert_dict, (old_basepoint, 'left'))		#just update all distances, since fix_distances algorithm won't work
		elif (old_basepoint, 'right') in cell_dict['boundary']['circuit']:
			add_distances(vert_dict, (old_basepoint, 'right'))
		else:
			print('neither basepoint option is on boundary; something went wrong. returning...')
			return False
	else:
		fix_distances(verts_to_fix, vert_dict)		#otherwise, we can use fix_distances

	print('split_path calling fix_cells')
	#if want to save time using cells_to_fix, would have to add in all the cells who have vertices whose distances changed with fix_distances as well.  Doesn't seem worth it for now.
	#fix_cells(cells_to_fix, cell_dict, vert_dict)
	cell_index_set = set(list(cell_dict))
	fix_cells(cell_index_set, cell_dict, vert_dict)

	#print('new cell index dict: ' + str(cell_index_dict))
					


#given a *set* of vertex indices and a vert_dict, assumes that vertices in verts_to_fix are the points in the interior of a subdiagram that replaced a previous subdiagram of the vKd, and that distances of other vertices in vert_dict were correct in the old diagram.  Vertices in verts_to_fix may or may not have distance entries in vert_dict.
#edits vert_dict to make the distances of verts_to_fix correct, as well as the distances to any vertices affected by distance to these vertices, resulting in correct distances for all vertices in the diagram
#assumes that verts_to_fix does not contain the basepoint
#edits verts_to_fix in the process
def fix_distances(verts_to_fix, vert_dict):
	new_boundary = set()				#first find all the vertices adjacent to verts_to_fix
	for v in verts_to_fix:
		for w in vert_dict[v]['adj']:
			if w not in verts_to_fix:
				new_boundary.add(w)

	#now put things from boundary into verts_to_fix if they didn't have a geodesic outside verts_to_fix, and add their adjacent vertices to new_boundary
	#iterate until no changes have been made, indicating everything on the boundary had a geodesic outside verts_to_fix
	boundary = set()
	while boundary != new_boundary:

		boundary = new_boundary	#update boundary
		new_boundary = set()	#empty new_boundary so we can fill it with new stuff

		for v in boundary:
			#determine if it is either the basepoint or has a neighbor outside of verts_to_fix with smaller distance
			if vert_dict[v]['distance'] == 0:
				keep_in_boundary = True
			else:
				keep_in_boundary = False
				for w in vert_dict[v]['adj']:
					if w not in verts_to_fix and vert_dict[w]['distance'] < vert_dict[v]['distance']:
						keep_in_boundary = True
						break

			#if so, it stays in boundary (so put it in new_boundary)
			if keep_in_boundary:
				new_boundary.add(v)
			#if not, it needs to be fixed as well, and its neighbors go in new_boundary
			else:
				verts_to_fix.add(v)
				for w in vert_dict[v]['adj']:
					if w not in verts_to_fix:
						new_boundary.add(w)

	#boundary now has all vertices adjacent to verts_to_fix, and they all have an old geodesic outside verts to fix.  They may not have the correct distance, but this can only happen if their new distance is smaller, and this can be dealt with more easily.  however, the smallest distance vertex on boundary must be correct.
	#sort boundary into sets of the same distance, and put in a dictionary with distance the keys
	boundary_dict = {}
	for v in boundary:
		v_dist = vert_dict[v]['distance']
		if v_dist in boundary_dict:
			boundary_dict[v_dist].add(v)
		else:
			boundary_dict[v_dist] = {v}

	#now we start from the smallest distance vertex on boundary and move in, updating distances.
	#iterate until there is nothing left in boundary (and therefore nothing left in verts_to_fix)
	while boundary_dict != {}:
		min_dist = min(list(boundary_dict))
		for v in boundary_dict[min_dist]:			#for minimum distance boundary vertices
			for w in vert_dict[v]['adj']:			#look at its neighbors
				if w in verts_to_fix:			#if one needs to be fixed,
					#set its distance and remove it from verts_to_fix
					vert_dict[w]['distance'] = min_dist + 1
					verts_to_fix.remove(w)
					#add it to boundary_dict for next loop
					if min_dist+1 in boundary_dict:
						boundary_dict[min_dist+1].add(w)
					else:
						boundary_dict[min_dist+1] = {w}
				elif vert_dict[w]['distance'] > min_dist + 1:		#if outside verts_to_fix but its distance is too big
					vert_dict[w]['distance'] = min_dist + 1		#set its distance
					#add it to boundary_dict for next loop
					if min_dist+1 in boundary_dict:
						boundary_dict[min_dist+1].add(w)
					else:
						boundary_dict[min_dist+1] = {w}

		#done with all min_dist vertices, so remove from boundary_dict
		del boundary_dict[min_dist]
