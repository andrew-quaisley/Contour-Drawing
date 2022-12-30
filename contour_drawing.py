from tkinter import *
import copy
import geometry as g
import vKd_calculations as calc
import drawing as d

#----- Main Program -----

#Do calculations first

#Set what vKd you want to draw here
vert_dict, basepoint = calc.BS_vert_dict(2,6)

#simple vKd for testing:
#vert_dict = {0:{'coords':(0,0), 'adj':[6, 1]}, 1:{'coords':(0, 0.5), 'adj':[0, 7, 2]}, 2:{'coords':(0,1), 'adj':[1, 3]}, 3:{'coords':(0.5, 1), 'adj':[2, 7, 4]}, 4:{'coords':(1,1), 'adj':[3, 5]}, 5:{'coords':(1, 0.5), 'adj':[4, 7, 6]}, 6:{'coords':(1, 0), 'adj':[5, 0]}, 7:{'coords':(0.5, 0.5), 'adj':[5, 3, 1]}}
#add_distances(vert_dict, 0)

cell_dict, cell_index_dict = calc.cell_dict(vert_dict)

#Then set up window and draw things
window = Tk()

#getting screen width and height of display
win_width = window.winfo_screenwidth() 
win_height = window.winfo_screenheight()

#set scale factor and location of origin (top left corner) on screen
#min_dimension = min(win_width, win_height)
x_scale = 0.95*win_width
y_scale = 0.85*win_height
origin = (0.025*win_width, 0.025*win_height)

#setting tkinter window size
resolution = str(win_width) + 'x' + str(win_height) + '+0+0'
window.geometry(resolution)

window.title('Contour Line Drawing Tool')
frame = Frame(window)
frame.pack(expand=True, fill=BOTH)

canvas = Canvas(frame, bg='white', width=win_width, height=win_height)

canvas.pack(expand=True, fill=BOTH)
window.update()


#draw stuff here

d.draw_vKd(vert_dict, canvas, origin, x_scale, y_scale)

d.draw_contour_lines_better(cell_dict, canvas, origin, x_scale, y_scale)

mode = 'view'			#mode options: 'view', 'select path', 'edit'
highlighted_object = False	#if some canvas object is highlighted (but not selected, i.e., the mouse is hovering over it), gets set to the tuple (canvas id, vertex/endpoint index, canvas coordinates of vertex/endpoint, [canvas coordinates of 2nd endpoint]) for that object
selected_objects = []		#list of ids of canvas objects that have been selected, in the order they were selected
path_to_split = []		#empty unless in path splitting mode, in which case it lists the vertex indices of the path that is currently selected
drag_coords = False		#the current canvas coordinates of a vertex being dragged.  False if no vertex is currently being dragged

distances_showing = True
graph_showing = True
contour_lines_showing = True

history = [(copy.deepcopy(vert_dict), copy.deepcopy(cell_dict), copy.deepcopy(cell_index_dict))]	#stores all info for undo/redo function
history_index = 0	#index of the current state of diagram in history


#-------- Bound functions for user interaction ----------


#### 'select path' mode functions ###

#highlights and unhighlights vertices/edges when the mouse is moved in 'select path' mode
def select_path_move(event):
	if mode == 'select path':
		global highlighted_object
		if path_to_split == []:		#if no edges/vertices are currently selected in the path
			object = 'vertex'
			var = d.closest_vertex(event.x, event.y, vert_dict, canvas, origin, x_scale, y_scale)	#find if a vertex should be highlighted now
		elif len(path_to_split) == 1:	#if there is just a starting vertex, we may want to be able to unselect it, so could highlight edge or vertex
			v = path_to_split[len(path_to_split)-1]
			var = d.closest_vertex(event.x, event.y, vert_dict, canvas, origin, x_scale, y_scale)
			if var and var[0] == v:
				object = 'vertex'
			else:
				object = 'edge'
				var = d.closest_edge_from(v, event.x, event.y, vert_dict, canvas, origin, x_scale, y_scale) #find if an edge should be
		else:				#otherwise, we're looking for an edge to highlight
			object = 'edge'
			v = path_to_split[len(path_to_split)-1]
			var = d.closest_edge_from(v, event.x, event.y, vert_dict, canvas, origin, x_scale, y_scale)

		if var:			#if something should be highlighted
			if highlighted_object:	#if some object is currently highlighted (must be a vertex in this case)
				if highlighted_object[1] == var[0]:	#if same vertex/edge is still highlighted
					return	#do nothing
				else:
					canvas.delete(highlighted_object[0])	#otherwise, delete current highlighting
			#now add new highlighting and update highlighted_object to v
			if object == 'vertex':
				v, v_coords = var
				highlighted_object = (canvas.create_oval(v_coords[0]-3, v_coords[1]-3, v_coords[0]+3, v_coords[1]+3, fill='pink'), v, v_coords)
			elif object == 'edge':
				w, v_coords, w_coords = var
				highlighted_object = (canvas.create_line(v_coords[0], v_coords[1], w_coords[0], w_coords[1], width=4, fill='pink'), w, w_coords, v_coords)
				canvas.tag_lower(highlighted_object[0])

		else:	#if no object should be highlighted
			if highlighted_object:	#if some object is highlighted
				canvas.delete(highlighted_object[0])	#delete current highlighting
				highlighted_object = False	#and update highlighted_object


#selects and unselects vertices/edges when they are clicked on in 'select path' mode
def select_path_click(event):
	global highlighted_object
	if mode == 'select path' and highlighted_object:	#if in correct mode and there is some highlighted object
		id = highlighted_object[0]
		v = highlighted_object[1]
		v_coords = highlighted_object[2]
		if path_to_split == []:				#if selecting first vertex, delete highlighting and update everything
			canvas.delete(id)
			highlighted_object = False
			path_to_split.append(v)
			selected_objects.append(canvas.create_oval(v_coords[0]-3, v_coords[1]-3, v_coords[0]+3, v_coords[1]+3, fill='magenta'))
		elif len(path_to_split) == 1 and path_to_split[0] == v or path_to_split[len(path_to_split)-2] == v:	#if unselecting most recent selection
			path_to_split.pop()			#keep highlighting but delete selection and update
			canvas.delete(selected_objects.pop())
		else:						#if selecting new edge, keep highlighting and update everything
			w_coords = highlighted_object[3]
			path_to_split.append(v)
			selected_objects.append(canvas.create_line(v_coords[0], v_coords[1], w_coords[0], w_coords[1], width=2, fill='magenta'))


#splits the current path and returns to 'view' mode
def select_path_split(event):
	global mode
	global path_to_split
	global highlighted_object
	global selected_objects
	global history
	global history_index
	if mode == 'select path':
		if len(path_to_split) >= 3:	#if we have a valid path to split
			#mode = 'view'
			calc.split_path(path_to_split, vert_dict, cell_dict, cell_index_dict)

			history_index = history_index + 1
			if len(history) > history_index:
				del history[history_index:]
			history.append((copy.deepcopy(vert_dict), copy.deepcopy(cell_dict), copy.deepcopy(cell_index_dict)))

			canvas.delete('all')
			d.draw_vKd(vert_dict, canvas, origin, x_scale, y_scale)
			d.draw_contour_lines_better(cell_dict, canvas, origin, x_scale, y_scale)
			if highlighted_object:
				highlighted_object = False
			selected_objects = []
			path_to_split = []
		else:
			print('selected path not long enough to split; must be at least two edges long')


#changes the mode based on the keyboard character pressed
def change_mode(event):
	global mode
	global path_to_split
	global highlighted_object
	global selected_objects
	global drag_coords

	new_mode_char = event.char
	if mode == 'select path' and new_mode_char != 'p':
		path_to_split = []
		if highlighted_object:
			canvas.delete(highlighted_object[0])
			highlighted_object = False
		for id in selected_objects:
			canvas.delete(id)
		selected_objects = []
	elif mode == 'edit' and new_mode_char != 'e':
		drag_coords = False

	if new_mode_char == 'v':
		mode = 'view'
	elif new_mode_char == 'p':
		mode = 'select path'
	elif new_mode_char == 'e':
		mode = 'edit'



#### 'edit' mode functions ####

#highlights nearest vertex, if it is close enough to be interacted with
def highlight_vertex(event):
	global highlighted_object
	if mode == 'edit':
		var = d.closest_vertex(event.x, event.y, vert_dict, canvas, origin, x_scale, y_scale)	#find if a vertex should be highlighted now

		if var:	#if vertex should be highlighted
			v, v_coords = var
			if highlighted_object:	#if some vertex is currently highlighted
				if highlighted_object[1] == v:	#if same vertex/edge is still highlighted
					return	#do nothing
				else:
					canvas.delete(highlighted_object[0])	#otherwise, delete current highlighting
			#now add new highlighting and update highlighted_object to v
			highlighted_object = (canvas.create_oval(v_coords[0]-3, v_coords[1]-3, v_coords[0]+3, v_coords[1]+3, fill='pink', tags='vertex'), v, v_coords)

		else:	#if no vertex should be highlighted
			if highlighted_object:	#if some object is highlighted
				canvas.delete(highlighted_object[0])	#delete current highlighting
				highlighted_object = False	#and update highlighted_object


#highlights nearest object (vertex or edge), if it is close enough to be interacted with
def highlight_object(event):
	return

#selects highlighted object with a click.  if holding control, keeps the selection on other selected objects
def select(event):
	return

#unselects selected object(s) with a click. if holding control and object highlighted, only unselects that object
def unselect(event):
	return

#drag a vertex around the canvas, if you're on it
def start_drag(event):
	global drag_coords
	if mode == 'edit' and highlighted_object and 'vertex' in canvas.gettags(highlighted_object[0]):	#if there is a highlighted vertex
		drag_coords = highlighted_object[2]	#set up initial coordinates of the dragged object
		v = highlighted_object[1]
		for w in vert_dict[v]['adj']:	#delete the contour segments of each cell v is on the boundary of, since they may change
			cell_index = cell_index_dict[(v,w)]
			if cell_index != 'boundary':
				for segment_id in cell_dict[cell_index]['segment ids']:
					canvas.delete(segment_id)
		continue_drag(event)	#snap vertex to the mouse

def continue_drag(event):
	global drag_coords
	if mode == 'edit' and drag_coords:
		new_coords = (canvas.canvasx(event.x), canvas.canvasy(event.y))
		x_change = new_coords[0] - drag_coords[0]
		y_change = new_coords[1] - drag_coords[1]

		(highlight_id, v, v_coords) = highlighted_object
		canvas.move(highlight_id, x_change, y_change)			#move highlighting
		canvas.move(vert_dict[v]['id'], x_change, y_change)	#move vertex
		canvas.move(vert_dict[v]['dist id'], x_change, y_change)	#move distance label
		
		for w in vert_dict[v]['adj']:	#move and update info for all the edges to v
			w_coords = (origin[0] + x_scale*vert_dict[w]['coords'][0], origin[1] + y_scale*vert_dict[w]['coords'][1])
			old_edge_id = vert_dict[v]['edge ids'][w]
			canvas.delete(old_edge_id)
			new_edge_id = canvas.create_line(new_coords[0], new_coords[1], w_coords[0], w_coords[1], tags='edge')
			vert_dict[v]['edge ids'][w] = new_edge_id
			vert_dict[w]['edge ids'][v] = new_edge_id

		drag_coords = new_coords	#update drag coordinates for next time

def end_drag(event):
	global drag_coords
	global highlighted_object
	if mode == 'edit' and drag_coords:
		(highlight_id, v, v_coords) = highlighted_object
		vert_dict[v]['coords'] = ((drag_coords[0] - origin[0])/x_scale, (drag_coords[1] - origin[1])/y_scale)	#update position of vertex in vert_dict
		highlighted_object = (highlight_id, v, drag_coords)	#update position of highlighted_object
		drag_coords = False	#no longer dragging
		
		#now need to fix all the 2-cells whose shapes have been changed and re-draw their contour lines
		cells_to_fix = set()
		for w in vert_dict[v]['adj']:
			cell_index = cell_index_dict[(v,w)]
			if cell_index != 'boundary':
				cells_to_fix.add(cell_index)

		calc.fix_cells(cells_to_fix, cell_dict, vert_dict)
		d.draw_contour_lines_better(cell_dict, canvas, origin, x_scale, y_scale, cells_to_fix)


#add a vertex to the graph, at the point clicked
def add_vertex(event):
	return

#start adding a new edge out of a vertex, by double-clicking on it while it's highlighted
def start_edge(event):
	return

#drag the new edge out of the vertex.  include ability to snap endpoint to vertices that could be endpoint of new edge if you get close enough to them
def drag_edge(event):
	return

#finish off edge by letting go.  if it was snapped to an endpoint it stays and updates all the dictionaries.  if not, it disappears
def finish_edge(event):
	return

#deletes highlighted object(s) and updates all the things
def delete_selection(event):
	return



#### Functions for all modes ####

#toggle what is shown on canvas
def toggle_distances(event):
	global distances_showing
	if graph_showing:
		if distances_showing:
			canvas.itemconfigure('distance', state='hidden')
			distances_showing = False
		else:
			canvas.itemconfigure('distance', state='normal')
			distances_showing = True

def toggle_graph(event):
	global graph_showing
	if graph_showing:
		canvas.itemconfigure('vertex', state='hidden')
		canvas.itemconfigure('edge', state='hidden')
		if distances_showing:
			canvas.itemconfigure('distance', state='hidden')
		graph_showing = False
	else:
		canvas.itemconfigure('vertex', state='normal')
		canvas.itemconfigure('edge', state='normal')
		if distances_showing:
			canvas.itemconfigure('distance', state='normal')
		graph_showing = True

def toggle_contour_lines(event):
	global contour_lines_showing
	if contour_lines_showing:
		canvas.itemconfigure('segment', state='hidden')
		contour_lines_showing = False
	else:
		canvas.itemconfigure('segment', state='normal')
		contour_lines_showing = True


#mouse panning

#also used for dragging vertices and endpoints of edges in edit mode
def move_start(event):
	print('move_start entered')
	canvas.scan_mark(event.x, event.y)

def move_move(event):
	global highlighted_object
	if highlighted_object:	#if there is an object highlighted,
		if mode == 'select path':			#if we are selecting a path, unhighlight and pan.  otherwise, don't pan (maybe should still pan if it's an edge, because edges are un-drag-able?)
			canvas.delete(highlighted_object[0])
			highlighted_object = False
			canvas.scan_dragto(event.x, event.y, gain=1)
	else:
		canvas.scan_dragto(event.x, event.y, gain=1)


#keyboard panning
def pan_up(event):
	canvas.scan_mark(0, 0)
	canvas.scan_dragto(0, 30, gain=1)

def pan_down(event):
	canvas.scan_mark(0, 0)
	canvas.scan_dragto(0, -30, gain=1)

def pan_left(event):
	canvas.scan_mark(0, 0)
	canvas.scan_dragto(30, 0, gain=1)

def pan_right(event):
	canvas.scan_mark(0, 0)
	canvas.scan_dragto(-30, 0, gain=1)


#windows zoom
def zoomer(event):
	#print('zoomer entered')
	global x_scale
	global y_scale
	global origin
	canvas_x = canvas.canvasx(event.x)
	canvas_y = canvas.canvasy(event.y)
	if (event.delta > 0):
		canvas.scale("all", canvas_x, canvas_y, 1.1, 1.1)
		x_scale = x_scale*1.1
		y_scale = y_scale*1.1
		origin = g.add(origin, g.mult(0.1, g.vector((canvas_x, canvas_y), origin)))
	elif (event.delta < 0):
		canvas.scale("all", canvas_x, canvas_y, 10/11, 10/11)
		x_scale = x_scale*10/11
		y_scale = y_scale*10/11
		origin = g.add(origin, g.mult(-1/11, g.vector((canvas_x, canvas_y), origin)))
	canvas.configure(scrollregion = canvas.bbox("all"))


#linux zoom
def zoomerP(event):
	global x_scale
	global y_scale
	global origin
	canvas_x = canvas.canvasx(event.x)
	canvas_y = canvas.canvasy(event.y)
	canvas.scale("all", canvas_x, canvas_y, 1.1, 1.1)
	x_scale = x_scale*1.1
	y_scale = y_scale*1.1
	origin = g.add(origin, g.mult(0.1, g.vector((canvas_x, canvas_y), origin)))
	canvas.configure(scrollregion = canvas.bbox("all"))

def zoomerM(event):
	global x_scale
	global y_scale
	global origin
	canvas_x = canvas.canvasx(event.x)
	canvas_y = canvas.canvasy(event.y)
	canvas.scale("all", canvas_x, canvas_y, 10/11, 10/11)
	x_scale = x_scale*10/11
	y_scale = y_scale*10/11
	origin = g.add(origin, g.mult(-1/11, g.vector((canvas_x, canvas_y), origin)))
	canvas.configure(scrollregion = canvas.bbox("all"))


#keyboard zoom
def zoomerIn(event):
	global x_scale
	global y_scale
	global origin
	#print('zoomerDown performed')
	canvas_x = canvas.canvasx(win_width/2)
	canvas_y = canvas.canvasy(win_height/2)
	canvas.scale("all", canvas_x, canvas_y, 1.1, 1.1)
	x_scale = x_scale*1.1
	y_scale = y_scale*1.1
	origin = g.add(origin, g.mult(0.1, g.vector((canvas_x, canvas_y), origin)))
	canvas.configure(scrollregion = canvas.bbox("all"))

def zoomerOut(event):
	global x_scale
	global y_scale
	global origin
	#print('zoomerUp performed')
	canvas_x = canvas.canvasx(win_width/2)
	canvas_y = canvas.canvasy(win_height/2)
	canvas.scale("all", canvas_x, canvas_y, 10/11, 10/11)
	x_scale = x_scale*10/11
	y_scale = y_scale*10/11
	origin = g.add(origin, g.mult(-1/11, g.vector((canvas_x, canvas_y), origin)))
	canvas.configure(scrollregion = canvas.bbox("all"))

#keyboard zoom in just one direction
def zoomerInX(event):
	global x_scale
	global origin
	canvas_x = canvas.canvasx(win_width/2)
	canvas_y = canvas.canvasy(win_height/2)
	canvas.scale("all", canvas_x, canvas_y, 1.1, 1)
	x_scale = x_scale*1.1
	origin = (origin[0] + 0.1*(origin[0] - canvas_x), origin[1])
	canvas.configure(scrollregion = canvas.bbox("all"))

def zoomerOutX(event):
	global x_scale
	global origin
	canvas_x = canvas.canvasx(win_width/2)
	canvas_y = canvas.canvasy(win_height/2)
	canvas.scale("all", canvas_x, canvas_y, 10/11, 1)
	x_scale = x_scale*10/11
	origin = (origin[0] - (origin[0] - canvas_x)/11, origin[1])
	canvas.configure(scrollregion = canvas.bbox("all"))

def zoomerInY(event):
	global y_scale
	global origin
	canvas_x = canvas.canvasx(win_width/2)
	canvas_y = canvas.canvasy(win_height/2)
	canvas.scale("all", canvas_x, canvas_y, 1, 1.1)
	y_scale = y_scale*1.1
	origin = (origin[0], origin[1] + 0.1*(origin[1] - canvas_y))
	canvas.configure(scrollregion = canvas.bbox("all"))

def zoomerOutY(event):
	global y_scale
	global origin
	canvas_x = canvas.canvasx(win_width/2)
	canvas_y = canvas.canvasy(win_height/2)
	canvas.scale("all", canvas_x, canvas_y, 1, 10/11)
	y_scale = y_scale*10/11
	origin = (origin[0], origin[1] - (origin[1] - canvas_y)/11)
	canvas.configure(scrollregion = canvas.bbox("all"))


#undo/redo functionality (currently only for paths that have been split, and only in view mode)
def undo(event):
	global history_index
	global vert_dict
	global cell_dict
	global cell_index_dict
	global highlighted_object
	global path_to_split
	global selected_objects
	if history_index > 0: #and mode == 'view'?
		history_index = history_index - 1
		vert_dict, cell_dict, cell_index_dict = history[history_index]
		vert_dict = copy.deepcopy(vert_dict)
		cell_dict = copy.deepcopy(cell_dict)
		cell_index_dict = copy.deepcopy(cell_index_dict)

		canvas.delete('all')
		d.draw_vKd(vert_dict, canvas, origin, x_scale, y_scale)
		d.draw_contour_lines_better(cell_dict, canvas, origin, x_scale, y_scale)

		if not graph_showing:
			canvas.itemconfigure('vertex', state='hidden')
			canvas.itemconfigure('edge', state='hidden')
			canvas.itemconfigure('distance', state='hidden')
		elif not distances_showing:
			canvas.itemconfigure('distance', state='hidden')
		if not contour_lines_showing:
			canvas.itemconfigure('segment', state='hidden')

		if mode == 'select path':
			path_to_split = []
			if highlighted_object:
				highlighted_object = False
			selected_objects = []


def redo(event):
	global history_index
	global vert_dict
	global cell_dict
	global cell_index_dict
	if len(history) - 1 > history_index: #and mode == 'view'?
		history_index = history_index + 1
		vert_dict, cell_dict, cell_index_dict = history[history_index]
		vert_dict = copy.deepcopy(vert_dict)
		cell_dict = copy.deepcopy(cell_dict)
		cell_index_dict = copy.deepcopy(cell_index_dict)

		canvas.delete('all')
		d.draw_vKd(vert_dict, canvas, origin, x_scale, y_scale)
		d.draw_contour_lines_better(cell_dict, canvas, origin, x_scale, y_scale)

		if not graph_showing:
			canvas.itemconfigure('vertex', state='hidden')
			canvas.itemconfigure('edge', state='hidden')
			canvas.itemconfigure('distance', state='hidden')
		elif not distances_showing:
			canvas.itemconfigure('distance', state='hidden')
		if not contour_lines_showing:
			canvas.itemconfigure('segment', state='hidden')

		if mode == 'select path':
			path_to_split = []
			if highlighted_object:
				highlighted_object = False
			selected_objects = []
	

#distinguish between panning canvas and dragging vertex
def click(event):
	if mode == 'edit' and highlighted_object and 'vertex' in canvas.gettags(highlighted_object[0]):
		start_drag(event)
	else:
		move_start(event)

def drag(event):
	if mode == 'edit' and drag_coords:
		continue_drag(event)
	else:
		move_move(event)

#distinguish between edit and select path mode functions
def motion(event):
	if mode == 'edit':
		highlight_vertex(event)
	elif mode == 'select path':
		select_path_move(event)

def release(event):
	if mode == 'edit':
		end_drag(event)
	elif mode == 'select path':
		select_path_click(event)
 


#----- Bind events and functions  -----

#change mode
canvas.bind('v', change_mode)
canvas.bind('p', change_mode)
canvas.bind('e', change_mode)

e = Event()

#'select path' mode:
#highlighting
canvas.bind('<Motion>', motion)	#also for highlighting in edit mode
#selecting
canvas.bind('<ButtonRelease-1>', release)	#also for ending dragging a vertex in edit mode
#splitting
canvas.bind('s', select_path_split)	#bind this to something better

#'edit' mode:
#highlighting vertices
#canvas.bind('<Motion>', highlight_vertex)
#drag vertices
canvas.bind('<ButtonPress-1>', click)		#also for panning canvas with mouse
canvas.bind('<B1-Motion>', drag)		#also for panning canvas with mouse
#canvas.bind('<ButtonRelease-1>', end_drag)

#showing and hiding different parts of the picture:
canvas.bind('d', toggle_distances)
canvas.bind('g', toggle_graph)
canvas.bind('l', toggle_contour_lines)

#Panning canvas using click and drag:
#canvas.bind("<ButtonPress-1>", move_start)
#canvas.bind("<B1-Motion>", move_move)
#Panning canvas with keyboard:
canvas.bind("<Up>", pan_up)
canvas.bind("<Down>", pan_down)
canvas.bind("<Left>", pan_left)
canvas.bind("<Right>", pan_right)
#linux scroll:
canvas.bind("<Button-4>", zoomerP)
canvas.bind("<Button-5>", zoomerM)
#windows scroll:
canvas.bind("<MouseWheel>",zoomer)
#keyboard scroll:
canvas.bind("<Shift-Up>", zoomerOut)
canvas.bind("<Shift-Down>", zoomerIn)
canvas.bind("<Control-Left>", zoomerInX)
canvas.bind("<Control-Right>", zoomerOutX)
canvas.bind("<Control-Down>", zoomerInY)
canvas.bind("<Control-Up>", zoomerOutY)

#undo/redo
canvas.bind('u', undo)
canvas.bind('r', redo)


canvas.focus_set()
window.mainloop()