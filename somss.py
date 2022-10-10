#!/usr/bin/env python3
"""     20.07.2022
	Folkenberg, Siro
        Nuclear Engineering and Thermal Physics, National Research Nuclear University-IATE, Russia.
	Geospatial (Spaceborne Optical Remote Sensing), Technical University of Kenya

      	================== SIMPLE ORBITAL MECHANICAL SYSTEM SIMULATION  =================================
						SOMSS Modules

	NB: All physical quantities in SI units
"""
import numpy
import matplotlib.pyplot
import tkinter
from tkinter import *
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

root = tkinter.Tk()
root.title("SOMSS")
root['bg'] = 'dark gray'
root.geometry('1360x800')
root.resizable(False, False)
wi = root.winfo_screenwidth()
he = root.winfo_screenheight()

newFile_image = tkinter.PhotoImage(file = './icons/new_file.gif')
openFile_image = tkinter.PhotoImage(file = './icons/open_file.gif')
saveFile_image = tkinter.PhotoImage(file = './icons/save_file.gif')


iframe = tkinter.Frame(root, bd = 2)
iframe.grid(row = 0, columnspan = 2, sticky='nesw')
iframe1 = tkinter.Frame(root, bg='dark gray')
iframe1.grid(row = 2, column = 0, sticky='nsew', pady=10, padx=20)
iframe2 = tkinter.Frame(root, bd=2, relief='sunken', borderwidth=3, bg='dark gray')
iframe2.grid(row = 2, column = 1, sticky='news', pady=10)
iframe3 = tkinter.Frame(root, bd=2, relief='sunken', borderwidth=3, bg = 'black')
iframe3.grid(row = 2, column = 2, sticky='news', pady=10)

menuBar = tkinter.Menubutton(iframe, text = 'File', font = ('Times New Roman', 9), relief='raised', width=8, padx=2, borderwidth=0)
menuBar.pack(side = 'left')
fileMenu = tkinter.Menu(menuBar, tearoff = 0)
menuBar['menu'] = fileMenu
fileMenu.add('command', label = 'New', compound = 'left', accelerator = 'Ctrl+N', underline = 0, image = newFile_image, font = ('Helvetica', 9))
fileMenu.add('command', label = 'Open', compound = 'left', accelerator = 'Ctrl+O', underline = 0, image = openFile_image, font = ('Helvetica', 9))
fileMenu.add('command', label = 'Save', compound = 'left', accelerator = 'Ctrl+S', underline = 0, image = saveFile_image, font = ('Helvetica', 9))
fileMenu.add('command', label = 'Save As', accelerator = 'Shift+Ctrls+S', underline = 0, image = saveFile_image, compound = 'left', font = ('Helvetica', 9))
fileMenu.add_separator()
fileMenu.add('command', label = 'Exit', underline = 0, compound = 'left', accelerator = 'Alt+F4', font = ('Helvetica', 9))

labelframe1 = tkinter.LabelFrame(iframe1, text="Status Window", bg='dark gray', fg='white', font = ('Helvetica bold', 9))
labelframe1.pack(fill="both", expand="yes")

canvFrame = tkinter.Frame(iframe2, bg = 'dark gray')
canvFrame.grid(row = 1, column = 1, sticky='e', padx=5, pady=5)
wn = tkinter.Canvas(canvFrame, width = 700, height = 640, bg = "black")
wn.pack(side='right')
wn1 = tkinter.Canvas(labelframe1, bg = "black", height = 170, width=340)
wn1.pack(fill="both", padx=5)

wnext = tkinter.Canvas(labelframe1, bg = "black", height = 130, width=340)
wnext.pack(fill="both", padx=5)

labelframe = tkinter.LabelFrame(iframe1, text="Extended Output", fg='white', bg='dark gray', font = ('Helvetica bold', 9))
labelframe.pack(fill="both", expand="yes", pady=10)

wn2 = tkinter.Text(labelframe, bg = "khaki", fg="blue", height=13, width=34, border=0)
wn2.grid(row=0, columnspan=4, sticky = 'snwe', padx=5, pady=5)
other_orbit = "\n>>> e = 0: A circular orbit; geosynchronous orbit; altitude 36,000 km((42164 km from the center of the Earth ); fairly CONST orbital vel 3074.6 m/s\n>>> 0<e<1: Elliptical orbit; perigee - orbital vel largest val; apogee - orbital vel lowest   "
unstable_orbit = "\n>>> Initial orbital vel < 1886 m/s: SAT-X inertia can no longer effectively counteract the grav due to the central mass. This velocity is not enough to propel SAT-X away from the central mass.\n>>> Grav pulls SAT-X towards the central mass, ultimately leading to an unstable orbit\n"
wn2.insert('1.0', 'METADATA\n\n')

class StayInOrbit:
	def __init__(self, satellite, earth, P_i, V_i, tm, tdelt):
		self.P_ini = P_i
		self.V_ini = V_i
		self.end = tm
		self.tdelt = tdelt
		self.e = earth
		self.sat = satellite
		self.G = 6.67408e-11
	def tm(self, start = 0):
		return numpy.linspace(start, self.end, int(self.end / self.tdelt) + 1)
	def Vx(self, P_x, P_y):
		return -(self.G * self.e * P_x) / ((P_x ** 2 + P_y ** 2) ** (3 / 2))
	def Vy(self, P_x, P_y):
		return -(self.G * self.e * P_y) / ((P_x ** 2 + P_y ** 2) ** (3 / 2))

	''' Runge - Kutta 4th order in solving
	the required derivatives '''
	def RK_xyv(self, P_x, P_y, V_x, V_y):
		RK1_x = V_x
		RK1_y = V_y
		RK1_Vx = self.Vx(P_x, P_y)
		RK1_Vy = self.Vy(P_x, P_y)
		RK2_x = V_x + (self.tdelt * RK1_Vx) / 2
		RK2_y = V_y + (self.tdelt * RK1_Vy) / 2
		RK2_Vx = self.Vx(P_x + (self.tdelt * RK1_x) / 2, P_y + (self.tdelt * RK1_y) / 2)
		RK2_Vy = self.Vy(P_x + (self.tdelt * RK1_x) / 2, P_y + (self.tdelt * RK1_y) / 2)
		RK3_x = V_x + (self.tdelt * RK2_Vx) / 2
		RK3_y = V_y + (self.tdelt * RK2_Vy) / 2
		RK3_Vx = self.Vx(P_x + (self.tdelt * RK2_x) / 2, P_y + (self.tdelt * RK2_y) / 2)
		RK3_Vy = self.Vy(P_x + (self.tdelt * RK2_x) / 2, P_y + (self.tdelt * RK2_y) / 2)
		RK4_x = V_x + self.tdelt * RK3_Vx
		RK4_y = V_y + self.tdelt * RK3_Vy
		RK4_Vx = self.Vx(P_x + self.tdelt * RK3_x, P_y + self.tdelt * RK3_y)
		RK4_Vy = self.Vy(P_x + self.tdelt * RK3_x, P_y + self.tdelt * RK3_y)
		return ([RK1_x, RK2_x, RK3_x, RK4_x], [RK1_y, RK2_y, RK3_y, RK4_y], [RK1_Vx, RK2_Vx, RK3_Vx, RK4_Vx], [RK1_Vy, RK2_Vy, RK3_Vy, RK4_Vy])
	def xy_step(self, P_x, P_y, RKx_dat, RKy_dat):
		P_x = P_x + (self.tdelt / 6) * (RKx_dat[0] + 2 * RKx_dat[1] + 2 * RKx_dat[2] + RKx_dat[3])
		P_y = P_y + (self.tdelt / 6) * (RKy_dat[0] + 2 * RKy_dat[1] + 2 * RKy_dat[2] + RKy_dat[3])
		return (P_x, P_y)
	def xyv_step(self, V_x, V_y, RK_Vx_dat, RK_Vy_dat):
		V_x = V_x + (self.tdelt / 6) * (RK_Vx_dat[0] + 2 * RK_Vx_dat[1] + 2 * RK_Vx_dat[2] + RK_Vx_dat[3])
		V_y = V_y + (self.tdelt / 6) * (RK_Vy_dat[0] + 2 * RK_Vy_dat[1] + 2 * RK_Vy_dat[2] + RK_Vy_dat[3])
		return (V_x, V_y)

	def run_sim(self):
		P_x, P_y = self.P_ini[0], self.P_ini[1]
		V_x, V_y = self.V_ini[0], self.V_ini[1]
		tdat = self.tm()
		Rx, Ry = numpy.zeros(len(tdat)), numpy.zeros(len(tdat))
		Rdx, Rdy = numpy.zeros(len(tdat)), numpy.zeros(len(tdat))
		for i in range(0, len(tdat)):
			RKx_dat, RKy_dat, RK_Vx_dat, RK_Vy_dat = self.RK_xyv(P_x, P_y, V_x, V_y)
			P_x, P_y = self.xy_step(P_x, P_y, RKx_dat, RKy_dat)
			V_x, V_y = self.xyv_step(V_x, V_y, RK_Vx_dat, RK_Vy_dat)
			Rx[i], Ry[i], Rdx[i], Rdy[i] = P_x, P_y, V_x, V_y
		return (Rx, Ry, Rdx, Rdy)

	''' Run animation  '''
	def run_anim(Rx, Ry, Rdx, Rdy):
		fig = matplotlib.pyplot.figure(figsize=(3,3))
		ax = matplotlib.pyplot.subplot()
		ax.grid(True)
		ax.set_facecolor('k')
		ax.grid(which = 'minor', linewidth = 0.5, linestyle=":", color='#EEEEEE')
		ax.grid(which = 'major', linewidth = 0.7)
		ax.minorticks_on()
		w = 700;h = 640;Xc = w/2;Yc = h/2
		Xs = (w - 100)/(2 * Rx.max() - Rx.min())
		Ys = (h - 100) / (2 * Ry.max() - Ry.min())
		Rx *= Xs;Ry *= -Ys;Rx += Xc;Ry += Yc
		earth_tkrad = 22;sat_tkrad = 10
		x0_earth, y0_earth, x_earth, y_earth = Xc - earth_tkrad, Yc - earth_tkrad, Xc + earth_tkrad, Yc + earth_tkrad
		x0_sat, y0_sat, x_sat, y_sat = Xc, Yc, Xc, Yc
		e_ob = wn.create_oval(x0_earth, y0_earth, x_earth, y_earth, fill = "blue", outline = "grey")
		s_ob = wn.create_oval(x0_sat, y0_sat, x_sat, y_sat, fill = "blue", outline = "grey")
		r_vect = wn.create_line(Xc, Yc, Xc, Yc, fill = "")
		disp_dat = wn1.create_rectangle(10, 150, 320, 40, fill = "black", outline = "blue")
		label1 = wn.create_text(Xc + 10, Yc + 10, text = "Earth", fill = "white", anchor = "nw", font = ('Helvetica bold', 8))
		label2 = wn.create_text(Rx[0] + 10, Ry[0] + 10, text = "SAT-X", fill = "white", anchor = "s", font = ('Helvetica bold', 5))
		sat_pos = []
		for i in range(0, len(Rdx)):
			sat_pos.append(Rx[i])
			sat_pos.append(Ry[i])
		sat_path2 = wn.create_polygon(sat_pos, outline = "red", dash = "..", fill = "")
		val = 0
		while val <= len(Rdx) + 1:
			if val == len(Rdx):
				val = 0
			wn.coords(s_ob, Rx[val] - sat_tkrad, Ry[val] - sat_tkrad, Rx[val] + sat_tkrad, Ry[val] + sat_tkrad)
			wn.coords(label2, Rx[val], Ry[val] - 10)
			wn.coords(r_vect, Xc, Yc, Rx[val], Ry[val])
			earth_sat_dist = numpy.sqrt(((Xc - Rx[val]) ** 2) + ((Yc - Ry[val]) ** 2)) * 2 / (Xs + Ys)
			if earth_sat_dist < 6378000:
				print("Satellite Orbit Unstable!")
				wn.delete(label1)
				wn.delete(label2)
				wn1.delete('all')
				wn1.create_text(15, 30, text = "Realtime Orbital Parameters: ", fill='blue', anchor='sw', font = ('Helvetica bold', 8))
				wn1.create_text(15, 100, text="Satellite Orbit Unstable!", fill="red", anchor="sw", font = ('Helvetica bold', 8))
				wn2.delete("3.0","end")
				wn2.insert('3.0', unstable_orbit)
				break
			else:
				sat_vel = numpy.sqrt(Rdx[val] ** 2 + Rdy[val] ** 2)
				sat_vel, earth_sat_dist = "%.5f" % sat_vel, "%.5f" % earth_sat_dist
				head_1 = wn1.create_text(15, 30, text = "Realtime Orbital Parameters: ", fill='blue', anchor='sw', font = ('Helvetica bold', 8))
				rad_mag = wn1.create_text(15, 80, text = "Orbital Radius: ", fill = "blue", anchor = "sw", font = ('Helvetica bold', 8))
				rad_mag1 = wn1.create_text(140, 80, text = " " + earth_sat_dist + " m", anchor = "sw", fill = "green", font = ('Helvetica bold', 8))
				V_t = wn1.create_text(15, 120, text = "Orbital Velocity: ", fill = "blue", anchor = "sw", font = ('Helvetica bold', 8))
				V_t1 = wn1.create_text(140, 120, text = " " + sat_vel + " m/s", fill = "green", anchor = "sw", font = ('Helvetica bold', 8))
				wn1.after(1, wn1.update())
				wn1.delete(rad_mag)
				wn1.delete(rad_mag1)
				wn1.delete(V_t)
				wn1.delete(V_t1)
				wn1.delete(head_1)
			val += 1
menuBar_3 = tkinter.Menubutton(iframe, text = 'Hypothetical Orbital System', font = ('Times New Roman', 9), relief='raised', width=30, borderwidth=0)
menuBar_3.pack(side = 'left')
orbitMenu = tkinter.Menu(menuBar_3, tearoff = 0)
menuBar_3['menu'] = orbitMenu

def orbital_sys():
	two_body_transient_toplevel = tkinter.Toplevel(root)
	progName = 'ORBITAL SYSTEM'
	two_body_transient_toplevel.title(progName)
	wWidth = two_body_transient_toplevel.winfo_reqwidth()
	wHeight = two_body_transient_toplevel.winfo_reqheight()
	pRight = int(two_body_transient_toplevel.winfo_screenwidth()/2 - wWidth/0.76)
	pDown = int(two_body_transient_toplevel.winfo_screenheight()/2 - wHeight/0.74)
	two_body_transient_toplevel.geometry('+{}+{}'.format(pRight, pDown))
	two_body_transient_toplevel.resizable(False, False)

	orbit_tabcont = tkinter.ttk.Notebook(two_body_transient_toplevel)
	orbit_tab1 = tkinter.Frame(orbit_tabcont, bg = 'dark gray')
	orbit_tab2 = tkinter.Frame(orbit_tabcont, bg = 'dark gray')
	orbit_tabcont.add(orbit_tab1, text = '2 Body System')
	frame = tkinter.LabelFrame(orbit_tab1, bg = 'dark gray', fg = 'black', text = 'INPUT PARAMETERS',\
	font = ('Helvetica bold', 8))
	frame.grid(row = 1, column = 0, sticky = 'sewn', padx = 15, pady = 15)
	tkinter.Label(frame, text = 'Orbital Period(s):', fg = 'khaki', bg = 'dark gray',\
	font = ('Helvetica', 8)).grid(row = 2, column = 0, sticky = 'w', padx = 5, pady = 5)
	tkinter.Label(frame, text = 'Orbital Velocity(m/s):', fg = 'khaki', bg = 'dark gray',\
	font = ('Helvetica', 8)).grid(row = 3, column = 0, sticky = 'w', padx = 5, pady = 5)
	tkinter.Label(frame, text = 'Orbiting Mass(kg):', fg = 'khaki', bg = 'dark gray',\
	font = ('Helvetica', 8)).grid(row = 4, column = 0, sticky = 'w', padx = 5, pady = 5)
	tkinter.Label(frame, text = 'Orbital Radius(s):', fg = 'khaki', bg = 'dark gray',\
	font = ('Helvetica', 8)).grid(row = 5, column = 0, sticky = 'w', padx = 5, pady = 5)

	var = StringVar()
	var_1 = StringVar()
	val= StringVar()
	val1=StringVar()
	val2=StringVar()
	val3=StringVar()
	e5 = tkinter.Entry(frame, textvariable = var_1, fg = 'blue', bg = 'khaki')
	e5.grid(row = 1, column = 0, columnspan = 7, padx = 5, sticky = 'we', pady = 5)
	e = tkinter.Entry(frame, fg = 'blue', bg = 'khaki', textvariable=val1)
	e.grid(row = 2, column = 1, columnspan = 7, sticky = 'we', padx = 5, pady = 5)
	e1 = tkinter.Entry(frame, fg = 'blue', bg = 'khaki', textvariable=val3)
	e1.grid(row = 3, column = 1, columnspan = 7, sticky = 'we', padx = 5, pady = 5)
	e2 = tkinter.Entry(frame, fg = 'blue', bg = 'khaki', textvariable=val2)
	e2.grid(row = 4, column = 1, columnspan = 7, sticky = 'we', padx = 5, pady = 5)
	e4 = tkinter.Entry(frame, fg = 'blue', bg = 'khaki', textvariable=val)
	e4.grid(row = 5, column = 1, columnspan = 7, sticky = 'we', padx = 5, pady = 5)
	e3 = tkinter.Entry(frame, fg = 'blue', textvariable = var, bg = 'khaki')
	e3.grid(row = 6, column = 0, columnspan = 7, sticky = 'we', padx = 5, pady = 5)

	button = tkinter.Button(frame, text = 'SAVE AS', fg = 'khaki', bg = 'gray', font = ('Helvetica', 7),\
	command = lambda: var.set(filedialog.asksaveasfilename(initialdir = '/', title = 'Select Folder and Define a data File Name',\
	filetypes = (('txt files','*.txt'), ('all files', '*.*')))))
	button.grid(row = 6, column = 7, columnspan = 1, sticky = 'we', padx = 5, pady = 5)

	button = tkinter.Button(frame, text = 'BROWSE', fg = 'khaki', bg = 'gray', font = ('Helvetica', 7),\
	command = lambda: var_1.set(filedialog.askopenfilename(initialdir = '/', title = 'Select MMSLAB text file',\
	filetypes = (('MMSLAB','*.txt'), ('all files','*.*')))))
	button.grid(row = 1, column = 7, columnspan = 1, sticky = 'we', padx = 5, pady = 5)
	button = tkinter. Button(frame, text = 'QUIT', fg = 'khaki', bg = 'gray', font = ('Helvetica', 7), command = two_body_transient_toplevel.destroy)
	button.grid(row = 7 , column = 0, sticky = 'w', padx = 5, pady = 5)

	''' main canvas  '''
	def mainCanv():
		wn1.delete('all')
		wn.delete('all')
		wnext.delete('all')
		wn2.delete("3.0","end")
		wn2.insert('3.0', other_orbit)
		R_orbital = float(val.get()) #4016400
		period_secs = float(val1.get()) #24*3600 Orbital period
		m_sat = float(val2.get()) #4990 satellite mass
		v_sat = float(val3.get()) #3074.6  #v_sat = 1773.9 orbital velocity
		period_factor = 3
		tdelt = 50
		M_earth = 5.9722e24
		ugrav = 6.67408e-11
		ugrav = str(ugrav)
		per = str(period_secs)
		radi = str(R_orbital)
		ms_e = str(M_earth)
		ms_s = str(m_sat)
		factor_p = str(period_factor)
		wnext.create_text(15, 20, text = "Predefined Orbital Parameters: ", fill = "gray", anchor = "sw", font = ('Helvetica bold', 8))
		wnext.create_text(15, 40, text = "Universal Gravitation: ", fill = 'blue', anchor = 'sw', font = ('Helvetica bold', 8))
		wnext.create_text(170, 40, text =' ' +ugrav+  " Nm\u00B2 \kg\u00B2", anchor = "sw", fill = "green", font = ('Helvetica bold', 8))
		wnext.create_text(15, 60, text = "Period: ", fill = 'blue', anchor = 'sw', font = ('Helvetica bold', 8))
		wnext.create_text(170, 60, text =' ' +per+  " s", anchor = "sw", fill = "green", font = ('Helvetica bold', 8))
		wnext.create_text(15, 80, text = 'Mass of Earth: ', fill = 'blue', anchor = 'sw', font = ('Helvetica bold', 8))
		wnext.create_text(170, 80, text =' ' +ms_e+  " kg", anchor = "sw", fill = "green", font = ('Helvetica bold', 8))
		wnext.create_text(15, 100, text = 'Satellite Mass: ', fill = 'blue', anchor = 'sw', font = ('Helvetica bold', 8))
		wnext.create_text(170, 100, text =' ' +ms_s+  " kg", anchor = "sw", fill = "green", font = ('Helvetica bold', 8))
		wnext.create_text(15, 120, text = 'Period Factor: ', fill = 'blue', anchor = 'sw', font = ('Helvetica bold', 8))
		wnext.create_text(170, 120, text =' ' +factor_p+ '' , anchor = "sw", fill = "green", font = ('Helvetica bold', 8))
		s_ob = StayInOrbit(m_sat, M_earth, [0, R_orbital], [v_sat, 0], (period_factor * period_secs), (tdelt))
		Rx, Ry, Rdx, Rdy = s_ob.run_sim()
		StayInOrbit.run_anim(Rx, Ry, Rdx, Rdy)
	def combo():
		mainCanv()

	button = tkinter.Button(frame, text = 'START', fg = 'khaki', bg = 'gray', font = ('Helvetica', 7), command = combo)
	button.grid(row = 7, column = 7, sticky = 'we', padx = 5, pady = 5)


	orbit_tabcont.add(orbit_tab2, text = 'n Body System')
	orbit_tabcont.grid(row = 0, column = 0)
	two_body_transient_toplevel.transient(root)
orbitMenu.add('command', label = '2D Rendition', compound = 'left', underline = 0, font = ('Helvetica', 9), command = orbital_sys)
orbitMenu.add('command', label = '3D Simple Rendition', compound = 'left', underline = 0, font = ('Helvetica', 9), command = "")

menuBar_3 = tkinter.Menubutton(iframe, text = 'On The Fly Python Shell', font = ('Times New Roman', 9), relief='raised', width=20, borderwidth=0)
menuBar_3.pack(side = 'left')
orbitMenu = tkinter.Menu(menuBar_3, tearoff = 0)
menuBar_3['menu'] = orbitMenu

''' a simple command line for chaging orbital parameters on-the-fly(incomplete)  '''
def shell():
        two_body_transient_toplevel = tkinter.Toplevel(root)
        progName = 'INTERACTIVE SHELL'
        two_body_transient_toplevel.title(progName)
        wWidth = two_body_transient_toplevel.winfo_reqwidth()
        wHeight = two_body_transient_toplevel.winfo_reqheight()
        pRight = int(two_body_transient_toplevel.winfo_screenwidth()/3.39)
        pDown = int(two_body_transient_toplevel.winfo_screenheight()/2 - wHeight/0.720)
        two_body_transient_toplevel.geometry('+{}+{}'.format(pRight, pDown))
        two_body_transient_toplevel.resizable(False, False)

        orbit_tabcont = tkinter.ttk.Notebook(two_body_transient_toplevel)
        orbit_tab1 = tkinter.Frame(orbit_tabcont, bg = 'dark gray')
        orbit_tab2 = tkinter.Frame(orbit_tabcont, bg = 'dark gray')
        orbit_tabcont.add(orbit_tab1, text = '       ')
        frame = tkinter.LabelFrame(orbit_tab1, bg = 'dark gray', fg = 'white', text = 'Interactive Shell',\
        font = ('Helvetica bold', 8))
        frame.pack(fill="both", expand="yes",padx=9,  pady=12)


        entry = tkinter.Text(frame, fg = 'green', bg = 'khaki', width = 68, height = 18.5, border = 0)
        entry.grid(row = 0, columnspan=4, sticky = 'snwe', padx=5, pady=5)
        shell_spit = tkinter.Text(frame, fg = 'white', bg = 'black', width = 68, height = 8, border = 0)
        shell_spit.grid(row = 1, columnspan = 4, sticky = 'nswe', padx = 5, pady = 5)

        button = tkinter.Button(frame, text = 'RUN', fg = 'khaki', bg = 'gray', font = ('Helvetica', 7), command = "")
        button.grid(row = 2, column = 3, sticky = 'e' )
        #tkinter.Label(frame, text = 'Raw data:', fg = 'khaki', bg = 'dark gray', font = ('Helvetica', 8),\
        #padx = 5).grid(row = 2, column = 0, sticky = 'w')
        var = tkinter.StringVar()
        entry = tkinter.Entry(frame, textvariable = var, fg = 'blue', bg = 'khaki')
        entry.grid(row = 2, column = 1, sticky = 'w')
        button = tkinter.Button(frame, text = 'BROWSE', fg = 'khaki', bg = 'gray', font = ('Helvetica', 7),\
        command = lambda: var.set(filedialog.askopenfilename(initialdir = '/', title = 'Select file',\
        filetypes = (('data file','*.txt'), ('all files','*.*')))))
        button.grid(row = 2, column = 0, sticky = 'w')


        orbit_tabcont.grid(row = 0, column = 0)
        two_body_transient_toplevel.transient(root)
orbitMenu.add('command', label = 'Interactive Shell', compound = 'left', underline = 0, font = ('Helvetica', 9), command = shell)

''' regression and correlation analysis widget(incomplete)  '''
def regres():
	regres_transient_toplevel = tkinter.Toplevel(root)
	progName = 'REGRESSION & CORRELATION ANALYSIS - VECTOR DATA'
	regres_transient_toplevel.title(progName)
	wWidth = regres_transient_toplevel.winfo_reqwidth()
	wHeight = regres_transient_toplevel.winfo_reqheight()
	pRight = int(regres_transient_toplevel.winfo_screenwidth()/2 - wWidth/0.76)
	pDown = int(regres_transient_toplevel.winfo_screenheight()/2 - wHeight/0.74 )
	regres_transient_toplevel.geometry('+{}+{}'.format(pRight, pDown))
	regres_transient_toplevel.resizable(False, False)

	isolated_frame = tkinter.Frame(regres_transient_toplevel)
	isolated_frame.grid(row = 0, column = 0, sticky = 'we')
	regres_frame = tkinter.Frame(regres_transient_toplevel)
	regres_frame.grid(row = 1, column = 0, sticky = 'sewn')
	regres_tabcont = tkinter.ttk.Notebook(regres_frame)
	regres_tab1 = tkinter.Frame(regres_transient_toplevel, bg = 'dark gray')
	regres_tab2 = tkinter.Frame(regres_transient_toplevel, bg = 'dark gray')
	regres_tab3 = tkinter.Frame(regres_transient_toplevel, bg = 'dark gray')
	regres_tab4 = tkinter.Frame(regres_transient_toplevel, bg = 'dark gray')

	regres_tabcont.add(regres_tab1, text = 'Linear Fitting')
	frame = tkinter.Frame(regres_tab1, relief = 'ridge', bg = 'dark gray')
	frame.grid(row = 1, column = 0, sticky = 'w', padx = 10, pady = 10)
	prev_frame = tkinter.LabelFrame(regres_tab1, relief = 'ridge', bg = 'dark gray', fg = 'black', text = 'PREVIEW', font = ('Helvetica bold', 8), width = 50)
	prev_frame.grid(row = 0, column = 1, sticky = 'sewn', padx = 15, pady = 15)
	var = StringVar()
	e6 = tkinter.Entry(isolated_frame, textvariable = "", bg = 'khaki', fg = 'blue', width = 40)
	e6.grid(row = 1, column = 0, columnspan = 7, sticky = 'sewn', pady = 5, padx = 5)
	button = tkinter.Button(isolated_frame, text = 'MMSLAB.txt', fg = 'khaki', bg = 'gray', font = ('Helvetica', 7), \
	command = lambda: var.set(filedialog.asksaveasfilename(initialdir = '/', title = 'Select Folder and Define MMSLAB file name',\
	filetypes = (('txt files','*.txt'), ('all files', '*.*')))))
	button.grid(row = 1, column = 7, columnspan = 1, sticky = 'we', padx = 5, pady = 5)

	e7 = tkinter.Entry(frame, textvariable = "", bg = 'khaki', fg = 'blue')
	e7.grid(row = 0, column = 0, columnspan = 7, sticky = 'sewn', pady = 5, padx = 5)
	button = tkinter.Button(frame, text = 'BROWSE', fg = 'khaki', bg = 'gray', font = ('Helvetica', 7), \
	command = lambda: var.set(filedialog.asksaveasfilename(initialdir = '/', title = 'Select Folder and Define MMSLAB file name',\
	filetypes = (('txt files','*.txt'), ('all files', '*.*')))))
	button.grid(row = 0, column = 7, columnspan = 1, sticky = 'we', padx = 5, pady = 5)


	regres_tabcont.add(regres_tab2, text = 'Polynomial Fitting')
	regres_tabcont.add(regres_tab3, text = 'Statistical Analysis')
	regres_tabcont.add(regres_tab4, text = 'Scatterplots')
	regres_tabcont.grid(row = 0, column = 0)
	regres_transient_toplevel.transient(root)
'''menuBar_5 = tkinter.Menubutton(iframe, text = 'Quick Corellation & Regression Analysis', font = ('Times New Roman', 9), relief = 'raised', borderwidth=0, width=35)
menuBar_5.pack(side = 'left')
corMenu = tkinter.Menu(menuBar_5, tearoff = 0)
menuBar_5['menu'] = corMenu'''

'''corMenu.add('command', label = 'Vector Data', compound = 'left', underline = 0, font = ('Helvetica', 9), command = regres)
corMenu.add('command', label = 'Raster Data', compound = 'left', underline = 0, font = ('Helvetica', 9), command = "")'''
def main():
	root.mainloop()
if __name__ == "__main__":
	main()
