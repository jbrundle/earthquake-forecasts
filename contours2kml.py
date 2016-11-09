'''
# Pyplot Contours --> to --> KML
# Mark R. Yoder, PhD
# Dept. of Physics, UC Davis
######
#
# this podule provides a simple interface to convert pyplot contours to kml.
# the basic method may be flawed a bit, but it seems to work pretty well... and it could use some serious rewriting and optimization
# (product of organic development and mixed requirements at the time of development).
# the initial approach here is to use a logical argument to determine how pyplot is drawing innter/outer polygons, pinching off
# polygins, etc. it might be the case that this information is contained more explicitly in the vertex objects (they're not (i don't think)
# just lists/tuples; they are an extended list-like object that contins information about how to get to the next/previous point.
# SO, thi seems to work pretty well most of the time, but it can make mistakes given complex polygons.
#
# notes and declarations:
# this code is for use on an 'as is' basis; it is not meant for mission critical applications. the contours -> kml conversion seems to be
# pretty reliable, but frankly most of it probably is not great code, so clean-up and contributions are very much invited.
#
# see test_xyz_to_kml() for a tutorial and working example.
# the basic meant and potatoes is a call like:
# kml_from_contours(cset=contours, colorbarname='napa_colorbar.png', open_file=True, close_file=True, contour_labels=None, top=top, bottom=bottom, fname_out=fname_out, alpha_kml=alpha_kml)
'''
import numpy
import pylab as plt
import matplotlib as mpl
import math

#######################################
# tests and tutorials:
#
def test_xyz_to_kml(fname_in='napa_etas.xyz', fname_out='napa_etas.kml', n_contours=15, fignum=0, bottom=0., top=1., alpha_kml=.8):
	# a test script. read in some xyz, make contours, output kml.
	# the front part of this script is data gathering; if you have pyplot contours and want kml, skip to the last two lines.
	# the moral of the script is to generate kml from contours like (you'll need to fill in some variable names, etc.):
	# kml_str = kml_from_contours(cset=contours, colorbarname='napa_colorbar.png', open_file=True, close_file=True, contour_labels=None, top=top, bottom=bottom, fname_out=fname_out, alpha_kml=alpha_kml)
	# a standard default execution will be:
	# kml_str = kml_from_contours(cset=contours_i_just_made, colorbarname='my_colorbar.png', open_file=True, close_file=True, contour_labels=None, top=1.0, bottom=0., fname_out='my_contours.kml', alpha_kml=.5)
	#
	# 1) read in the data file. we'll make contours from these data.
	# for now, assume file is in a format like: [ [x,y,z], ...]
	ary_xyz = []
	with open(fname_in, 'r') as f:
		for rw in f:
			if rw.startswith('#!'):
				# this gives us parameters/meta data.
				#prams = [rws.split('=') for rws in rw[2:].split()]
				prams = {x[0]:(None if str(x[1].lower()).startswith('none') else float(x[1])) for x in [rws.split('=') for rws in rw[2:].split()]}
				#return prams
				#prams = {x.split('=')[0]:(float(x.split('=')[1]) if x[1]!=None else None) for
			if rw[0]=='#': continue
			#
			ary_xyz += [[float(x) for x in rw.split()]]
		#
	#
	# ... and in this case, the data are properly gridded, so we could take any number of shortcuts to properly shape the array 
	# (aka, set(), max()/min(), etc. if there are issues, we'll just put it together old-school (aka, spin the list and
	# move to a new row/col when we encounter new coordinate values.
	#
	# first, sort (it's probably already sorted):
	ary_xyz.sort(key=lambda x: (x[1], x[0]))	# so sorted by y,x; reading sequentially will be rows of x and columns of y.
	Xs, Ys, Zs = (numpy.array(col) for col in  numpy.array(zip(*ary_xyz)))
	#
	# side lengths (number of unique entries in Xs, Ys).
	len_X = len(set(Xs))
	len_Y = len(set(Ys))
	#
	#
	for x in [Xs, Ys, Zs]: x.shape=(len_Y, len_X)
	#
	# ok, preparation part finished. this is where most people will start; get some contours using pyplot.contour() or .contourf() and convert to kml.
	#
	plt.figure(fignum)
	plt.ion()
	plt.clf()
	contours = plt.contourf(Xs, Ys, Zs, n_contours)
	#
	# now, write kml. there are two main function calls that can be made. 1 fetches a kml string, the other writes the string to file.
	# they call one another as necessary, so only one call is required.
	# write_kml_file(kml_str=None, cset=None, fout='conts.kml', colorbarname='scale.png', markers=None, top=1.0, bottom=0.0)
	# def kml_from_contours(cset=None, colorbarname='scale.png', open_file=True, close_file=True, contour_labels=None, top=1.0, bottom=0.0, fname_out=None, markers=None)
	#
	# some comments on parameters:
	# open_file, close_file : open/close the KML file. if we're appending to a file in progress, set "open" to False, and the script won't write the open-header.
	#   if we're going to write more stuff (aka, add some place markers), set "close" to False. the default is to create a self-contained kml file, 
	#   so open_file=True, close_file=True
	# top, bottom: top/bottom contours to keep. use fractions, like .9, .1 to keep the 10th to 90th percent contours. i think you can also give integers if you know
	#  the contour count, and the script will recognize integers vs floats (if not, it's an easy community contribution).
	# to just return a string, set fname_out=None, otherwise specify a filename for the kml output. 
	# alpha_kml: an alpha setting for the kml output (aka, fill color density).
	kml_str = kml_from_contours(cset=contours, colorbarname='napa_colorbar.png', open_file=True, close_file=True, contour_labels=None, top=top, bottom=bottom, fname_out=fname_out, alpha_kml=alpha_kml)
	
	return contours
	
	
####################################################################
####################################################################
# actual working code:
###################################

def arrayFromSet(cs=None, alpha='9d'):
	# get an array (list) from a contourset (matplotlib.contour.ContourSet object returned by contourf() call).
	# ... and invariably, there is a more efficient way to do this, so OpenSourcers are welcome to have at this...
	# format:
	# [ [z0, [[x,y], [x,y],...], [z1, [paths1]], 
	#
	#
	if alpha==None: alpha='9d'
	if type(alpha)==type(1) or type(alpha)==type(1.0): alpha='%02x' % alpha	# if alpha is a number, convert to string-hex
	#																									# (aka, '%02x' % 220 --> 'dc')
	#
	levels=cs.levels	# lower bound of contours
	layers=cs.layers	# upper bound of contours (at least for contours <0)
	dz= layers[0] - levels[0]
	collects=cs.collections
	carray = []
	for i in xrange(0, len(collects)):
		bgra_array = 255*collects[i].get_facecolor()[0]
		#strclr = '7d%02x%02x%02x' % ( bgra_array[2] , bgra_array[1] , bgra_array[0] )
		strclr = '%s%02x%02x%02x' % (alpha, bgra_array[2] , bgra_array[1] , bgra_array[0] )
		#
		for trace in collects[i].get_paths():
			#carray+=[[levels[i], layers[i], '%s' % strclr, []]]
			# objContour(low=None, high=None, alpha=None, R=None, G=None, B=None, verts=None, RGB=None)
			carray += [objContour(low=cs.levels[i], high=cs.layers[i], RGB=strclr, verts=[])]
		
			for lng, lat in trace.vertices:
				#carray[-1][-1]+=[[lng, lat]]
				carray[-1].verts+=[[lng, lat]]
				#
			#
		#
	#
	#
	return carray
#
def fixed_conts_to_file(cset=None, contours_string=None, file_out='contour_shapes.txt', top=1.0, bottom=0.0):
	if contours_string==None: contours_string = fixed_conts_to_str(cset=cset, top=top, bottom=bottom)
	#
	fout = open(file_out, 'w')
	fout.write(contours_string)
	fout.close()
	#
	return contours_string
#	
def fixed_conts_to_str(cset=None, top=1.0, bottom=0.0):
	#
	levels = cset.levels
	fixed_conts = conts_to_plot_lists(cset=cset)
	#				# contours are then like [ level0, level1, ...]
	#				# level: [poly0, poly1, poly2...]
	#				# poly: [X, Y]
	#
	n_levels = len(fixed_conts)-1	# and, of course, not exactly number of levels, but the max index.
	#top      = max(n_levels - int(math.ceil(top*float(n_levels))), 0)
	#bottom   = min(int(math.ceil(bottom*float(n_levels))), n_levels)	# (and just in case we accidentally calc <0 or >len)
	#top      = int(math.ceil(top*float(n_levels)))
	#bottom   = max(n_levels-int(math.ceil(bottom*float(n_levels))), 0)	# (and just in case we accidentally calc <0 or >len)
	top      = int(math.ceil(top*float(n_levels)))
	bottom   = int(math.floor(bottom*float(n_levels)))	# (and just in case we accidentally calc <0 or >len)
	#
	#
	outstr='# contour coordinates\n'
	outstr+='# this functionality is adapted somewhat hastily from some diagnostic functions, so it is a bit rough.\n'
	outstr+='# #!contour indicates a contour level (change). #!trace indicates a trace, or contour shape (borrowed \n'
	outstr+='# from the matplotlib vernacular). matplotlib produces contour shapes in a funny way, so there is a process \n'
	outstr+='# to fix this. in the process, some of the shapes are grouped into sets where there were inner/outer sets and \n'
	outstr+='# spatially independent shapes. this accounts for the notation "trace n,m", in which m indicates shapes within \n'
	outstr+='# a group of shapes n. all shapes within #!contour blocks belong to the same contour; aka, the n,m indices could \n'
	outstr+='# be flattened (and probably will in later iterations of this code.\n\n'
	#
	# levels:
	outstr+='#!contour levels:\t'
	for x in levels: outstr+='%f\t' % x
	outstr=outstr[:-1]+'\n'
	#
	#for i in xrange(len(cs)):
	for i in xrange(bottom, top):
		#
		outstr+='#!contour\t%d\n' % i
		#itrace=0
		for itraces, traces in enumerate(fixed_conts[i]):
		#for trace in cs[i]:
			for itrace, trace in enumerate(traces):
				#print trace
				#print len(trace), len(trace[0])
				#continue
				#tmpLL=[]		# KML is not always rendering properly. maybe incomplete polygons?
				outstr+='#!trace\t%d,%d\n' % (itraces, itrace)
				#
				for i_ll in xrange(len(trace[0])):
				#for lng,lat in trace.vertices:
					lng = trace[0][i_ll]
					lat = trace[1][i_ll]
				#for lng, lat in zip(*trace):
					#kmlstr+='%s,%s,0\n' % (round(lng,7),round(lat,7))
					outstr+='%f\t%f,0\n' % (lng, lat)    # and z=0
					#tmpLL+=[[lng, lat]]
				#if tmpLL[0]!=tmpLL[-1]:
				#	print "completing a polygon, %d." % fixedpolys
				#	fixedpolys+=1
				#	kmlstr+='%f,%f,0\n' % (tmpLL[0][0], tmpLL[0][1])
	#
	return outstr
#
def conts_to_str(cset=None, top=1.0, bottom=0.0):
	'''
	# top/bottom: max/min contours to be used. "top" is the highest level.
	'''
	# conts_to_plot_lists_raw(cset=None)
	
	cs=cset.collections
	levels = cset.levels
	#
	# what's our range?
	if top==None: top=1.0
	if bottom==None: bottom=0.0
	#
	# if top, bottom are >1, assume they're the number/index of the contour range being selected.
	# otherwise, it's a percentage.
	# first contour is the highest value.
	#
	n_levels = len(levels)-1
	#top      = max(n_levels - int(math.ceil(top*float(n_levels))), 0)
	#bottom   = min(int(math.ceil(bottom*float(n_levels))), n_levels)	# (and just in case we accidentally calc <0 or >len)
	#top      = int(math.ceil(top*float(n_levels)))
	#bottom   = max(n_levels-int(math.ceil(bottom*float(n_levels))), 0)	# (and just in case we accidentally calc <0 or >len)
	top      = int(math.ceil(top*float(n_levels)))
	bottom   = int(math.floor(bottom*float(n_levels)))	# (and just in case we accidentally calc <0 or >len)

	#
	outstr='# contour coordinates\n'
	#
	# levels:
	outstr+='#!contour levels:\t'
	for x in levels: outstr+='%f\t' % x
	outstr=outstr[:-1]+'\n'
	#
	#for i in xrange(len(cs)):
	
	for i in xrange(bottom, top):
		outstr+='#contour\t%d\n' % i
		itrace=0
		for itrace, trace in enumerate(cs[i].get_paths()):
		#for trace in cs[i]:
			#tmpLL=[]		# KML is not always rendering properly. maybe incomplete polygons?
			outstr+='#trace\t%d\n' % itrace
			for lng,lat in trace.vertices:
			#for lng, lat in zip(*trace):
				#kmlstr+='%s,%s,0\n' % (round(lng,7),round(lat,7))
				outstr+='%f\t%f,0\n' % (lng, lat)    # and z=0
				#tmpLL+=[[lng, lat]]
			#if tmpLL[0]!=tmpLL[-1]:
			#	print "completing a polygon, %d." % fixedpolys
			#	fixedpolys+=1
			#	kmlstr+='%f,%f,0\n' % (tmpLL[0][0], tmpLL[0][1])
	#
	return outstr


def innerouterpolys(polylist):
	#print "do nothing yet."
	# but, eventually: separate polygons that are "inside" two two polys (aka, and "outside" poly), and
	# associate their inner polys with them.
	# use shapely.geometry (sgp) module?? maybe not, as per compatiblity.
	#
	# note, pollys are coming like [ [[Xs], [Ys]], [] ], 
	# so the jth vertex of the ith poly is (x,y) = (polylist[i][0][j], polylist[i][1][j])
	polylistplus=[]		# indexed entries: [polyindex, [list of inners], [verts] ]
	#outerpolys=[]			# these will be lists. each entry is like: [[outer],[true-inner],[true-inner],..]
	#
	#for ply in polylist:
	for i in xrange(len(polylist)):
		polylistplus += [[i, [], polylist[i]]]
		#
		# shorthand:
		# which poly each polygon it is inside, and how many polys it is inside (aka, len(that list)).
		#x0,y0=poly[i][0][0], poly[i][1][0]	# we only need to test one point since we're starting with contours (don't cross).
		#x0,y0=polylistplus[-1][2][0][0], polylistplus[-1][2][1][0]	# we only need to test one point since we're starting with contours (don't cross).
		x0, y0 = polylist[i][0][0], polylist[i][1][0]
		#print "x0, y0: ", x0, y0
		# in general, we'd need to test all points to be inside.
		#
		# for each polygon in this level:
		# is the ith ("top") polygon inside the j'th poly?
		for j in xrange(len(polylist)):
			if j==i: continue 
			X,Y = polylist[j][0][:], polylist[j][1][:]
			#if x0>=max(X) or x0<=min(X) or y0>max(Y) or y0<min(Y): 
			#	print "outside max/min..."
			#	continue
			#
			if X[0]!=X[-1]:
				X+=[X[0]]
				Y+=[Y[0]]	# complete the poly...
			#
			N=len(X)
			ncrossings = 0
			# how many poly boundaries do we cross if we draw a line out of the poly in one direction.
			# equivalently (and in computer language), how many segments at y1 < y <y2 (or upside down)
		# are to the right of the point (or to the left, or up/down -- pick one).
			for k in xrange(1,N):
				k1 = k-1
				#k2 = (k+1)%N	# note the k%N: the poly does not have to be closed
				k2 = k	# but it should be, or you can count a crossing twice and get a bogus answer.
				x1, y1 = X[k1], Y[k1]
				x2, y2 = X[k2], Y[k2]
				'''
				x1,y1 = polylist[j][0][k1], polylist[j][1][k1]
				print "j,len(polylist):", j, len(polylist)
				print "k2, len(polylist[j][0])", k2, len(polylist[j][0])
				print "k2, len(polylist[j][1])", k2, len(polylist[j][1])
				x2,y2 = polylist[j][0][k2], polylist[j][1][k2]
				'''
				#if y0>=min(y1, y2) and y0<=max(y1, y2) and x0<max(x1, x2):
				#
				if x0>=min(x1, x2) and x0<max(x1, x2):	# note, one must be <= and the other < or, if we're on a grid -- and
					fx = (y2-y1)/(x2-x1)*(x0-x1)			# we're always on a grid, we'll count each crossing (left seg then right).
					if fx>(y0-y1):							# that's why the test script worked (x1,y1 != x2, y2) and the production
						ncrossings += 1						# failed...
			# clean up a bit:
			X=None
			Y=None
			#
			#print i,j,j,ncrossings, ncrossings%2
			if ncrossings%2==1:
				# i'th poly is inside j'th poly...
				polylistplus[-1][1] += [j]	
			#
		#
		# so now, we have a list of polygons and the polys they are inside.
		# for outerPolys, len(innersList)%2==0. if polyA is inside polyB and nA-nB=1, polyA is an inner-poly to polyB
		#
		#for rw in polylistplus:
		#	print rw
		outerpolys=[]
		for ply1 in polylistplus:
			#print "inner-len: ", len(ply1[1]), len(ply1[1])%2
			if len(ply1[1])%2==0:
				#print "***", len(ply1[1])
				# it's an outer poly...
				outerpolys+=[[ply1[2]]]
				for ply2 in polylistplus:
					# find its inners:
					if ply2==ply1: continue
					if len(ply2[1])%2==0: continue	# skip outers...
					#
					#print len(ply2[1]), (len(ply1[1])+1)
					if ply1[0] in ply2[1] and len(ply2[1])==(len(ply1[1])+1):
						# the outer poly's index is in ply2's "inside-list", then ply2 is inside ply...
						# AND, ply2 is one "deeper" (inside exactly one more poly) thatn ply
						outerpolys[-1]+=[ply2[2]]
	#
	#return polylist
	return outerpolys
#
def conts_to_plot_lists(cset=None):
	# this corrects a mpl/kml mismatch. mpl likes to plot inner polygons by drawing a line from the outside
	# ring to the inner poly, draw the inenr poly, then back to the outer ring over the first line. mpl
	# interprets the line as infinitely thin and basically ignores it. KML will not ignore that line, and so
	# you get a big fat mess. This script separates the outer/inner groups into outer and inner poly rings.
	#
	# BUT, (2012 09 05) this needs to be revisited. inner-inner polys are not plotting correctly. aka,
	# the z-profile like: /-\_/-\_/-\ (like an s0 within an s1) need to be considered. in such a case,
	# we get: <outer><inner><inner></inner></inner></outer>, but there should be 1 outer/inner poly then a
	# separate poly entirely. the second "inner" needs to be pulled out as a separate poly. (see similar note below)
	#
	# so, polys that are inside an odd number of polys are "inner" polys; polys inside an even number are "outside" 
	# any poly (as per knot theory), or otherwise constitute an "outer" poly. poly1 is an "inner" poly of poly2 if:
	#   1: it is inside poly2
	#   2: its "inner index n2 = n1 - 1
	#
	cs=cset.collections
	outlist=[]	# the length of this will be the number of levels: (or contours?)
					# outlist -> [ [level0: [[cont0x], [cont0y]], [[cont1x], [cont1y]] ], [level1: ] ???
	#levels=[]
	#contlevel=0
	#
	for i in xrange(len(cs)):
		# level-level:
		outlist+=[[]]
		#contcount=0
		# each level will have multiple polygons:
		for trace in cs[i].get_paths():
			# each "trace" will be a polygon (there might be multiple polygons per level), which might have "internal" polygons.
			tmpLL=[]		# KML is not always rendering properly. maybe incomplete polygons?
			#outlist[-1]+=[[[],[]]]
			newpoly=[ [[], []] ]	# newpoly[poly-index][x=0 or y=1][x_i or y_i]. note, outer poly is first element.
			# draw the polygons (first the outside one, then the inside ones -- presumably)
			llindex=0
			for lng,lat in trace.vertices:
				#
				# fix mangled contours:
				#
				newpoly[-1][0]+=[lng]
				newpoly[-1][1]+=[lat]
				#
				# have we closed the polygon (in the middle?):
				# what is happening is that the inner polys are being drawn screwy ways.
				# the errors are not a total mess; basically, the link between then inner and outer poly is drawn
				# from the last/first point rather than clostest two points.
				# the first poly is probably the outer-most poly. subsequent closed sequences are inside polys.
				# however, see comments above. the inner-inner->outer polys must be separated. do the following procedure
				# first, then pass that list of polys for further parsing.
				#
				# SO, make a list of all closed objects in a single poly object, then insert them into the outer sequence
				# (noting the inner poly has to close and then return to the starting point on the outer ring.
				# is newpoly[-1] == newpoly[0] ?
				if newpoly[-1][0][-1] == newpoly[-1][0][0] and newpoly[-1][1][-1] == newpoly[-1][1][0] and len(newpoly[-1][0])>1:
					# the polygon has closed.
					# but has it really, or is this just a "bay" or a "nub" for which the closing boundaries intersect.
					# if this poly is ended and a new poly folows, the next point is a line-segment connecting the next
					# (inner) poly. if it is a bay or a nub, the shape will close on itself again. I THINK, that we can say
					# this poly is truly closed IF this end coordinate does not appear again in the sequence.
					# there is probably a smart, compiled way to do this, but for now let's just loop it:
					ispinched=True
					for i in xrange(llindex+1, len(trace.vertices)):
						if trace.vertices[i][0]==lng and trace.vertices[i][1]==lat:
							# we've found our "end" point at least one mor time, so keep building the poly.
							# obviously, there is a faster way to do this, but let's just see if it works first.
							ispinched=False
							break
					#
					if ispinched==True:
						# pinch it off and start a new poly.
						newpoly+=[ [[],[]] ]
						tmpLL=[]
						#contcount+=1
				#
				llindex+=1
			#
			# now, clean up a bit:
			# if it's a bogus poly -- probably a line back to a prev. poly (each poly is like: poly0, line (segment) to poly1,
			#		poly1, either: {line back to poly0, line to poly2, nothing?}
			# i think the problem we're having is multi-inner polys aka, profile like 1 - 0 - 1. (ring=1, dip-ring=0, bump in middle=1)
			# ___|-\_/-\_/-|___ . the inner-inner bit is causing a problem. to fix this, i think, we have to draw all the polys
			# and explicitly determine which polys they are inside. (<outer><inner><inner></innter></inner></outer> ??)
			# note: a poly-ring inside an even  number of polys is an outer ring; a poly inside an odd number of rings is an inner ring.
			if len(newpoly[-1][0])<2: newpoly.pop()
			if newpoly[-1][0][-1] != newpoly[-1][0][0] and newpoly[-1][1][-1] != newpoly[-1][1][0]:
				newpoly.pop()
			#
			# at this point, newpoly is like: [[likely-outermost],[maybe-inner], [maybe-inner], ... ]
			# where we infer, or have in the past, that the first row represents the outer-most poly; subsequent polys are inners.
			# this is not the cast. find the 2xinner (aka, outer) polys within as per comments above.
			#
			# break up newpoly into inner-outer groups
			newpolylist = innerouterpolys(newpoly)
			#if len(newpoly)>2:
			#	print "len(newpoly), len(newpolylist): %d, %d" % (len(newpoly), len(newpolylist))
				#print newpoly
				#print newpolylist
			#newpolylist = [newpoly]
			# return the full list and interpret i=0 as the outer ring; i>0 as the inner boundaries.
			#
			#outlist[-1]+=[newpoly]	# add this trace(set) to the last level of outlist.
			for ply in newpolylist:
				#outlist[-1]+=[newpoly]
				outlist[-1]+=[ply]
			#
			#contcount+=1
		#contlevel+=1
	return outlist

#
def conts_to_plot_lists_raw(cset=None):
	'''
	# cset: return value from contour() or contourf()
	#this is probably pretty raw and will show the problems and mismatch between pyplot and kml contours. see conts_to_plot_lists() (not raw)
	'''
	cs=cset.collections
	outlist=[]	# the length of this will be the number of levels: (or contours?)
					# outlist -> [ [level0: [[cont0x], [cont0y]], [[cont1x], [cont1y]] ], [level1: ] ???
	#levels=[]
	#contlevel=0
	#
	for i in xrange(len(cs)):
		# level-level:
		outlist+=[[]]
		#contcount=0
		# each level will have multiple polygons:
		for trace in cs[i].get_paths():
			# each "trace" will be a polygon (there might be multiple polygons per level), which might have "internal" polygons.
			tmpLL=[]		# KML is not always rendering properly. maybe incomplete polygons?
			#outlist[-1]+=[[[],[]]]
			newpoly=[ [[], []] ]	# newpoly[poly-index][x=0 or y=1][x_i or y_i]. note, outer poly is first element.
			# draw the polygons (first the outside one, then the inside ones -- presumably)
			for lng,lat in trace.vertices:
				#
				# fix mangled contours:
				#
				newpoly[-1][0]+=[lng]
				newpoly[-1][1]+=[lat]
				#
			#
			outlist[-1]+=[newpoly]	# add this trace(set) to the last level of outlist.
			#
			#contcount+=1
		#contlevel+=1
	return outlist	
#
def contsToPlotListsOuter(cset=None):
	# simple polygons, not distinguishing inner/outer.
	# @cset : return from contour() or contourf()
	#if cset==None: cset=self.conts
	#
	#
	cs=cset.collections
	outlist=[]	# the length of this will be the number of levels: (or contours?)
					# outlist -> [ [level0: [[cont0x], [cont0y]], [[cont1x], [cont1y]] ], [level1: ] ???
	contlevel=0
	#
	for i in xrange(len(cs)):
		# level-level:
		outlist+=[[]]
		contcount=0
		for trace in cs[i].get_paths():
			tmpLL=[]		# KML is not always rendering properly. maybe incomplete polygons?
			outlist[-1]+=[[[],[]]]
			for lng,lat in trace.vertices:
				#
				# fix mangled contours:
				#
				outlist[-1][-1][0]+=[lng]
				outlist[-1][-1][1]+=[lat]
				tmpLL+=[[lng, lat]]
				#
				# have we closed the polygon (in the middle?):
				# what is happening is that the inner polys are being drawn screwy ways.
				# the errors are not a total mess; basically, the link between then inner and outer poly is drawn
				# from the last/first point rather than clostest two points.
				# the first poly (i think) is the outer poly. subsequent closed sequences are inner polys (excluded regions).
				# SO, make a list of all closed objects in a single poly object, then insert them into the outer sequence
				# (noting the inner poly has to close and then return to the starting point on the outer ring.
				if outlist[-1][-1][0][-1] == outlist[-1][-1][0][0] and outlist[-1][-1][1][-1] == outlist[-1][-1][1][0] and len(outlist[-1][-1][0])>1:
					# the polygon has closed. pinch it off and start a new poly.
					outlist[-1]+=[[[],[]]]
					tmpLL=[]
					contcount+=1

			if len(outlist[-1][-1][0])<2: outlist[-1].pop()
			if outlist[-1][-1][0][-1] != outlist[-1][-1][0][0] and outlist[-1][-1][1][-1] != outlist[-1][-1][1][0]:
				outlist[-1].pop()
			# (outlist[contour level][polynum][x or y][x,y index] )
			
			
			'''
			#if tmpLL[0]!=tmpLL[-1]:
			while tmpLL[0]!=tmpLL[-1]:
			#	print "completing a polygon, %d." % fixedpolys
				#fixedpolys+=1
				#outlist[-1][-1][0] += [tmpLL[0][0]]
				#outlist[-1][-1][1] += [tmpLL[0][1]]
				outlist[-1][-1][0].pop()
				outlist[-1][-1][1].pop()
				tmpLL.pop()
			'''
				#
			contcount+=1
		contlevel+=1
	return outlist

def plotPolyList(polys=None, markerstr='.-', fignum=0, pllevels=None):
	#
	if type(pllevels).__name__ in ('float', 'int'): pllevels=[pllevels]
	plt.figure(fignum)
	plt.clf()
	plt.ion()
	nlevels=len(polys)
	
	#
	ilevel=0
	for level in polys:
		#
		if ilevel not in pllevels:
			ilevel+=1
			continue
			
			 
		lvlclr=plotcolor(ilevel, nlevels*2)
		print ilevel, lvlclr
		for xy in level:
			for ring in xy:
			#fixedpolys = checkPoly(inpoly=xy, fignum=None)
			#for fxy in fixedpolys:
				plt.plot(ring[0], ring[1], markerstr, color=lvlclr)
			#plt.plot(xy[0], xy[1], markerstr, color=lvlclr)
		ilevel+=1
	return polys			

def plotPolyListOuter(polys=None, markerstr='.-', fignum=0, pllevels=None):
	''''
	# this was, i think, primarily a development and diagnostic script, basically to separate and plot outer polygons.
`	# note the @polys input is like the return value from contsToPlotListsOuter()
	#
	'''
	#		
	# plot the simpler "just outers" version of the polygon list.
	#
	if type(pllevels).__name__ in ('float', 'int'): pllevels=[pllevels]
	plt.figure(fignum)
	plt.clf()
	plt.ion()
	nlevels=len(polys)
	
	#
	ilevel=0
	for level in polys:
		#
		if ilevel not in pllevels:
			ilevel+=1
			continue
			
			 
		lvlclr=plotcolor(ilevel, nlevels*2)
		print ilevel, lvlclr
		for xy in level:
			#fixedpolys = checkPoly(inpoly=xy, fignum=None)
			#for fxy in fixedpolys:
			#	plt.plot(fxy[0], fxy[1], markerstr, color=lvlclr)
			plt.plot(xy[0], xy[1], markerstr, color=lvlclr)
		ilevel+=1
	return polys

def get_marker_kml_str(markers):
	#kmlstr='<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2">\n<Document>\n'
	kmlstr=''
	for ev in markers:
		thisdt=str(mpd.num2date(ev[0]))
		datestr=thisdt.split('.')[0]
		kmlstr+='<Placemark><name>%s, %.2f</name><Point><coordinates>%f,%f,0</coordinates></Point></Placemark>' % (datestr, ev[3], ev[2], ev[1])
	
	return kmlstr
	
def write_kml_file(kml_str=None, cset=None, fout='conts.kml', colorbarname='scale.png', markers=None, top=1.0, bottom=0.0):
	#
	# depricated: wrapping this into kml_from_contours()
	#
	if kml_str==None: kml_str=kml_from_contours(cset=cset, colorbarname=colorbarname, open_file=True, close_file=False, top=top, bottom=bottom, fname_out=None)
	#
	#
	# this bit can be used to write marker objects into the kml file (nominally earthquakes or other point-like locations).
	if markers !=None:
		# we've been given some earthquakes to plot as well.
		kmlstr+='\n'
		kmlstr+=get_marker_kml_str(markers)
	kmlstr+='</Document>\n'
	kmlstr+='</kml>'
	#
	with open(fout, 'w') as f:
		f.write(kmlstr)
	#
	return None

def kml_from_contours(cset=None, colorbarname='scale.png', open_file=True, close_file=True, contour_labels=None, top=1.0, bottom=0.0, fname_out=None, fname_out_mode='w', markers=None, alpha_kml=.8):
	# get a KML string from a contourset (matplotlib.contour.ContourSet object returned by contourf() call).
	# (will it be necessary to write directly to file? is there a string length limit problem?)
	# open_file: does this string "open the file", aka, open the kml file (i think that's what this means). if True, we get
	# a stand-alond kml; if False, we get a fragment that we can stick into an existing file.
	# close_file: same idea: is this the end of the file? if true, add a footer; if not, it's a fragment.
	# contour_labels: a list [] of contour labels like ['freezing', 'cold', 'cool', 'luke-warm', 'warm', 'hot', 'scalding'].
	#                 we should be able to parse labels to contours more or less directly.
	# top/bottom: top and bottome ranges (percentages/fractions) to output. aka, if (top,bottom) = (.5, 1.), export the 
	# top 50% of contours to kml.
	# fname_out: filename for output (will call write_kml_file) note this addition depricates write_kml_file()
	#
	cs=cset.collections
	#
	polys = conts_to_plot_lists(cset=cset)		# this function fixes multi-closed polys. use for actual kml contours.
															# still has same polys[cont. level][contour index][x=0,y=1][x,y index] format
															# the polys array will differ from cset.collections in the number of polys
															# per level; the number of levels should be the same.
	#
	# contour range to render:
	if top==None: top=1.0
	if bottom==None: bottom=0.0
	# convert alpha_kml (if it's a float) to hex:
	if isinstance(alpha_kml, float):
		# there's probably a better way to do this, but this seems to work ok:
		# hex(integer) returns something like '0xaa', like hex(128) --> '0x80', where the bit to the right of 'x' is the hex value.
		alpha_kml = hex(int(255*alpha_kml)).split('x')[1]
	#
	# if top, bottom are >1, assume they're the number/index of the contour range being selected.
	# otherwise, it's a percentage.
	# first contour is the highest value.
	#
	levels = cset.levels
	n_levels = len(levels)-1
	kml_top      = int(math.ceil(top*float(n_levels)))
	kml_bottom   = int(math.floor(bottom*float(n_levels)))	# (and just in case we accidentally calc <0 or >len)
	#
	print "kml bottom, top: ", kml_bottom, kml_top
	#
	#resolution = 5
	#LevelsNumber = 5 * resolution
	if contour_labels==None:
		#contour_labels = ['No','Low','Guarded','Elevated','High','Severe']
		contour_labels = ['very low', 'low', 'medium', 'high', 'quite high']
	#
	resolution = int(len(polys)/len(contour_labels))
	#startindex=resolution
	startindex=0
	kmlstr=''
	#
	if open_file==True:
		kmlstr='<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2">\n<Document>\n'
	#
	# styles come from the contours collection (cs):
	#for i in xrange(0,len(cs)):
	for i in xrange(kml_bottom, kml_top):
		bgra_array = 255*cs[i].get_facecolor()[0]
		kmlstr+='<Style id="l%d">\n' % i
		kmlstr+='<LineStyle><color>00000000</color></LineStyle>\n'
		kmlstr+='<PolyStyle><color>%s%02x%02x%02x</color></PolyStyle>\n' % (alpha_kml, bgra_array[2] , bgra_array[1] , bgra_array[0] )
		#kmlstr+='<PolyStyle><color>ff%02x%02x%02x</color></PolyStyle>\n' % ( bgra_array[2] , bgra_array[1] , bgra_array[0] )
		kmlstr+='</Style>\n'
	#
	#
	kmlstr+='<ScreenOverlay id="scale">\n'
	kmlstr+='<name>Color Scale</name>\n'
	kmlstr+='<Icon><href>scale.png</href></Icon>\n'
	kmlstr+='<Icon><href>%s</href></Icon>\n' % colorbarname
	kmlstr+='<overlayXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n'
	kmlstr+='<screenXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n'
	kmlstr+='<size x="0" y="0" xunits="pixels" yunits="pixels"/>\n'
	kmlstr+='</ScreenOverlay>\n'
	#
	#fixedpolys=0
	# now, stake out the contours from the polys array.
	# note: len(cs) = len(polys) (same number of levels)
	#	   len(cs[i]) != len(polys[i]) (in fact, len(polys[i]>=len(cs[i]) ), which is to say
	#	   polys may have more distinct polygons per contour level. of course, the individual polys
	#	   will have differend length as well.
	#for i in xrange(startindex, len(cs)):
	#for i in xrange(startindex, len(polys)):
	for i in xrange(kml_bottom, kml_top):
		# each i is a contour level.
		kmlstr+='<Placemark>\n'
		#
		warningindex=int(float(i)/float(resolution))
		if warningindex>(len(contour_labels)-1): warningindex=len(contour_labels)-1	# in the event that we have off-integer numbers.
		kmlstr+='<name>%s Risk</name>\n' % contour_labels[warningindex]
		kmlstr+='<styleUrl>#l%d</styleUrl>\n' % i
		kmlstr+='<MultiGeometry>\n'
		 
		for ii in xrange(len(polys[i])):
			# each ii is a polygon (set).
			kmlstr+='<Polygon>\n'
			kmlstr+='<extrude>0</extrude>\n'
			kmlstr+='<altitudeMode>clampToGround</altitudeMode>\n'
			kmlstr+='<outerBoundaryIs>\n'
			kmlstr+='<LinearRing>\n'
			kmlstr+='<coordinates>\n'
			
			#tmpLL=[]		# KML is not always rendering properly. maybe incomplete polygons?
			#for lng,lat in trace.vertices:
			for ill in xrange(len(polys[i][ii][0][0])):
				# first set is the outerBoundary
				# noting that each polygon is stored as [[x], [y]], so len(polys[i][ii])=2 always.
				# len(polys[i][ii][0])=len(polys[i][ii][1]) is the length or size (num. verts.) of the polygon.
				lng=polys[i][ii][0][0][ill]
				lat=polys[i][ii][0][1][ill]
				#
				kmlstr+='%f,%f,0\n' % (lng, lat)
				#tmpLL+=[[lng, lat]]
			#if tmpLL[0]!=tmpLL[-1]:
			#	print "completing a polygon, %d." % fixedpolys
			#	fixedpolys+=1
			#	kmlstr+='%f,%f,0\n' % (tmpLL[0][0], tmpLL[0][1])

			kmlstr+='</coordinates>\n</LinearRing>\n</outerBoundaryIs>\n'
			#
			# inner polys?
			for iInner in xrange(1,len(polys[i][ii])):
				thispoly=polys[i][ii][iInner]
				# if any, these will be inner polys, each like [[x], [y]]
				kmlstr+='<innerBoundaryIs>\n<LinearRing>\n<coordinates>\n'
				# coords like 'lon,lat,alt\n'
				# -77.05668055019126,38.87154239798456,100
				for ill in xrange(len(thispoly[0])):
					kmlstr+='%f,%f,0\n' % (thispoly[0][ill], thispoly[1][ill])
				kmlstr+='</coordinates>\n</LinearRing>\n</innerBoundaryIs>\n'
			#
			kmlstr+='</Polygon>\n'
		#
		kmlstr+='</MultiGeometry>\n'
		kmlstr+='</Placemark>\n'
	#
	if markers !=None:
		# we've been given some earthquakes to plot as well.
		kmlstr+='\n'
		kmlstr+=get_marker_kml_str(markers)
	#
	if close_file==True:
		kmlstr+='</Document>\n'
		kmlstr+='</kml>'
	#
	if fname_out!=None:
		#z=write_kml_file(kml_str=kmlstr, cset=cset, fout=fname_out, colorbarname=colorbarname, markers=markers, top=top, bottom=bottom)
		with open(fname_out, fname_out_mode) as f:
			f.write(kmlstr)
	#
	return kmlstr

#
def kml_from_contours_raw(cset=None, colorbarname='scale.png', contour_labels=['No','Low','Guarded','Elevated','High','Severe']):
	# get a KML string from a contourset (matplotlib.contour.ContourSet object returned by contourf() call).
	# (will it be necessary to write directly to file? is there a string length limit problem?)
	#
	# (use this script to visualize the mpl/kml compatibility problem).
	#
	cset_collections=cset.collections
	#
	startindex=0
	#
	#
	kmlstr='<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://www.opengis.net/kml/2.2">\n<Document>\n'
	#
	for i in xrange(0,len(cset_collections)):
		bgra_array = 255*cset_collections[i].get_facecolor()[0]
		kmlstr+='<Style id="l%d">\n' % i
		kmlstr+='<LineStyle><color>00000000</color></LineStyle>\n'
		kmlstr+='<PolyStyle><color>7d%02x%02x%02x</color></PolyStyle>\n' % ( bgra_array[2] , bgra_array[1] , bgra_array[0] )
		kmlstr+='</Style>\n'
	#
	#
	kmlstr+='<ScreenOverlay id="scale">\n'
	kmlstr+='<name>Color Scale</name>\n'
	#kmlstr+='<Icon><href>scale.png</href></Icon>\n'
	kmlstr+='<Icon><href>%s</href></Icon>\n' % colorbarname
	kmlstr+='<overlayXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n'
	kmlstr+='<screenXY x="0" y="1" xunits="fraction" yunits="fraction"/>\n'
	kmlstr+='<size x="0" y="0" xunits="pixels" yunits="pixels"/>\n'
	kmlstr+='</ScreenOverlay>\n'
	#
	fixedpolys=0
	for i in xrange(startindex, len(cset_collections)):
		kmlstr+='<Placemark>\n'
		kmlstr+='<name>%s</name>\n' % contour_labels[i/int(labels_len)]
		kmlstr+='<styleUrl>#l%d</styleUrl>\n' % i
		kmlstr+='<MultiGeometry>\n'
		 
		for trace in cset_collections[i].get_paths():
			kmlstr+='<Polygon>\n'
			kmlstr+='<extrude>0</extrude>\n'
			kmlstr+='<altitudeMode>clampToGround</altitudeMode>\n'
			kmlstr+='<outerBoundaryIs>\n'
			kmlstr+='<LinearRing>\n'
			kmlstr+='<coordinates>\n'
			
			tmpLL=[]		# KML is not always rendering properly. maybe incomplete polygons?
			for lng,lat in trace.vertices:
				#kmlstr+='%s,%s,0\n' % (round(lng,7),round(lat,7))
				
				# fixing contours mangled by pyplot:
				# an older attempted fixing routine that didn't quite work out.
				# keep this version as the "raw" form; see other KML function that 
				# uses a repaired set of polygons.
				#
				thislng = lng	# ... but i don't think that rounding error turned out to be the problem...
				thislat = lat
				#
				kmlstr+='%f,%f,0\n' % (thislng, thislat)
				tmpLL+=[[thislng, thislat]]
			#if tmpLL[0]!=tmpLL[-1]:
			#	print "completing a polygon, %d." % fixedpolys
			#	fixedpolys+=1
			#	kmlstr+='%f,%f,0\n' % (tmpLL[0][0], tmpLL[0][1])

			kmlstr+='</coordinates>\n'
			kmlstr+='</LinearRing>\n'
			kmlstr+='</outerBoundaryIs>\n'
			kmlstr+='</Polygon>\n'
		#
		kmlstr+='</MultiGeometry>\n'
		kmlstr+='</Placemark>\n'
	#
	kmlstr+='</Document>\n'
	kmlstr+='</kml>'
	#
	return kmlstr

def makeColorbar(z_vals=None, colorbarname=None, reffMag=5.0, fnum=None, fontcolor='k', warnings = ['No','Low','Guarded','Elevated','High','Severe'], color_map=None, label_str=''):
	#
	if colorbarname==None: colorbarname='scale.png'
	if color_map==None: color_map = mpl.cm.spectral
	#
	startindex=0
	#
	plt.figure(fnum)
	plt.close()
	#
	fig = plt.figure(num=fnum, figsize=(1,3))
	plt.clf()
	axe = fig.add_axes([0.05, 0.05, 0.2, 0.9])
	#
	fig.figurePatch.set_alpha(0)
	#
	matplotlib.rc('font', family='monospace', weight='black')
	#
	#ratefactorExp=mc-reffMag				# so the raw rate probably should be with mc_catalog = mc_etas.
	ratefactorExp=0.0											# but, we'll want to report some rate of m>reffMag
	norm = mpl.colors.Normalize(vmin=self.Z2d.min(), vmax=self.Z2d.max())
	timefactExp = math.log10(year2secs)

	#
	# let's add intermediate ticks: half and half-half:
	ticdiff0=self.Z2d.max()-self.Z2d.min()
	ticdiff=ticdiff0/4.0
	#
	print "max set: ", self.Z2d.max(), timefactExp, ratefactorExp
	#tics = [self.Z2d.min(), self.Z2d.max()]
	tics = [self.Z2d.min()]
	tcklbls=['%.2f'% (self.Z2d.min() + timefactExp + ratefactorExp)]
	while tics[-1]<=self.Z2d.max():
		tics+=[tics[-1]+ticdiff]
		tcklbls+=['%.2f'% (tics[-1]+ timefactExp + ratefactorExp)]
	#print t1, t2
	#
	#cb1 = mpl.colorbar.ColorbarBase(axe, norm=norm, ticks=tics, format="%g%%", orientation='vertical')
	cb1 = mpl.colorbar.ColorbarBase(axe, norm=norm, ticks=tics, format="%.2f", orientation='vertical', cmap=color_map)
	#cb1.set_ticklabels([t1, t2], update_ticks=True)
	cb1.set_ticklabels(tcklbls, update_ticks=True)
	#
	cb1.set_label('%s' % label_str, color=fontcolor)
	#cb1.set_label('ETAS rate')
	plt.savefig(colorbarname)	# vertical
	#
	#suffix=''
	i=len(colorbarname)-1
	#while colorbarname[i]!='.':
	#	#suffix.insert(0, colorbarcopy[-i])
	#	# this will get us the last '.'
	#	i-=1
	i=colorbarname.rfind('.')
	hname=colorbarname[:i] + '-h' + colorbarname[i:]
	#
	im = ipp.open(colorbarname, 'r')
	im.load()							# just to be sure. open() only creates a pointer; data are not loaded until analyzed.
	im=im.rotate(-90)
	#im.show()
	im.save(hname)
	#

	#
	return [colorbarname, hname]

