#!/usr/bin/env python3
# -*- coding: iso-8859-1 -*-
import  re
import  math



number = re.compile( "^[0-9\.\-eE]+$" )

stack = []
formt = "%.10lg\n"
memor = {}



def __print():
	global	stack
	print()
	for i in range( len( stack ) ):
		print( formt%( stack[i] ) )
	print()



def calc( text ):
	global	stack, formt, memor
	try:
#		print( stack )
		for itm in text.strip().split():
#			print( itm )
			if( number.match( itm ) and itm != "-" ):
				stack.append( float( itm ) )
			elif( itm[0] == ">" ):
				memor[itm[1:]] = stack.pop()
			elif( itm[0] == "<" and itm[1:] in memor ):
				stack.append( memor[itm[1:]] )
			elif( itm == "+" ):
				if( len( stack ) < 2 ):
					raise Exception( "ERROR: '+' wrong number of parameters" )
				op1 = stack.pop()
				op2 = stack.pop()
				stack.append( op2 + op1 )
			elif( itm == "-" ):
				if( len( stack ) < 2 ):
					raise Exception( "ERROR: '-' wrong number of parameters" )
				op1 = stack.pop()
				op2 = stack.pop()
				stack.append( op2 - op1 )
			elif( itm == "*" ):
				if( len( stack ) < 2 ):
					raise Exception( "ERROR: '*' wrong number of parameters" )
				op1 = stack.pop()
				op2 = stack.pop()
				stack.append( op2 * op1 )
			elif( itm == "/" ):
				if( len( stack ) < 2 ):
					raise Exception( "ERROR: '/' wrong number of parameters" )
				op1 = stack.pop()
				op2 = stack.pop()
				if( op1 == 0.0 ):
					stack.append( op2 )
					stack.append( op1 )
					raise Exception( "ERROR: '/' division by ZERO" )
				stack.append( op2 / op1 )
			elif( itm == "^" or itm == "**" ):
				if( len( stack ) < 2 ):
					raise Exception( "ERROR: '^' wrong number of parameters" )
				op1 = stack.pop()
				op2 = stack.pop()
				stack.append( math.pow( op2, op1 ) )
			elif( itm == "exp" ):
				if( len( stack ) < 1 ):
					raise Exception( "ERROR: 'exp' wrong number of parameters" )
				op1 = stack.pop()
				stack.append( math.exp( op1 ) )
			elif( itm == "ln" ):
				if( len( stack ) < 1 ):
					raise Exception( "ERROR: 'ln' wrong number of parameters" )
				op1 = stack.pop()
				stack.append( math.log( op1 ) )
			elif( itm == "alog" ):
				if( len( stack ) < 1 ):
					raise Exception( "ERROR: 'alog' wrong number of parameters" )
				op1 = stack.pop()
				stack.append( math.pow( 10.0, op1 ) )
			elif( itm == "log" ):
				if( len( stack ) < 1 ):
					raise Exception( "ERROR: 'log' wrong number of parameters" )
				op1 = stack.pop()
				stack.append( math.log10( op1 ) )
			elif( itm == "sin" ):
				if( len( stack ) < 1 ):
					raise Exception( "ERROR: 'sin' wrong number of parameters" )
				op1 = stack.pop()
				stack.append( math.sin( op1 ) )
			elif( itm == "cos" ):
				if( len( stack ) < 1 ):
					raise Exception( "ERROR: 'cos' wrong number of parameters" )
				op1 = stack.pop()
				stack.append( math.cos( op1 ) )
			elif( itm == "tan" ):
				if( len( stack ) < 1 ):
					raise Exception( "ERROR: 'tan' wrong number of parameters" )
				op1 = stack.pop()
				stack.append( math.tan( op1 ) )
			elif( itm == "asin" ):
				if( len( stack ) < 1 ):
					raise Exception( "ERROR: 'asin' wrong number of parameters" )
				op1 = stack.pop()
				stack.append( math.asin( op1 ) )
			elif( itm == "acos" ):
				if( len( stack ) < 1 ):
					raise Exception( "ERROR: 'acos' wrong number of parameters" )
				op1 = stack.pop()
				stack( math.acos( op1 ) )
			elif( itm == "atan" ):
				if( len( stack ) < 1 ):
					raise Exception( "ERROR: 'atan' wrong number of parameters" )
				op1 = stack.pop()
				stack.append( math.atan( op1 ) )
			elif( itm == "sqrt" ):
				if( len( stack ) < 1 ):
					raise Exception( "ERROR: 'atan' wrong number of parameters" )
				op1 = stack.pop()
				stack.append( math.sqrt( op1 ) )
			elif( itm == "neg" ):
				if( len( stack ) < 1 ):
					raise Exception( "ERROR: 'neg' wrong number of parameters" )
				op1 = stack.pop()
				stack.append( - op1 )
			elif( itm == "inv" ):
				if( len( stack ) < 1 ):
					raise Exception( "ERROR: 'inv' wrong number of parameters" ) 
				op1 = stack.pop()
				if( op1 == 0.0 ):
					stack.append( 0.0 )
					raise Exception( "ERROR: 'inv' division by ZERO" )
				stack.append( 1.0 / op1 )
			elif( itm == "deg" ):
				if( len( stack ) < 1 ):
					raise Exception( "ERROR: 'deg' wrong number of parameters" )
				op1 = stack.pop()
				stack.append( op1 * 180.0 / math.pi )
			elif( itm == "rad" ):
				if( len( stack ) < 1 ):
					raise Exception( "ERROR: 'rad' wrong number of parameters" )
				op1 = stack.pop()
				stack.append( op1 / 180.0 * math.pi )
			elif( itm == "pi" ):
				stack.append( math.pi )
			elif( itm == "clr" ):
				stack = []
			elif( itm == "dup" ):
				if( len( stack ) < 1 ):
					raise Exception( "ERROR: 'dup' wrong number of parameters" )
				stack.append( stack[-1] )
			elif( itm == "drop" ):
				if( len( stack ) < 1 ):
					raise Exception( "ERROR: 'drop' wrong number of parameters" )
				stack.pop()
			elif( itm == "sci" ):
				if( len( stack ) < 1 ):
					raise Exception( "ERROR: 'sci' wrong number of parameters" )
				op1 = stack.pop()
				formt = "%%.%dlg\n"%( op1 )
			elif( itm == "fix" ):
				if( len( stack ) < 1 ):
					raise Exception( "ERROR: 'fix' wrong number of parameters" )
				op1 = stack.pop()
				formt = "%%.%dlf\n"%( op1 )
			elif( itm == "eng" ):
				if( len( stack ) < 1 ):
					raise Exception( "ERROR: 'eng' wrong number of parameters" )
				op1 = stack.pop()
				formt = "%%.%dle\n"%( op1 )
			elif( itm == "swap" ):
				if( len( stack ) < 2 ):
					raise Exception( "ERROR: 'swap' wrong number of parameters" )
				op1 = stack.pop()
				op2 = stack.pop()
				stack.append( op1 )
				stack.append( op2 )
			elif( itm == "rot" ):
				if( len( stack ) < 3 ):
					raise Exception( "ERROR: 'rot' wrong number of parameters" )
				op1 = stack.pop()
				op2 = stack.pop()
				op3 = stack.pop()
				stack.append( op2 )
				stack.append( op1 )
				stack.append( op3 )
			elif( itm == "vars" ):
				for key in sorted( memor ):
					print( key, formt%( memor[key] ) )
			elif( itm == "_c" ):
				stack.append( 299792458.0 )
			elif( itm == "_na" ):
				stack.append( 6.0221415e23 )
			elif( itm == "_h" ):
				stack.append( 6.6260693e-34 )
			elif( itm == "_kb" ):
				stack.append( 1.3806505e-23 )
			elif( itm == "_r" ):
				stack.append( 8.314472 )
			elif( itm == "_qe" ):
				stack.append( 1.60217653e-19 )
			elif( itm == "_me" ):
				stack.append( 9.1093826e-31 )
			elif( itm == "_a0" ):
				stack.append( 0.5291772108 )
			elif( itm == "_ha" ):
				stack.append( 2625.49962955 )
	except Exception as error:
		print( error )
#	__print()





if( __name__ == "__main__" ):
	import	tkinter

	def cback( evnt ):
		global	stack, formt, stk, cmd
		buf = ""
		stk.delete( "@0,0", "end" )
		txt = cmd.get()
		if( txt == "vars" ):
			buf += "-- Variables --\n"
			for key in sorted( memor ):
				buf +=  "( " + key + " )\n" + formt%( memor[key] )
		elif( txt == "cons" ):
			buf = """-- Constants --
_c  : speed of light [m/s]
_na : Avogadro [1/mol]
_h  : Planck [J.s]
_kb : Boltzmann [J/K]
_r  : _kb * _na [J/(mol.K)]
_qe : e charge [C]
_me : e mass [kg]
_a0 : Bohr [A]
_ha : Hartree [kJ/mol]
"""
		else:
			calc( txt )
			for i in range( len( stack ) ):
				buf += formt%( stack[i] )
		stk.insert( "@0,0", buf )
		cmd.delete( 0, tkinter.END )

	app = tkinter.Tk()
	app.title( "RPN Calculator" )
	app.minsize( width = 320, height = 480 )
	app.maxsize( width = 320, height = 480 )
	app.resizable( tkinter.NO, tkinter.NO )
	frm = tkinter.Frame( app )
	frm.pack( side = tkinter.TOP, expand = 1, fill = tkinter.BOTH )
	cmd = tkinter.Entry( frm, font = ( "Courier New", "18" ), justify = tkinter.LEFT )
	cmd.pack( side = tkinter.BOTTOM, anchor = tkinter.SW, expand = 1, fill = tkinter.X )
	cmd.bind( '<Return>', cback )
	stk = tkinter.Text( frm, font = ( "Courier New", "18" ), bg = "#e6e6e6" )
	stk.pack( side = tkinter.TOP, expand = 1, fill = tkinter.BOTH )
	#app.wm_attributes( "-topmost", True )
	cmd.focus_force()
	#import os
	#os.system( "/usr/bin/osascript -e 'tell app \"Finder\" to set frontmost of process \"Python\" to true'" )
	app.mainloop()
