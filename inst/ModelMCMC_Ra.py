#include <Rcpp.h>
#include <Python.h>
	#####################  Librerias
from numpy import isnan,savetxt,genfromtxt, array, log, unique, exp, append,concatenate,zeros, repeat,linspace
from scipy.stats import uniform as unif
from numpy.random import seed



	#import csv
	##################### Read calibration curve

def plumMCMC(dirt,corename,T_mod,num_sup,det_lim,iterations , by,shape1_m,mean_m,shape_acc,mean_acc,fi_mean,fi_acc,As_mean,As_acc,resolution,seeds,thi,burnin,bqkg,Cs,Sdate,CSTdate):
	seed(int(seeds))
	plomo="/"+corename+".csv"
	fimean=fi_mean
	shapefi=fi_acc
	ASmaean=As_mean
	shapeAS=As_acc
	shape2_m= (shape1_m*(1-mean_m) )/mean_m
	scale_acc=mean_acc/shape_acc
	scale_fi=fimean/shapefi
	scale_As=ASmaean/shapeAS
	Data=genfromtxt (dirt+'Results '+corename+plomo, delimiter = ',')
	print(Data)
	##################### Data definition 210Pb
	if bqkg:
		Bqkg_cons=10.
	else:
		Bqkg_cons=500./3.#1000.
	density=Data[:,1]   * Bqkg_cons
	activity=Data[:,2]
	sd_act=Data[:,3]
	thic=Data[:,4]
	depth=Data[:,0]
	supp=Data[:,5]
	sd_supp=Data[:,6]
	Ran=len(supp)


	activity=activity*density
	sd_act=sd_act*density

	lam=0.03114

      


	dep_time_data=append(depth-thic,depth)
	dep_time_data=list(set(dep_time_data))
	X1, X0= [],[]
	for i1 in range(len(depth)):
		for k1 in range(len(dep_time_data)):
			if depth[i1]== dep_time_data[k1]:
				X1=append(X1,int(k1))
			if (depth-thic)[i1]== dep_time_data[k1]:
				X0=append(X0,int(k1))

	m=1
	breaks=array(m*by)
	while m*by< depth[-1]:
		m += 1
		breaks=append(breaks,m*by)


	#################### Functions

	def support(param):
		tmp3=True
		for i in param:
			if i <= 0.:
				tmp3=False
		if param[1+Ran]>=1.:
			tmp3=False
		if times([depth[-1]],param)[-1]>last_t(param[0]):
			tmp3=False
		return tmp3

	def last_t(fi):
		return(1./lam)*log(fi/(lam*det_lim))

	def ln_like_data(param):
		Asup=param[1:Ran+1]*density
		loglike = 0.
		tmp2=param[0]/lam
		ts=times(dep_time_data,param)
		for i in range(len(activity)):
			A_i= Asup[i] + tmp2 *(exp(-lam*ts[int(X0[i])] ) - exp(-lam*ts[int(X1[i])]) )
			Tau=.5*(sd_act[i]**(-2.))
			loglike = loglike + Tau*((A_i-activity[i])**2.)
		return loglike

	def ln_like_T(param):
		Asup=param[1:Ran+1]*density
		loglike = 0.
		tmp2=param[0]/lam
		ts=times(dep_time_data,param)
		for i in range(len(activity)):
			A_i= Asup[i] + tmp2 *(exp(-lam*ts[int(X0[i])] ) - exp(-lam*ts[int(X1[i])]) )#exp(tmp2 - lam*ts[int(X0[i])] ) - exp(tmp2 - lam*ts[int(X1[i])])
			Tau=.5*(sd_act[i]**(-2.))
			loglike = loglike + 3.5*log(4. + Tau*((A_i-activity[i])**2.) )
		return loglike


	def ln_like_supp(param):
		logsupp=0.
		for i in range(len(supp)):
			Tau=.5*(sd_supp[i]**-2.)
			logsupp = logsupp + Tau*((param[1+i]-supp[i])**2.)
		return logsupp

	def times(x,param):
		w=param[Ran+1]
		a=param[Ran+2:]
		t1=m-1
		ms1=array([a[m-1]])
		while t1 >0 :
		   ms1= append(ms1, w*ms1[-1]+(1-w)*a[t1-1])
		   t1 -= 1
		ms=ms1[::-1]
		ages=array([])
		y_last=append([0],array([sum(ms[:i+1]*by) for i in range(len(ms))] ) )
		for i in range(len(x)):
		   k1=0
		   while  breaks[k1]< x[i]:
			   k1 += 1
		   ages=append(ages,y_last[k1]+(ms[k1]*(by-(breaks[k1]-x[i]))))
		return ages

	def pendi(param):
		w=param[Ran+1]
		a=param[Ran+2:]
		t1=m-1
		ms1=array([a[m-1]])
		while t1 >0 :
		   ms1= append(ms1, w*ms1[-1]+(1-w)*a[t1-1])
		   t1 -= 1
		ms=ms1[::-1]
		return ms



	def ln_prior_supp(param):
		prior=0.
		prior= prior -  (  (shapefi-1.)*log(param[0])-(param[0]/scale_fi) )# prior for fi
		for k in range(Ran):
			prior= prior -  (  (shapeAS-1.)*log(param[1+k])-(param[1+k]/scale_As) )# prior for supp
		prior= prior -  ( ((1./by)-1.)*log(param[1+Ran])- log(by)+  ((1./by)*(shape1_m-1.))*log(param[1+Ran]) + (shape2_m - 1.)*log(1.-param[1+Ran]**(1./by) ) )# prior for w	#
		for ms in range(m):
		   prior= prior -  (  (shape_acc-1.)*log(param[ms+2+Ran])-(param[ms+2+Ran]/scale_acc) )
		return prior
		
	if Cs==True:
		def Cslike(param):		
			return 0.
	else:	
		def Cslike(param):
			tcs=times([Cs],param)
			Tau=.5*(1.**-2.)
			return Tau*(((Sdate-tcs)-CSTdate)**2.)
	


	if T_mod:
		log_data=ln_like_T
	else:
		log_data=ln_like_data

	def obj(param):
		objval= ln_like_supp(param) + ln_prior_supp(param) + log_data(param)+Cslike(param)
		return objval




	#################### Initial valules
	print("Seaching initial values")
	fi_ini_1= unif.rvs(size=1,loc=50, scale=200)  #200.
	fi_ini_2= unif.rvs(size=1,loc=250, scale=150) #100.
	supp_ini_1= unif.rvs(size=Ran,loc=15, scale=30) #5.
	supp_ini_2= unif.rvs(size=Ran,loc=1, scale=15) #20.
	w_ini = unif.rvs(size=1,loc=.2,scale=.3) #.3
	w_ini0 = unif.rvs(size=1,loc=.3,scale=.3)  #.7
	m_ini_1=unif.rvs(size=m,loc=0, scale=15)  #  repeat(array(3.1),m,axis=0)
	m_ini_2=unif.rvs(size=m,loc=0, scale=15)  # repeat(array(.5),m,axis=0)
	#print("here")
	x=append(append(append(fi_ini_1,supp_ini_1),w_ini), m_ini_1)
	xp=append(append(append(fi_ini_2,supp_ini_2),w_ini0), m_ini_2)

	while not support(x):
	 	m_ini_1=unif.rvs(size=m,loc=0, scale=3)
	 	x=append(append(append(fi_ini_1,supp_ini_1),w_ini), m_ini_1)


	while not support(xp):
	 	m_ini_2=unif.rvs(size=m,loc=0, scale=3)
	 	xp=append(append(append(fi_ini_2,supp_ini_2),w_ini0), m_ini_2)



	################### MCMC
	################## New MCMC test

	thi = int((len(x)))*thi #thi = 25, 50, 100
	burnin=len(xp) *burnin  #burin 10000 20000
	print("Total iterations,")
	print(burnin + iterations*thi)


	leadchrono = pytwalk(n=len(x),U=obj,Supp=support, ww=[ 0.0, 0.4918, 0.4918, 0.0082+0.082, 0.0])
	i, k ,k0, n=0 , 0, 0, len(x)
	U , Up = obj(x), obj(xp)
	por=int(iterations/10.)
	Output = zeros((iterations+1, n+1))
	Output[ 0, 0:n] = x.copy()
	Output[ 0, n] = U
	por2=int(burnin/5.)
	while i< iterations:
		onemove=leadchrono.onemove(x, U, xp, Up)
		k+= 1
		if (all([k<burnin,k % por2==0]) ):
			print("burn-in progress")
			print int(100*(k+.0)/burnin)
		if (unif.rvs() < onemove[3] ):
			x, xp, ke, A, U, Up =onemove
			k0+=1
			if all([k % thi ==0 , k>int(burnin)]):
				Output[i+1,0:n] = x.copy()
				Output[i+1,n] = U
				if any([i % por==0, i==0]) :
					print int(100*(i+.0)/iterations),"%"
				   #print((time.clock()-tiempomedir)/60)
				i+= 1
		else:
			if all([k % thi ==0 , k>int(burnin)]):
				Output[i+1,0:n] = x.copy()
				Output[i+1,n] = U
				if any([i % por==0, i==0]) :
					print int(100*(i+.0)/iterations),"%"
				   #print((time.clock()-tiempomedir)/60)
				i+= 1



	#Output=array(Output)
	print("Acceptance rate")
	print((100.0 * k0)/(burnin + iterations*thi))


	##################
	"""
	out 0 -> Fi
	out 1 -> Supported Activity
	out 2 -> w
	out 3-n -> dates
	out -1 -> Energy
	"""
	savetxt(dirt+'Results '+corename+'/Results_output.csv', Output,delimiter=',')
	estim=[]
	for i in range((iterations-1)):
		estim.append(times(breaks,Output[(i+1),:-1])  )
	estim=array(estim)
	savetxt(dirt+'Results '+corename+'/dates.csv', estim  )
	intervals=[]

	for i in range(len(estim[1,])):
		sort=sorted(estim[:,(i)])
		mean=sum(sort)/len(sort)
		disc=int(len(sort)*.025)+1
		disc1=int(len(sort)*.975)
		sort=sort[disc:disc1]
		intervals.append([breaks[i],sort[0],mean,sort[-1]])

	savetxt(dirt+'Results '+corename+'/intervals.csv', intervals ,delimiter=',')
	depths=array([append([0.0],breaks)])
	savetxt(dirt+'Results '+corename+'/depths.csv', depths,delimiter=',')

	grafdepts=linspace(0,breaks[-1],resolution)
	grafdepts2=grafdepts+(grafdepts[1]-grafdepts[0])/2
	grafdepts2=grafdepts2[0:(len(grafdepts2)-1)]

	grafage=linspace(0,(max(estim[:,-1])+.10),resolution)
	y=[]
	for i in range(len(depths[0,:])-1):
		logvect=array(grafdepts2>depths[0,i])*array(grafdepts2<=depths[0,i+1])
		for k in range(len(logvect)):
			if logvect[k]:
				if i!=0:
					y1=estim[:,i-1]+((estim[:,i]-estim[:,i-1])/(depths[0,i+1]-depths[0,i]))*(grafdepts[k]-depths[0,i])
					porc=[]
					for posi in range(len(grafage)-1):
						porc.append(sum(array(y1>=grafage[posi])*array(y1<grafage[posi+1]) ))
					y.append(porc/(max(porc)+0.0 ))
				else:
					y1=((estim[:,i]/depths[0,i+1])*(grafdepts[k]) )
					porc=[]
					for posi in range(len(grafage)-1):
						porc.append(sum(array(y1>=grafage[posi])*array(y1<grafage[posi+1]) ))
					y.append(porc/(max(porc)+0.0 ))

	savetxt(dirt+'Results '+corename+'/Graphs.csv', array(y),delimiter=',')
	slopes=[]
	for i in range(iterations-1):
		slopes.append(pendi(Output[(i+1),:-1])  )
	savetxt(dirt+'Results '+corename+'/Slopes.csv', array(slopes),delimiter=',')






#D="/home/aquinom/Documents/Joan -Caro/SAI2016/"
#plomo="SAI16Ra.csv"

#plumMCMC(D,plomo,False,0,.01,100,2.5,2.,.7,1.5,10,1,1,1,1,10,123)






########################################################
###  This is the python implementation of the t-walk ### 
###  By Andres Christen.                             ###
###  A generic self adjusting MCMC                   ###
###  see:  http://www.cimat.mx/~jac/twalk/           ###
###  see also twalktutorial.py                       ###
########################################################



from numpy.random import uniform, normal
from numpy import ones, zeros, cumsum, shape, mat, cov, mean, ceil, matrix, sqrt
from numpy import floor, exp, log, sum, pi, savetxt, loadtxt

from time import time, localtime, strftime


try:
	from pylab import plot, hist, xlabel, ylabel, title
except:
	print "pytwalk: WARNING: pylab module not available, Ana, TS and Hist methods will fail."

#### Some auxiliar functions and constants
## square of the norm.
def SqrNorm(x):
	return sum(x*x) 

log2pi = log(2*pi)
log3 = log(3.0)

def Remain( Tr, it, sec1, sec2):
	""" Remaining time Information messages:
        total iterations Tr, current it, start time, current time, as returned by time() (floats)."""

	# how many seconds remaining
	ax = int( (Tr - it) *  ((sec2 - sec1)/it) )


	if (ax < 1):

		return " "

	if (ax < 60):

		return "Finish in approx. %d seconds." % (ax,)

	if (ax <= 360):

		return "Finish in approx. %d minutes and %d seconds." % ( ax // 60, ax % 60)

	if (ax > 360):

		ax += sec2  # current time plus seconds remaining=end time
		return "Finish by " + strftime("%a, %d %b %Y, %H:%M.", localtime(ax))






class pytwalk:
	"""This is the t-walk class.

	Initiates defining the dimension= n and -log of the objective function= U,
	Supp defines the support, returns True if x within the support, eg:

	Mytwalk = pytwalk( n=3, U=MyMinusLogf, Supp=MySupportFunction).

	Then do: Mytwalk.Run?

	Other parameter are:
	ww= the prob. of choosing each kernel, aw, at, n1phi (see inside twalk.py)
	with default values as in the paper, normally NOT needed to be changed."""

	def __init__( self, n, U=(lambda x: sum(0.5*x**2)), Supp=(lambda x: True),
		ww=[0.0000, 0.4918, 0.4918, 0.0082, 0.0082], aw=1.5, at=6.0, n1phi=4.0):

		self.n = n
		self.U = U
		self.Supp = Supp
		self.Output = zeros((1, n+1)) ### No data (MCMC output) yet
		self.T = 1
		self.Acc = zeros(6)  ### To save the acceptance rates of each kernel, and the global acc. rate

		#### Kernel probabilities
		self.Fw = cumsum(ww)
		
		#### Parameters for the propolsals
		self.aw = aw  ### For the walk move
		self.at = at ### For the Traverse move

		#n1phi = 5 ### expected value of parameters to move
		self.pphi = min( n, n1phi)/(1.0*n) ### Prob. of choosing each par.
		
		self.WAIT = 30


	def _SetUpInitialValues( self, x0, xp0):
		"""Private method."""

		### Check x0 and xp0 in the support

		if any(abs(x0 -xp0) <= 0):
			print "pytwalk: ERROR, not all entries of initial values different."
			return [ False, 0.0, 0.0]

		if not(self.Supp(x0)):
			print "pytwalk: ERROR, initial point x0 out of support."
			return [ False, 0.0, 0.0]
		u = self.U(x0)

		if not(self.Supp(xp0)):
			print "pytwalk: ERROR, initial point xp0 out of support."
			return [ False, u, 0.0]
		up = self.U(xp0)
		
		return [ True, u, up]



	def Run( self, T, x0, xp0):
		"""Run the twalk.
		
		   Run( T, x0, xp0),
		   T = Number of iterations.
		   x0, xp0, two initial points within the support,
		   ***each entry of x0 and xp0 most be different***. 
		"""
		
		sec = time()
		print "pytwalk: Running the twalk with %d iterations." % (T,), strftime("%a, %d %b %Y, %H:%M.", localtime(sec))

		### Check x0 and xp0 are in the support
		[ rt, u, up] = self._SetUpInitialValues( x0, xp0)

		if (not(rt)):
			return 0
		

		### send an estimation for the duration of the sampling if 
		### evaluating the ob. func. twice (in self._SetUpInitialValues) takes more than one second

		sec2 = time() # last time we sent a message
		print "       " + Remain( T, 2, sec, sec2)

		x = x0     ### Use x and xp by reference, so we can retrive the last values used
		xp = xp0

		### Set the array to place the iterations and the U's ... we donot save up's
		self.Output = zeros((T+1, self.n+1))
		self.T = T+1
		self.Acc = zeros(6)
		kercall = zeros(6) ## Times each kernel is called
				
		#### Make local references for less writing
		n = self.n
		Output = self.Output
		U = self.U
		Supp = self.Supp
		Acc = self.Acc
		Fw = self.Fw
		
		Output[ 0, 0:n] = x.copy()
		Output[ 0, n] = u

		j1=1
		j=0

		### Sampling
		for it in range(T):
		
			y, yp, ke, A, u_prop, up_prop = self.onemove( x, u, xp, up)

			kercall[ke] += 1
			kercall[5] += 1 
			if (uniform() < A):  
				x = y.copy()   ### Accept the propolsal y
				u = u_prop
				xp = yp.copy()   ### Accept the propolsal yp
				up = up_prop
				
				Acc[ke] += 1
				Acc[5] += 1


			### To retrive the current values
			self.x = x
			self.xp = xp
			self.u = u
			self.up = up

			Output[it+1,0:n] = x.copy()
			Output[it+1,n] = u

			### Estimate the remaing time, every 2**j1 iterations
			if ((it % (1 << j1)) == 0):

				j1 += 1
				j1 = min( j1, 10)  # check the time at least every 2^10=1024 iterations
				ax = time()
				if ((ax - sec2) > (1 << j)*self.WAIT): # Print an estimation every WAIT*2**j 

					print "pytwalk: %10d iterations so far. " % (it,) + Remain( T, it, sec, ax)
					sec2 = ax
					j += 1
					j1 -= 1 # check the time as often 



		
		if (Acc[5] == 0):
			print "pytwalk: WARNING,  all propolsals were rejected!"
			print strftime("%a, %d %b %Y, %H:%M:%S.", localtime(time()))
			return 0
		else:
			print "pytwalk: finished, " + strftime("%a, %d %b %Y, %H:%M:%S.", localtime(time()))
			

		for i in range(6):
			if kercall[i] != 0:
				Acc[i] /= kercall[i]
		return 1


	def  onemove( self, x, u, xp, up):
		"""One move of the twalk.  This is basically the raw twalk kernel.
		   It is usefull if the twalk is needed inside a more complex MCMC.
		
		   onemove(x, u, xp, up),
		   x, xp, two points WITHIN the support ***each entry of x0 and xp0 must be different***.
		   and the value of the objective at x, and xp
		   u=U(x), up=U(xp).
		   
		   It returns: [y, yp, ke, A, u_prop, up_prop]
		   y, yp: the proposed jump
		   ke: The kernel used, 0=nothing, 1=Walk, 2=Traverse, 3=Blow, 4=Hop
		   A: the M-H ratio
		   u_prop, up_prop: The values for the objective func. at the proposed jumps 
		"""

		#### Make local references for less writing
		n = self.n
		U = self.U
		Supp = self.Supp
		Fw = self.Fw
		
		ker = uniform() ### To choose the kernel to be used
		ke = 1
		A = 0
		
		## Kernel nothing exchange x with xp, not used
		if ((0.0 <= ker) & (ker < Fw[0])): 
			ke = 0
			y = xp.copy()
			up_prop = u
			yp = x.copy()
			u_prop = up
			### A is the MH acceptance ratio
			A = 1.0;  #always accepted


		## The Walk move
		if ((Fw[0] <= ker) & (ker < Fw[1])):
			
			ke = 1

			dir = uniform()

			if ((0 <= dir) & (dir < 0.5)):  ## x as pivot
		
				yp = self.SimWalk( xp, x)

				y = x.copy()
				u_prop = u

				if ((Supp(yp)) & (all(abs(yp - y) > 0))):
					up_prop = U(yp)
					A = exp(up - up_prop)
				else:
					up_prop = None
					A = 0; ##out of support, not accepted
						
			else:  ## xp as pivot

				y = self.SimWalk( x, xp)

				yp = xp.copy()
				up_prop = up

				if ((Supp(y)) & (all(abs(yp - y) > 0))):
					u_prop = U(y)
					A = exp(u - u_prop)
				else:
					u_prop = None
					A = 0; ##out of support, not accepted


		#### The Traverse move
		if ((Fw[1] <= ker) & (ker < Fw[2])):

			ke = 2
			dir = uniform()

			if ((0 <= dir) & (dir < 0.5)):  ## x as pivot

				beta = self.Simbeta()
				yp = self.SimTraverse( xp, x, beta)

				y = x.copy()
				u_prop = u
				
				if Supp(yp):				
					up_prop = U(yp)
					if (self.nphi == 0):
						A = 1 ###Nothing moved
					else:
						A = exp((up - up_prop) +  (self.nphi-2)*log(beta))
				else:
					up_prop = None
					A = 0 ##out of support, not accepted
			else:			## xp as pivot

				beta = self.Simbeta()
				y = self.SimTraverse( x, xp, beta)

				yp = xp.copy()
				up_prop = up

				if Supp(y):
					u_prop = U(y)
					if (self.nphi == 0):
						A = 1 ###Nothing moved
					else:
						A = exp((u - u_prop) +  (self.nphi-2)*log(beta))
				else:
					u_prop = None
					A = 0 ##out of support, not accepted

		### The Blow move
		if ((Fw[2] <= ker) & (ker < Fw[3])): 

			ke = 3
			dir = uniform()

			if ((0 <= dir) & (dir < 0.5)):  ## x as pivot
				yp = self.SimBlow( xp, x)
				
				y = x.copy()
				u_prop = u
				if ((Supp(yp)) & all(yp != x)):
					up_prop = U(yp)
					W1 = self.GBlowU( yp, xp,  x)
					W2 = self.GBlowU( xp, yp,  x) 
					A = exp((up - up_prop) + (W1 - W2))
				else:
					up_prop = None
					A = 0 ##out of support, not accepted
			else:  ## xp as pivot
				y = self.SimBlow( x, xp)

				yp = xp.copy()
				up_prop = up
				if ((Supp(y)) & all(y != xp)):
					u_prop = U(y)
					W1 = self.GBlowU(  y,  x, xp)
					W2 = self.GBlowU(  x,  y, xp)
					A = exp((u - u_prop) + (W1 - W2))
				else:
					u_prop = None
					A = 0 ##out of support, not accepted
		

		### The Hop move
		if ((Fw[3] <= ker) & (ker < Fw[4])): 

			ke = 4
			dir = uniform()

			if ((0 <= dir) & (dir < 0.5)):  ## x as pivot
				yp = self.SimHop( xp, x)
				
				y = x.copy()
				u_prop = u
				if ((Supp(yp)) & all(yp != x)):
					up_prop = U(yp)
					W1 = self.GHopU( yp, xp,  x)
					W2 = self.GHopU( xp, yp,  x) 
					A = exp((up - up_prop) + (W1 - W2))
				else:
					up_prop = None
					A = 0 ##out of support, not accepted
			else:  ## xp as pivot
				y = self.SimHop( x, xp)

				yp = xp.copy()
				up_prop = up
				if ((Supp(y)) & all(y != xp)):
					u_prop = U(y)
					W1 = self.GHopU(  y,  x, xp)
					W2 = self.GHopU(  x,  y, xp)
					A = exp((u - u_prop) + (W1 - W2))
				else:
					u_prop = None
					A = 0 ##out of support, not accepted
		
		return [y, yp, ke, A, u_prop, up_prop]



#################################################################################
##### Auxiliars for the kernels

	### Used by the Walk kernel
	def SimWalk( self, x, xp):
		aw = self.aw
		n = self.n
		
		phi = (uniform(size=n) < self.pphi) ### parametrs to move
		self.nphi = sum(phi)
		z = zeros(n)

		for i in range(n):
			if phi[i]:
				u = uniform()
				z[i] = (aw/(1+aw))*(aw*u**2.0 + 2.0*u - 1.0)

		return x + (x - xp)*z

	#### Used by the Traverse kernel
	def Simbeta(self):
		at = self.at
		if (uniform() < (at-1.0)/(2.0*at)):
			return exp(1.0/(at+1.0)*log(uniform()))
		else:
			return exp(1.0/(1.0-at)*log(uniform()))

	def SimTraverse( self,  x, xp, beta):
		n = self.n
	
		phi = (uniform(size=n) < self.pphi)
		self.nphi = sum(phi)

		rt = x.copy()
		for i in range(n):
			if (phi[i]):
				rt[i] = xp[i] + beta*(xp[i] - x[i])
			
		return rt


	### Used by the Blow kernel
	def SimBlow( self, x, xp):
		n = self.n
	
		self.phi = (uniform(size=n) < self.pphi)
		self.nphi = sum(self.phi)
	
		self.sigma = max(self.phi*abs(xp - x))

		rt = x.copy()
		for i in range(n):
			if (self.phi[i]):
				rt[i] = xp[i] + self.sigma * normal()
			
		return rt


	def GBlowU( self, h, x, xp):
		nphi = self.nphi
		self.sigma = max(self.phi*abs(xp - x)) #recalculate sigma, but same phi	
		if (nphi > 0):
			return (nphi/2.0)*log2pi + nphi*log(self.sigma) + 0.5*SqrNorm(h - xp)/(self.sigma**2)
		else:
			return 0


	### Used by the Hop kernel
	def SimHop( self, x, xp):
		n = self.n
	
		self.phi = (uniform(size=n) < self.pphi)
		self.nphi = sum(self.phi)
	
		self.sigma = max(self.phi*abs(xp - x))/3.0

		rt = x.copy()
		for i in range(n):
			if (self.phi[i]): 
				rt[i] = x[i] + self.sigma * normal()

		return rt


	def GHopU( self, h, x, xp): ## It is actually equal to GBlowU!
		nphi = self.nphi
		self.sigma = max(self.phi*abs(xp - x))/3.0 ##Recalculate sigma, but same phi

		if (nphi > 0):
			return (nphi/2.0)*log2pi  + nphi*log(self.sigma) + 0.5*SqrNorm(h - xp)/(self.sigma**2)
		else:
			return 0



#################################################################################
#####  Output analysis auxiliar methods

	def IAT( self, par=-1, start=0, end=0, maxlag=0):
		"""Calculate the Integrated Autocorrelation Times of parameters par
		   the default value par=-1 is for the IAT of the U's"""
		if (end == 0):
			end = self.T

		if (self.Acc[5] == 0):
			print "twalk: IAT: WARNING,  all propolsals were rejected!"
			print "twalk: IAT: Cannot calculate IAT, fixing it to the sample size."
			return self.T

		iat = IAT( self.Output, cols=par, maxlag=maxlag, start=start, end=end)
		
		return iat
	

	def TS( self, par=-1, start=0, end=0):
		"""Plot time sries of parameter par (defualt = log f) etc."""
		if par == -1:
			par = self.n
		
		if (end == 0):
			end = self.T

		if (par == self.n):
			plot( range( start, end), -1*self.Output[ start:end, par])
			ylabel("Log of Objective")
		else:
			plot( range( start, end), self.Output[ start:end, par])
			ylabel("Parameter %d" % par)
		xlabel("Iteration")


	def Ana( self, par=-1, start=0, end=0):
		"""Output Analysis, TS plots, accepatnce rates, IAT etc."""
		if par == -1:
			par = self.n

		if (end == 0):
			end = self.T

		print "Acceptance rates for the Walk, Traverse, Blow and Hop kernels:" + str(self.Acc[1:5])
		print "Global acceptance rate: %7.5f" % self.Acc[5]
		
		iat = self.IAT( par=par, start=start, end=end)
		print "Integrated Autocorrelation Time: %7.1f, IAT/n: %7.1f" % (iat, iat/self.n)
		
		self.TS( par=par, start=start, end=end)


	def Hist( self, par=-1, start=0, end=0, g=(lambda x: x[0]), xlab=None, bins=20, normed=False):
		"""Basic histigrams and output analysis.  If par=-1, use g.
		   The function g provides a transformation to be applied to the data,
		   eg g=(lambda x: abs(x[0]-x[1]) would plot a histogram of the distance
		   between parameters 0 and 1, etc."""

		if (end == 0):
			end = self.T

		if (par == -1):
			ser = zeros(end-start)
			for it in range(end-start):
				ser[it] = g(self.Output[ it+start, :])
			if (xlab == None):
				xlab = "g"
		else:
			ser = self.Output[ start:end, par]
			if (xlab == None):
				xlab = "parameter %d" % (par,)
			
		xlabel(xlab)
		print "Mean for %s= %f" % ( xlab, mean(ser))
		return hist( ser, bins=bins, normed=normed)



	def Save( self, fnam, start=0, thin=1):
		"""Saves the Output as a text file, starting at start (burn in), with thinning (thin)."""

		print "Saving output, all pars. plus the U's in file", fnam
		
		savetxt( fnam, self.Output[ start::thin, ])



	def Load( self, fnam, start=0, thin=1):
		"""Loads the Output from a text file, typically written with the Save method.
		   It will overwrite any other twalk output.  Updates the dimension n and the sample size T."""
		
		print "Loading output from file", fnam
		
		self.Output = loadtxt(fnam)
		self.T, self.n = self.Output.shape
		self.n -= 1

		
##### A simple Random Walk M-H
	def RunRWMH( self, T, x0, sigma):
		"""Run a simple Random Walk M-H"""

		print "pytwalk: This is the Random Walk M-H running with %d iterations." % T
		### Local variables
		x = x0.copy()
		if not(self.Supp(x)):
			print("pytwalk: ERROR, initial point x0 out of support.")
			return 0

		u = self.U(x)
		n = self.n

		### Set the array to place the iterations and the U's
		self.Output = zeros((T+1, n+1))
		self.Acc = zeros(6)
				
		#### Make local references for less writing
		Output = self.Output
		U = self.U
		Supp = self.Supp
		Acc = self.Acc
		
		Output[ 0, 0:n] = x.copy()
		Output[ 0, n] = u

		y = x.copy()
		for it in range(T):
			y = x + normal(size=n)*sigma ### each entry with sigma[i] variance 
			if Supp(y):        ### If it is within the support of the objective
				uprop = U(y)   ### Evaluate the objective
				if (uniform() < exp(u-uprop)):  
					x = y.copy()   ### Accept the propolsal y
					u = uprop
					Acc[5] += 1

			Output[it+1,0:n] = x
			Output[it+1,n] = u
		
		if (Acc[5] == 0):
			print("pytwalk: WARNING,  all propolsals were rejected!")
			return 0

		Acc[5] /= T;
		return 1








############################################################################################
#### Auxiliary functions to calculate Integrated autocorrelation times of a time series 


####  Calculates an autocovariance 2x2 matrix at lag l in column c of matrix Ser with T rows
####  The variances of each series are in the diagonal and the (auto)covariance in the off diag.
def AutoCov( Ser, c, la, T=0):
	if (T == 0):
		T = shape(Ser)[0]  ### Number of rows in the matrix (sample size)

	return cov( Ser[0:(T-1-la), c], Ser[la:(T-1), c], bias=1)
	
	
	

#### Calculates the autocorrelation from lag 0 to lag la of columns cols (list)
#### for matrix Ser
def AutoCorr( Ser, cols=0, la=1):
	T = shape(Ser)[0]  ### Number of rows in the matrix (sample size)

	ncols = shape(mat(cols))[1] ## Number of columns to analyse (parameters)

	#if ncols == 1:
	#	cols = [cols]
		
	### Matrix to hold output
	Out = matrix(ones((la+1)*ncols)).reshape( la+1, ncols)
		
	for c in range(ncols):
		for l in range( 1, la+1):  
			Co = AutoCov( Ser, cols[c], l, T) 
			Out[l,c] = Co[0,1]/(sqrt(Co[0,0]*Co[1,1]))
	
	return Out
	

### Makes an upper band matrix of ones, to add the autocorrelation matrix
### gamma = auto[2*m+1,c]+auto[2*m+2,c] etc. 
### MakeSumMat(lag) * AutoCorr( Ser, cols=c, la=lag) to make the gamma matrix
def MakeSumMat(lag):
	rows = (lag)/2   ### Integer division!
	Out = mat(zeros([rows,lag], dtype=int))
	
	for i in range(rows): 
		Out[i,2*i] = 1
		Out[i,2*i+1] = 1
	
	return Out


### Finds the cutting time, when the gammas become negative
def Cutts(Gamma):
	cols = shape(Gamma)[1]
	rows = shape(Gamma)[0]
	Out = mat(zeros([1,cols], dtype=int))
	Stop = mat(zeros([1,cols], dtype=bool))
	
	if (rows == 1):
		return Out
		
	i = 0
	###while (not(all(Stop)) & (i < (rows-1))):
	for i in range(rows-1):
		for j in range(cols):  ### while Gamma stays positive and decreasing
			if (((Gamma[i+1,j] > 0.0) & (Gamma[i+1,j] < Gamma[i,j])) & (not Stop[0,j])):
				Out[0,j] = i+1 ## the cutting time for colomn j is i+i
			else:
				Stop[0,j] = True
		i += 1
	
	
	return Out


####  Automatically find a maxlag for IAT calculations
def AutoMaxlag( Ser, c, rholimit=0.05, maxmaxlag=20000):
	Co = AutoCov( Ser, c, la=1)
	rho = Co[0,1]/Co[0,0]  ### lag one autocorrelation
	
	### if autocorrelation is like exp(- lag/lam) then, for lag = 1
	lam = -1.0/log(abs(rho)) 
	
	### Our initial guess for maxlag is 1.5 times lam (eg. three times the mean life)
	maxlag = int(floor(3.0*lam))+1
	
	### We take 1% of lam to jump forward and look for the
	### rholimit threshold
	jmp = int(ceil(0.01*lam)) + 1
	
	T = shape(Ser)[0]  ### Number of rows in the matrix (sample size)

	while ((abs(rho) > rholimit) & (maxlag < min(T/2,maxmaxlag))):
		Co = AutoCov( Ser, c, la=maxlag)
		rho = Co[0,1]/Co[0,0]
		maxlag = maxlag + jmp
		###print("maxlag=", maxlag, "rho", abs(rho), "\n")
		
	maxlag = int(floor(1.3*maxlag));  #30% more
	
	if (maxlag >= min(T/2,maxmaxlag)): ###not enough data
		fixmaxlag = min(min( T/2, maxlag), maxmaxlag)
		print "AutoMaxlag: Warning: maxlag= %d > min(T/2,maxmaxlag=%d), fixing it to %d" % (maxlag, maxmaxlag, fixmaxlag)
		return fixmaxlag
	
	if (maxlag <= 1):
		fixmaxlag = 10
		print "AutoMaxlag: Warning: maxlag= %d ?!, fixing it to %d" % (maxlag, fixmaxlag)
		return fixmaxlag
		
	print "AutoMaxlag: maxlag= %d." % maxlag
	return maxlag
	
	
### Find the IAT
def IAT( Ser, cols=-1,  maxlag=0, start=0, end=0):

	ncols = shape(mat(cols))[1] ## Number of columns to analyse (parameters)
	if ncols == 1:
		if (cols == -1):
			cols = shape(Ser)[1]-1 ### default = last column
		cols = [cols]
	
	if (end == 0):
		end = shape(Ser)[0]

	if (maxlag == 0):
		for c in cols:
			maxlag = max(maxlag, AutoMaxlag( Ser[start:end,:], c))

	#print("IAT: Maxlag=", maxlag)

	#Ga = MakeSumMat(maxlag) * AutoCorr( Ser[start:end,:], cols=cols, la=maxlag)
	
	Ga = mat(zeros((maxlag/2,ncols)))
	auto = AutoCorr( Ser[start:end,:], cols=cols, la=maxlag)
	
	### Instead of producing the maxlag/2 X maxlag MakeSumMat matrix, we calculate the gammas like this
	for c in range(ncols):
		for i in range(maxlag/2):
			Ga[i,c] = auto[2*i,c]+auto[2*i+1,c]
	
	cut = Cutts(Ga)
	nrows = shape(Ga)[0]
		
	ncols = shape(cut)[1]
	Out = -1.0*mat(ones( [1,ncols] ))
	
	if any((cut+1) == nrows):
		print("IAT: Warning: Not enough lag to calculate IAT")
	
	for c in range(ncols):
		for i in range(cut[0,c]+1):
			Out[0,c] += 2*Ga[i,c]
	
	return Out



############################################################################################












