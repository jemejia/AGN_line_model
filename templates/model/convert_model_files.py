import numpy as np
import sys as s


c=29979245800
wl_min=1200
wl_max=10000
model_in=np.loadtxt(s.argv[1])
wl_in=model_in[:,0]
L_in=model_in[:,2]
#wl_in=wl_in[::-1]
#L_in=L_in[::-1]

#amin=np.argmin( np.abs(wl_in - wl_min) )
#amax=np.argmin( np.abs(wl_in - wl_max) )


wl_out=wl_in
#wl_out=wl_in[amin:amax +1]
L_out=1e16*L_in*c/(wl_out*wl_out)

model_out=np.empty(  ( len(wl_out),2 )  )

model_out[:,0]=wl_out
model_out[:,1]=L_out

np.savetxt("cont_model_"+ s.argv[2] + ".txt", model_out)
