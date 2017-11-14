import pickle
import theano
import theano.tensor as T
import numpy as np
import lasagne
import matplotlib.pyplot as plt
from sklearn.svm import SVR
import random


#Find the max ratio of UnC to C so we can bring the values in a similar region to assist convergence
W, UnC,C = pickle.load(open('sirius.p','rb'))
m = UnC.max()/C.max()
W, UnC,C = pickle.load(open('vega.p','rb'))
m = max(m,UnC.max()/C.max())




X = []
y = []




#Load the values and put them in the training and validation arrays
#I multiply by 1e9 to bring the values between 0-100
W, UnC,C = pickle.load(open('sirius.p','rb'))

for a,b,c in zip(W,UnC,C):
	X.append(np.array([a,(b/m)*1e9]))
	y.append(c*1e9)

W, UnC,C = pickle.load(open('vega.p','rb'))


for a,b,c in zip(W,UnC,C):
	X.append(np.array([a,(b/m)*1e9]))
	y.append(c*1e9)





#Reshape to corrent format
#Each input has to be a vector
X = np.array(X).reshape((1409,1,2))
y = np.array(y).reshape((1409,1))


#This is the neural network

lin = lasagne.layers.InputLayer((1,2))
d1 = lasagne.layers.DenseLayer(lin,num_units=100,nonlinearity=lasagne.nonlinearities.tanh)
d2 = lasagne.layers.DenseLayer(d1,num_units=100)
out = lasagne.layers.DenseLayer(d2,num_units=1)



# Uncomment this to load the values after trained
# values = pickle.load(open('values.pkl','rb'))
# lasagne.layers.set_all_param_values(out, values)



target = T.fvector('i')

#Create the training functions

prediction = lasagne.layers.get_output(out).flatten()
get_out = theano.function([lin.input_var],prediction,allow_input_downcast=True)
cost = ((prediction-target)**2).sum()
params = lasagne.layers.get_all_params(out)
updates = lasagne.updates.sgd(cost,params,1e-2)

train = theano.function([lin.input_var,target],cost,updates=updates,allow_input_downcast=True)



#train loop, and randomize order each time

for q in range(10):
	c = list(zip(X,y))
	random.shuffle(c)
	X,y = zip(*c)
	cost = 0
	q = []
	for x,Y in zip(X,y):
		cost += train(x,Y)/len(y)
		q.append(get_out(x))
	print(cost)
	




values = lasagne.layers.get_all_param_values(out)
pickle.dump(values,open('values.pkl','wb'))






