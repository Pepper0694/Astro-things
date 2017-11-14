import theano
import theano.tensor as T
import numpy as np
from lasagne.layers import *
import lasagne
import pickle
import matplotlib.pyplot as plt

X = []
Y = []
Tg = []

a,b,c = pickle.load(open('sirius.pkl','rb'))

X.append(np.array(a).reshape((1,1230)))
Y.append(np.array(b).flatten())
Tg.append(np.array(c).flatten())

a,b,c = pickle.load(open('vegaIR.pkl','rb'))

plt.plot(a,b)
plt.show()

X.append(np.array(a).reshape((1,1230)))
Y.append(np.array(b).flatten())
Tg.append(np.array(c).flatten())









l_in = InputLayer((1,1230))
dense = DenseLayer(incoming=l_in, num_units=1230,nonlinearity=lasagne.nonlinearities.tanh)
d2 = DenseLayer(incoming=dense, num_units=1230,nonlinearity=lasagne.nonlinearities.tanh)
out = lasagne.layers.BiasLayer(incoming=dense)


f = DenseLayer(incoming=l_in, num_units=1230,nonlinearity=lasagne.nonlinearities.tanh)


Sirius = T.ivector('y')
observed = T.ivector('ssdf')



adder = lasagne.layers.get_output(out).flatten()
mult = lasagne.layers.get_output(f).flatten()


all_params = lasagne.layers.get_all_params(out)

get_s = theano.function([l_in.input_var],[adder,mult])

prediction = mult*observed + adder

get_out = theano.function([l_in.input_var,observed],prediction)

cost = T.mean(lasagne.objectives.squared_error(prediction,Sirius))


updates = lasagne.updates.sgd(cost, all_params,learning_rate=1)


train = theano.function(
    [l_in.input_var, observed, Sirius],
    cost, updates=updates,allow_input_downcast=True)



get_cost = theano.function([l_in.input_var, observed, Sirius], cost,allow_input_downcast=True)




try:
	while True:
		mc = 0
		for x,y,t in zip(X,Y,Tg):
			mc += 0.5*train(X[0],y,(1e10)*t)

		print(mc)
except:
	values = get_s(X[0])
	pickle.dump(values,open('values.p','wb'))
	for x,y,t in zip(X,Y,Tg):

		pred = values[0].flatten() + values[1].flatten()*y
		plt.plot(x.flatten(),y)
		plt.show()
		plt.clf()
		plt.plot(x.flatten(),pred)
		plt.show()
		plt.clf()



