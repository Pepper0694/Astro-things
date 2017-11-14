import theano
import theano.tensor as T
import numpy as np
from lasagne.layers import *
import lasagne
import pickle
import matplotlib.pyplot as plt
from lasagne.layers import BiasLayer


X,Y,A = pickle.load(open('sirius.p','rb'))

A = np.array(A)


X = X.reshape((1,1230))










l_in = InputLayer((1,1230))
dense = DenseLayer(incoming=l_in, num_units=1230,nonlinearity=lasagne.nonlinearities.tanh)
d2 = dense = DenseLayer(incoming=dense, num_units=1230,nonlinearity=lasagne.nonlinearities.tanh)
out = lasagne.layers.BiasLayer(incoming=dense)

Sirius = T.ivector('y')
observed = T.ivector('ssdf')



prediction = lasagne.layers.get_output(out).flatten()


get_out = theano.function([l_in.input_var],prediction)


all_params = lasagne.layers.get_all_params(out)


cost = T.sum(lasagne.objectives.squared_error(observed+prediction,Sirius))


updates = lasagne.updates.momentum(cost, all_params,learning_rate=0.1)


train = theano.function(
    [l_in.input_var, observed, Sirius],
    cost, updates=updates,allow_input_downcast=True)



get_cost = theano.function([l_in.input_var, observed, Sirius], cost,allow_input_downcast=True)
print('start training')


#try:
#	while True:
#		train(X,Y,A)
#		print(get_cost(X,Y,A))

#except:
#	values = lasagne.layers.get_all_param_values(out)
#	pickle.dump(values,open('Regression.pkl','wb'))





